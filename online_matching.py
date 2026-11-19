"""Random-arrival online fractional hedonic matching, with exact offline references."""
import argparse
from collections import defaultdict, deque
from functools import lru_cache
import gzip
from hashlib import sha256
from itertools import permutations
import json
from math import lcm
from pathlib import Path
from random import Random
from statistics import mean, quantiles

PROTOCOL_HASH = "26ac16bdfea94fa13bfb3eef1ef651ed97c1c4e493ae697814ed01a26e101cef"
N, BATCH, ORDERS, WAIT = 9, 48, 32, 3
POLICIES = ("fractional", "integral", "greedy_matching", "greedy_coalition")
SCALE = lcm(*range(1, N + 1))


def key_seed(*keys):
    return int.from_bytes(sha256(repr(keys).encode()).digest()[:8], "big")


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def assignment(matrix):
    """Exact integer Hungarian method for a square maximum-weight assignment."""
    size = len(matrix)
    potential_row = [0] * (size + 1)
    potential_col = [0] * (size + 1)
    owner = [0] * (size + 1)
    previous = [0] * (size + 1)
    for row in range(1, size + 1):
        owner[0] = row
        column = 0
        visited = set()
        best = [None] * (size + 1)
        while True:
            visited.add(column)
            current = owner[column]
            delta, next_column = None, None
            for dest in range(1, size + 1):
                if dest in visited:
                    continue
                cost = -matrix[current - 1][dest - 1] - potential_row[current] - potential_col[dest]
                if best[dest] is None or cost < best[dest]:
                    best[dest], previous[dest] = cost, column
                if delta is None or best[dest] < delta:
                    delta, next_column = best[dest], dest
            for dest in range(size + 1):
                if dest in visited:
                    potential_row[owner[dest]] += delta
                    potential_col[dest] -= delta
                elif best[dest] is not None:
                    best[dest] -= delta
            column = next_column
            if owner[column] == 0:
                break
        while column:
            parent = previous[column]
            owner[column] = owner[parent]
            column = parent
    solution = [0] * size
    for dest in range(1, size + 1):
        solution[owner[dest] - 1] = dest - 1
    return solution


def integral_solver(weights):
    @lru_cache(maxsize=None)
    def best(mask):
        if not mask:
            return 0, ()
        first = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << first)
        value, edges = best(rest)
        for second in bits(rest):
            tail, subset = best(rest ^ (1 << second))
            proposal = tail + weights[first][second]
            if proposal > value:
                value, edges = proposal, ((first, second),) + subset
        return value, edges
    return best


def coalition_welfares(weights):
    welfare = [0] * (1 << N)
    sums = [0] * (1 << N)
    for mask in range(1, 1 << N):
        first = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << first)
        sums[mask] = sums[rest] + sum(weights[first][i] for i in bits(rest))
        welfare[mask] = 2 * sums[mask] * (SCALE // mask.bit_count())
    return welfare


def offline_coalitions(welfare):
    @lru_cache(maxsize=None)
    def best(mask):
        if not mask:
            return 0
        first = mask & -mask
        remain = mask ^ first
        subset, optimum = remain, 0
        while True:
            block = subset | first
            optimum = max(optimum, welfare[block] + best(mask ^ block))
            if not subset:
                return optimum
            subset = (subset - 1) & remain
    return best


def prefix_solver(weights, tie_seed, fraction):
    @lru_cache(maxsize=None)
    def choices(mask):
        vertices = list(bits(mask))
        rng = Random(key_seed("ties", tie_seed, mask))
        perturbed = [[0] * len(vertices) for _ in vertices]
        for i, src in enumerate(vertices):
            for j in range(i, len(vertices)):
                dest = vertices[j]
                mark = rng.randrange(10000)
                # Summed tie noise is less than 100,000 in either objective.
                value = 100000 * weights[src][dest] + mark
                perturbed[i][j] = perturbed[j][i] = value
        if fraction:
            outgoing = assignment(perturbed)
            incoming = [0] * len(vertices)
            for i, j in enumerate(outgoing):
                incoming[j] = i
            return {v: (vertices[outgoing[i]], vertices[incoming[i]])
                    for i, v in enumerate(vertices)}
        # Paper's integral ablation, with the same symmetric preference tie marks.
        # Replacing two singleton diagonal marks with a paired edge keeps tie
        # handling comparable to the assignment objective.
        marginal = [[2 * perturbed[i][j] - perturbed[i][i] - perturbed[j][j]
                     for j in range(len(vertices))] for i in range(len(vertices))]
        _, edges = integral_solver(marginal)((1 << len(vertices)) - 1)
        chosen = {v: (v, v) for v in vertices}
        for i, j in edges:
            chosen[vertices[i]] = (vertices[j], vertices[j])
            chosen[vertices[j]] = (vertices[i], vertices[i])
        return chosen
    return choices


def ordered_matching(weights, order, choices, draws):
    available = set(order[:WAIT])
    seen_mask = sum(1 << v for v in available)
    value, fractional_prefixes, matches = 0, 0, []
    for time, arriving in enumerate(order[WAIT:], WAIT + 1):
        seen_mask |= 1 << arriving
        options = choices(seen_mask)
        fractional_prefixes += any(a != b for a, b in options.values())
        partner = options[arriving][int(draws[time][0] >= 0.5)]
        available.add(arriving)
        if arriving == partner:
            keepout = 1 / 3 + 2 * WAIT * (WAIT - 1) * (WAIT - 2) / (
                3 * (time - 1) * (time - 2) * (time - 3))
            if draws[time][1] < keepout:
                available.remove(arriving)
        elif partner in available:
            available.remove(arriving)
            available.remove(partner)
            assert weights[arriving][partner] >= 0
            value += weights[arriving][partner] * SCALE
            matches.extend((arriving, partner))
        for vertex in sorted(available - {arriving, partner}):
            self_probability = options[vertex].count(vertex) / 2
            if draws[time][vertex + 2] < self_probability / (time - 2 + self_probability):
                available.remove(vertex)
    assert len(set(matches)) == len(matches)
    return value, fractional_prefixes


def greedy_matching(weights, order):
    free, value = [], 0
    for arriving in order:
        if free:
            partner = max(free, key=lambda v: (weights[arriving][v], -v))
            if weights[arriving][partner] > 0:
                value += weights[arriving][partner] * SCALE
                free.remove(partner)
                continue
        free.append(arriving)
    return value


def greedy_coalition(welfare, order):
    groups = []
    for arriving in order:
        bit = 1 << arriving
        if groups:
            best = max(range(len(groups)), key=lambda i: welfare[groups[i] | bit] - welfare[groups[i]])
            if welfare[groups[best] | bit] > welfare[groups[best]]:
                groups[best] |= bit
                continue
        groups.append(bit)
    return sum(welfare[g] for g in groups)


def generated(family, replicate):
    rng = Random(key_seed("matrix", family, replicate))
    labels = [0, 0, 1, 1, 1, 2, 2, 2, 2] if family == "unequal_clusters" else [i // 3 for i in range(N)]
    rng.shuffle(labels)
    weights = [[0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            same = labels[i] == labels[j]
            if family == "unit_triangles":
                w = int(same)
            elif family == "independent":
                w = rng.choice((-3, -1, 0, 1, 2, 4))
            elif family == "unequal_clusters":
                w = rng.choice((1, 2, 3, 4, 5) if same else (-4, -2, -2, 0, 3))
            else:
                w = rng.choice((2, 3, 3, 4) if same else (-3, -2, -1, 0, 0, 1))
            weights[j][i] = weights[i][j] = w
    return weights


def snap_panels(path):
    adjacency = defaultdict(set)
    with gzip.open(path, "rt") as archive:
        for line in archive:
            i, j = map(int, line.split())
            adjacency[i].add(j)
            adjacency[j].add(i)
    assert len(adjacency) == 4039 and sum(map(len, adjacency.values())) == 176468
    rng = Random(3162026)
    anchors = sorted(adjacency)
    rng.shuffle(anchors)
    used, out = set(), []
    for anchor in anchors:
        if anchor in used:
            continue
        active, reached, visited = deque([anchor]), [anchor], {anchor}
        while active and len(reached) < N:
            here = active.popleft()
            neighbours = sorted(adjacency[here] - used - visited)
            rng.shuffle(neighbours)
            for neighbour in neighbours:
                active.append(neighbour)
                reached.append(neighbour)
                visited.add(neighbour)
                if len(reached) == N:
                    break
        if len(reached) == N:
            used.update(reached)
            out.append([[0 if i == j else 2 if j in adjacency[i] else -1 for j in reached]
                        for i in reached])
            if len(out) == BATCH:
                break
    assert len(out) == BATCH and len(used) == N * BATCH
    return out


def evaluate(weights, family, replicate):
    welfare = coalition_welfares(weights)
    full = (1 << N) - 1
    offline = offline_coalitions(welfare)(full)
    matching = integral_solver(weights)(full)[0] * SCALE
    assert 0 < matching <= offline
    sums = {policy: 0 for policy in POLICIES}
    fractional_prefixes = 0
    fraction = prefix_solver(weights, (family, replicate), True)
    integral = prefix_solver(weights, (family, replicate), False)
    rng = Random(key_seed("arrivals", family, replicate))
    for _ in range(ORDERS):
        order = rng.sample(range(N), N)
        if family == "triangular_high_degree_first":
            order.sort(key=lambda i: -sum(max(w, 0) for w in weights[i]))
        draws = [[rng.random() for _ in range(N + 2)] for _ in range(N + 1)]
        fractional, prefix_count = ordered_matching(weights, order, fraction, draws)
        integer, count = ordered_matching(weights, order, integral, draws)
        assert count == 0
        fractional_prefixes += prefix_count
        greedy_pair_value = greedy_matching(weights, order)
        greedy_group_value = greedy_coalition(welfare, order)
        assert max(fractional, integer, greedy_pair_value) <= matching and greedy_group_value <= offline
        for name, score in zip(POLICIES, [fractional, integer, greedy_pair_value, greedy_group_value]):
            sums[name] += score
    return {"ratios": {name: score / (ORDERS * offline) for name, score in sums.items()},
            "offline_matching_ratio": matching / offline,
            "fractional_prefix_rate": fractional_prefixes / (ORDERS * (N - WAIT)),
            "positive_edges": sum(w > 0 for i, row in enumerate(weights) for w in row[i + 1:])}


def summarise(rows, family):
    differences = [r["ratios"]["fractional"] - r["ratios"]["integral"] for r in rows]
    rng = Random(key_seed("bootstrap", family))
    sampling = [mean(rng.choices(differences, k=len(rows))) for _ in range(5000)]
    p = quantiles(sampling, n=1000)
    return {"ratio": {policy: mean(r["ratios"][policy] for r in rows) for policy in POLICIES},
            "offline_matching_ratio": mean(r["offline_matching_ratio"] for r in rows),
            "fractional_minus_integral": mean(differences), "ci95": [p[24], p[974]],
            "fractional_prefix_rate": mean(r["fractional_prefix_rate"] for r in rows),
            "positive_edges": mean(r["positive_edges"] for r in rows),
            "graph_level": rows}


def partitions(vertices):
    if not vertices:
        yield []
        return
    first, *tail = vertices
    for blocks in partitions(tail):
        yield [[first], *blocks]
        for block in range(len(blocks)):
            yield [current + [first] if i == block else current for i, current in enumerate(blocks)]


def audit():
    exact_assignment, exact_partition, prefix_degree = 0, 0, 0
    for size in range(1, 7):
        for index in range(6):
            rng = Random(key_seed("audit", size, index))
            matrix = [[rng.randint(-6, 8) for _ in range(size)] for _ in range(size)]
            score = lambda order: sum(matrix[i][j] for i, j in enumerate(order))
            assert score(assignment(matrix)) == max(map(score, permutations(range(size))))
            exact_assignment += 1
        for index in range(4):
            weights = generated("independent", size * 100 + index)
            half_edges = prefix_solver(weights, ("audit", index), True)((1 << size) - 1)
            for i, pair in half_edges.items():
                assert len(pair) == 2
                for j in pair:
                    assert pair.count(j) == half_edges[j].count(i)
                prefix_degree += 1
            w = coalition_welfares(weights)
            evaluator = lambda partition: sum(w[sum(1 << v for v in block)] for block in partition)
            all_partitions = list(partitions(list(range(size))))
            assert offline_coalitions(w)((1 << size) - 1) == max(map(evaluator, all_partitions))
            assert integral_solver(weights)((1 << size) - 1)[0] * SCALE == max(
                evaluator(p) for p in all_partitions if all(len(group) <= 2 for group in p))
            exact_partition += 1
    return {"assignment_against_all_permutations": exact_assignment,
            "welfare_and_matching_against_enumerated_partitions": exact_partition,
            "fractional_row_degree_and_symmetry_checks": prefix_degree}


def run(snap):
    protocol = sha256(Path(__file__).with_name("MATCHING_PROTOCOL.md").read_bytes()).hexdigest()
    assert protocol == PROTOCOL_HASH
    cohorts = {}
    panels = {"signed_triangles": [generated("signed_triangles", i) for i in range(BATCH)],
              "independent": [generated("independent", i) for i in range(BATCH)],
              "unequal_clusters": [generated("unequal_clusters", i) for i in range(BATCH)],
              "unit_triangles": [generated("unit_triangles", i) for i in range(BATCH)]}
    panels["triangular_high_degree_first"] = panels["signed_triangles"]
    snap_hash = None
    if snap:
        panels["snap_ego_topology"] = snap_panels(snap)
        snap_hash = sha256(snap.read_bytes()).hexdigest()
    for family, matrices in panels.items():
        cohorts[family] = summarise([evaluate(w, family, i) for i, w in enumerate(matrices)], family)
    primary = cohorts["signed_triangles"]
    return {"protocol_sha256": protocol, "snap_archive_sha256": snap_hash,
            "graphs_per_cohort": BATCH, "internal_replicates_per_graph": ORDERS,
            "agents_per_graph": N, "audit": audit(), "cohorts": cohorts,
            "primary_gate_passed": primary["fractional_minus_integral"] >= 0.01 and primary["ci95"][0] > 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snap", type=Path)
    parser.add_argument("--write", type=Path)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.audit:
        print(json.dumps(audit(), indent=2))
        return
    result = json.dumps(run(args.snap), indent=2) + "\n"
    if args.verify:
        if args.verify.read_text() != result:
            raise ValueError("Output differs from saved matching results")
        print("Verified matching results")
    elif args.write:
        args.write.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
