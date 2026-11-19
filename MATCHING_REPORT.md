# When fractional prefix matching helps online welfare

Replacing the integral prefix matching with the fractional optimum improved the average fraction of exact offline welfare by **3.43 percentage points** on the prespecified 48 signed triangular-cluster games. The paired graph bootstrap interval was **[2.71, 4.18]** points. On panels from an independent public SNAP topology, the same ablation gained **2.23 points** with a conditional interval of **[1.61, 2.87]**.

| Cohort | Fractional | Integral | Greedy pair | Greedy coalition | Fractional minus integral |
|---|---:|---:|---:|---:|---:|
| Signed triangles, primary | 41.15% | 37.71% | 52.12% | 91.26% | +3.43 pp |
| Independent signed weights | 48.44% | 47.67% | 62.86% | 82.13% | +0.77 pp |
| Unequal dense clusters | 47.39% | 46.09% | 63.95% | 87.77% | +1.30 pp |
| Repeated equal triangles | 40.82% | 34.47% | 50.00% | 100.00% | +6.35 pp |
| Signed triangles, degree-first stress | 41.16% | 36.58% | 54.81% | 82.33% | +4.58 pp |
| SNAP topology proxy | 48.00% | 45.77% | 68.46% | 95.02% | +2.23 pp |

Both paper-style algorithms used the same waiting period, information and availability randomness. Greedy pair and coalition start assigning immediately and were stronger practical comparators on these finite families. The offline pair-only ceiling was **63.75%** of unrestricted offline welfare on the main cohort. Restricting coalition size is a meaningful distinction, independently of decision quality. [Frozen protocol](MATCHING_PROTOCOL.md) · [Results](matching_results.json) · [SODA 2026 original](https://arxiv.org/html/2505.06163v2).

## Exact implementation

This implements Bullinger, Romen and Schlenga's random-arrival Algorithm 1 with k = 3 and its Appendix B.3 integral substitution. Both algorithms sample a prefix-optimum partner and pair only with an available participant. Availability draws match the paper's stochastic transitions. Tie resolution depends on the seen set and seed, independent of order and availability.

Fill the slack at each vertex of a fractional matching with a zero-weight self-loop. A symmetric fractional solution then gives a doubly stochastic matrix. Conversely, symmetrising any doubly stochastic matrix gives a feasible symmetric fractional matching once its diagonal is removed. Symmetric weights make both objectives identical, accounting for twice-counted edges. The Birkhoff–von Neumann theorem ensures a maximum-weight permutation is a directed optimum. The Hungarian algorithm computes it exactly and symmetrisation gives a half-integral solution. Symmetric integer tie marks cannot change the primary optimum. Independent subset recurrences compute integral matching and arbitrary coalition partitions. Small instances check assignment against all permutations, both offline optima against partition enumeration, and fractional row constraints and symmetry.

SNAP supplies an independent anonymised topology. The benchmark constructs weights by mapping edges to +2 and other pairs to −1. The 48 disjoint panels derive from one social dataset and the interval describes sensitivity to panel sampling under these assumptions. The saved outputs report aggregate panel performance. [Stanford SNAP source](https://snap.stanford.edu/data/ego-Facebook.html).

## Replay and follow-up

```sh
curl -fL https://snap.stanford.edu/data/facebook_combined.txt.gz -o /tmp/facebook_combined.txt.gz
python3 online_matching.py --audit
python3 online_matching.py --snap /tmp/facebook_combined.txt.gz --verify matching_results.json
```

Python 3.10+ standard library. The archive argument enables the external topology cohort. The full version has replayed byte for byte. Intervals resample whole graphs, with all 32 simulated arrival orders grouped by graph.

The comparison isolates the fractional relaxation. The practical heuristic remains stronger on these small families. A useful next question for Bullinger is whether adapting pair versus larger-group choice can preserve a worst-case guarantee. This study tests aggregate fractional welfare. The earlier stability work uses different utilities and an incentive-based goal.
