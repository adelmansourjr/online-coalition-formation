# Online Coalition Formation

Implementation and evaluation of the random-arrival fractional prefix-matching algorithm studied by Bullinger, Romen and Schlenga. The comparison includes integral matching, immediate greedy matching, larger coalitions and exact offline reference solutions.

[Report](MATCHING_REPORT.md) · [Protocol](MATCHING_PROTOCOL.md) · [Results](matching_results.json)

Use Python 3.10+ with the standard library.

```sh
python3 online_matching.py --audit
curl -fL https://snap.stanford.edu/data/facebook_combined.txt.gz -o facebook_combined.txt.gz
python3 online_matching.py --snap facebook_combined.txt.gz --verify matching_results.json
```

The public topology panel uses anonymised SNAP graph structure transformed into stated synthetic utility weights. The report preserves the stronger greedy baselines and limits of the paper-style methods.
