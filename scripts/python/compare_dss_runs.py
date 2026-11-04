#!/usr/bin/env python3
"""
MethylBench - compare the primary DSS/DMRcate run with its sensitivity runs
(Suppl. Table S7): number and direction of DMCs and DMRs per platform and tier,
pairwise DMC Jaccard and base-pair-weighted DMR Jaccard indices between
platforms, and the agreement of each platform with itself between runs.

Reads the tables written by 13_DMR_DSS_analysis.R:
  Tier{1,2}_DML_significant_<METHOD>.tsv   (chr, pos, diff, fdr)
  Tier{1,2}_DMR_<METHOD>.tsv               (seqnames, start, end, ..., meandiff)

Usage:
  python3 scripts/python/compare_dss_runs.py \
    --run primary=$OUT/dmr_dss --run unsmoothed=$OUT/dmr_dss_unsmoothed \
    [--run unpaired=$OUT/dmr_dss_unpaired] -o $OUT/dss_sensitivity_summary.tsv
The first --run is the reference for the between-run comparison.
"""
import argparse, csv, os, sys
from collections import defaultdict
from itertools import combinations

SEQ = ["ONT", "TWIST", "WGEC", "RRBS"]


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_dmc(d, tier, m):
    p = os.path.join(d, f"Tier{tier}_DML_significant_{m}.tsv")
    if not os.path.exists(p):
        return None
    rows = read_tsv(p)
    ids = {f"{r['chr']}:{r['pos']}" for r in rows}
    hyper = sum(1 for r in rows if float(r["diff"]) > 0)
    return ids, hyper


def load_dmr(d, tier, m):
    p = os.path.join(d, f"Tier{tier}_DMR_{m}.tsv")
    if not os.path.exists(p):
        return None
    rows = read_tsv(p)
    iv = defaultdict(list)
    hyper = 0
    for r in rows:
        iv[r["seqnames"]].append((int(float(r["start"])), int(float(r["end"]))))
        hyper += float(r["meandiff"]) > 0
    return merge(iv), len(rows), hyper


def merge(iv):
    out = {}
    for c, xs in iv.items():
        xs.sort()
        m = [list(xs[0])]
        for s, e in xs[1:]:
            if s <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = m
    return out


def bp(iv):
    return sum(e - s + 1 for xs in iv.values() for s, e in xs)


def bp_inter(a, b):
    tot = 0
    for c in set(a) & set(b):
        x, y, i, j = a[c], b[c], 0, 0
        while i < len(x) and j < len(y):
            lo, hi = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
            if lo <= hi:
                tot += hi - lo + 1
            if x[i][1] < y[j][1]:
                i += 1
            else:
                j += 1
    return tot


def jac_set(a, b):
    u = len(a | b)
    return len(a & b) / u if u else float("nan")


def jac_bp(a, b):
    i = bp_inter(a, b)
    u = bp(a) + bp(b) - i
    return i / u if u else float("nan")


def rng(v):
    v = [x for x in v if x == x]
    return f"{min(v):.3f}-{max(v):.3f}" if v else "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="NAME=/path/to/dmr_dss_dir")
    ap.add_argument("-o", "--out", default="dss_sensitivity_summary.tsv")
    a = ap.parse_args()
    runs = []
    for s in a.run:
        if "=" not in s:
            sys.exit(f"ERROR: --run needs NAME=/path (got '{s}')")
        n, d = s.split("=", 1)
        if not os.path.isdir(d):
            sys.exit(f"ERROR: folder not found: {d}")
        runs.append((n, d))

    data = {}
    for n, d in runs:
        for tier in (1, 2):
            for m in SEQ + (["EPIC"] if tier == 2 else []):
                data[(n, tier, m)] = (load_dmc(d, tier, m), load_dmr(d, tier, m))

    rows = []
    ref = runs[0][0]
    for tier in (1, 2):
        methods = SEQ + (["EPIC"] if tier == 2 else [])
        for n, _ in runs:
            for m in methods:
                dmc, dmr = data[(n, tier, m)]
                r = {"run": n, "tier": f"Tier{tier}", "method": m}
                if dmc:
                    r["n_DMC"] = len(dmc[0])
                    r["pct_DMC_hyper"] = round(100 * dmc[1] / len(dmc[0]), 1) if dmc[0] else ""
                if dmr:
                    r["n_DMR"] = dmr[1]
                    r["pct_DMR_hyper"] = round(100 * dmr[2] / dmr[1], 1) if dmr[1] else ""
                if n != ref:
                    rdmc, rdmr = data[(ref, tier, m)]
                    if dmc and rdmc:
                        r[f"DMC_jaccard_vs_{ref}"] = round(jac_set(dmc[0], rdmc[0]), 3)
                    if dmr and rdmr:
                        r[f"DMR_bp_jaccard_vs_{ref}"] = round(jac_bp(dmr[0], rdmr[0]), 3)
                rows.append(r)
            # pairwise ranges between platforms within the run
            seqpairs = [p for p in combinations(methods, 2) if "EPIC" not in p]
            epicpairs = [p for p in combinations(methods, 2) if "EPIC" in p]
            def pj(kind, pairs):
                vals = []
                for x, y in pairs:
                    A, B = data[(n, tier, x)], data[(n, tier, y)]
                    if kind == "DMC" and A[0] and B[0]:
                        vals.append(jac_set(A[0][0], B[0][0]))
                    if kind == "DMR" and A[1] and B[1]:
                        vals.append(jac_bp(A[1][0], B[1][0]))
                return vals
            rows.append({"run": n, "tier": f"Tier{tier}", "method": "pairs seq-seq",
                         "DMC_jaccard_range": rng(pj("DMC", seqpairs)), "DMR_bp_jaccard_range": rng(pj("DMR", seqpairs))})
            if epicpairs:
                rows.append({"run": n, "tier": f"Tier{tier}", "method": "pairs EPIC-seq",
                             "DMC_jaccard_range": rng(pj("DMC", epicpairs)), "DMR_bp_jaccard_range": rng(pj("DMR", epicpairs))})

    cols = ["run", "tier", "method", "n_DMC", "pct_DMC_hyper", "n_DMR", "pct_DMR_hyper",
            "DMC_jaccard_range", "DMR_bp_jaccard_range"]
    cols += sorted({k for r in rows for k in r} - set(cols))
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(r.get(c, "")) for c in cols))
    print(f"\nWritten: {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
