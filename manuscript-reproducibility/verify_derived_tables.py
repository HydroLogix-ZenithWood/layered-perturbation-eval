"""Verify supplied manuscript-derived tables using Python's standard library.

This checks aggregation and table consistency, not original model inference or
cell-ID intersections. All paths are relative to this file unless --data is set.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

BASE = Path(__file__).resolve().parent
METRICS = ("rho_T", "rho_G", "delta_rho", "coverage_auc_T", "coverage_auc_G",
           "delta_coverage_auc", "random_order_expected_auc")


def require(condition, label):
    if not condition:
        raise ValueError(label)


def grouped(rows, fields):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[f] for f in fields)].append(row)
    return result


def ranks(values):
    return [sum(x > y for y in values) + (sum(x == y for y in values) + 1) / 2
            for x in values]


def spearman(a, b):
    x, y = ranks(a), ranks(b)
    x, y = [v - mean(x) for v in x], [v - mean(y) for v in y]
    denominator = math.sqrt(sum(v * v for v in x) * sum(v * v for v in y))
    require(denominator > 0, "Undefined rank correlation")
    return sum(u * v for u, v in zip(x, y)) / denominator


def coverage(predictor, outcome):
    # Fractional inclusion of the complete tie block; outcomes never break ties.
    result = []
    for k, cutoff in enumerate(sorted(predictor, reverse=True), 1):
        high = [i for i, v in enumerate(predictor) if v > cutoff]
        tied = [i for i, v in enumerate(predictor) if v == cutoff]
        fraction = (k - len(high)) / len(tied)
        result.append((sum(outcome[i] for i in high)
                       + fraction * sum(outcome[i] for i in tied)) / k)
    return result


def quantile(values, probability):
    ordered = sorted(values)
    h = (len(ordered) - 1) * probability
    lo, hi = math.floor(h), math.ceil(h)
    return ordered[lo] + (h - lo) * (ordered[hi] - ordered[lo])


def verify(data):
    differences = []

    def close(actual, expected, label):
        delta = abs(float(actual) - float(expected))
        require(math.isfinite(delta) and delta < 1e-10,
                f"{label}: {actual} differs from {expected}")
        differences.append(delta)

    def read(name):
        with (data / name).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    manifest = json.loads((BASE / "input_manifest.json").read_text())
    for entry in manifest["files"]:
        path = data / Path(entry["file"]).name
        require(hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"],
                f"Input hash mismatch: {path.name}")

    summary = read("Figure4_estimand_summary.csv")
    require(summary == read("Table_S24_Figure4_estimands.csv"), "Table S24 mismatch")
    fields = ("source", "reference")
    target_groups = grouped(read("Figure4_target_mean_scores.csv"), fields)
    draw_groups = grouped(read("Figure4_draw_metrics.csv"), fields)
    curve_groups = grouped(read("Figure4_target_mean_coverage_curves.csv"), fields)
    summary_groups = grouped(summary, fields)
    require(len(target_groups) == 6, "Expected six Figure 4 conditions")
    require(target_groups.keys() == draw_groups.keys() == curve_groups.keys()
            == summary_groups.keys(), "Figure 4 condition identities differ")
    for key, target_rows in target_groups.items():
        require(len(target_rows) == len({r["target"] for r in target_rows}) == 17,
                f"{key}: gallery must contain 17 unique targets")
        draws = draw_groups[key]
        require(len(draws) == 100 and {int(r["iteration"]) for r in draws} == set(range(100)),
                f"{key}: expected 100 unique draws")
        t, g, v = ([float(r[name]) for r in target_rows] for name in ("T", "G", "V"))
        ct, cg = coverage(t, v), coverage(g, v)
        rt, rg = spearman(t, v), spearman(g, v)
        values = dict(zip(METRICS, (rt, rg, rg - rt, mean(ct), mean(cg),
                                    mean(cg) - mean(ct), mean(v))))
        entries = summary_groups[key]
        require(len(entries) == 7 and {r["metric"] for r in entries} == set(METRICS),
                f"{key}: unexpected metric identities")
        for row in entries:
            metric = row["metric"]
            require(int(row["target_count"]) == 17 and int(row["draw_count"]) == 100,
                    f"{key}: inconsistent denominators")
            distribution = [float(r[metric]) for r in draws]
            for field, actual in (("target_mean_statistic", values[metric]),
                                  ("within_draw_mean", mean(distribution)),
                                  ("within_draw_p025", quantile(distribution, .025)),
                                  ("within_draw_p975", quantile(distribution, .975))):
                close(actual, row[field], f"{key}/{metric}/{field}")
        for row in draws:
            close(float(row["rho_G"]) - float(row["rho_T"]), row["delta_rho"], "draw delta rho")
            close(float(row["coverage_auc_G"]) - float(row["coverage_auc_T"]),
                  row["delta_coverage_auc"], "draw delta coverage")
        require(len(curve_groups[key]) == 34, f"{key}: expected two 17-point curves")
        for row in curve_groups[key]:
            k = int(row["coverage_k"])
            require(row["predictor"] in {"T", "G"} and 1 <= k <= 17, "Invalid curve key")
            close((ct if row["predictor"] == "T" else cg)[k - 1], row["outcome_mean"],
                  f"{key}/curve/{row['predictor']}/{k}")
            close(mean(v), row["uniform_order_expectation"], "Uniform-order expectation")

    current = read("microglia_current_reference_overlap_recomputed.csv")
    current_groups = grouped(current, ("experiment", "source", "readout"))
    overlaps = read("Table_S25_microglia_overlap_summary.csv")
    require(len(current) == 79200 and len(current_groups) == len(overlaps) == 20,
            "Unexpected current overlap record counts")
    for row in overlaps:
        records = current_groups[(row["experiment"], row["source"], row["readout"])]
        keys = [(r["target"], r["iteration"], r["direction"]) for r in records]
        require(len(keys) == len(set(keys)), "Duplicate target/direction/draw key")
        shared = [int(r["n_shared_native_start_reference_cells"]) for r in records]
        starts = [int(r["n_start_cells"]) for r in records]
        require(all(0 <= a <= b for a, b in zip(shared, starts)), "Invalid overlap counts")
        expected = {"K": len({r["target"] for r in records}),
                    "n_target_direction_draw_rows": len(records),
                    "n_rows_with_overlap": sum(x > 0 for x in shared),
                    "fraction_rows_with_overlap": mean(x > 0 for x in shared),
                    "min_shared": min(shared), "max_shared": max(shared),
                    "mean_shared": mean(shared), "min_start": min(starts), "max_start": max(starts)}
        for name, value in expected.items():
            close(value, row[name], f"Current overlap/{row['experiment']}/{row['readout']}/{name}")
        # unique_start_cells needs original identities and is not recomputed here.

    earlier = read("Table_S26_earlier_query_exclusion.csv")
    require(len(earlier) == 324, "Unexpected earlier query count")
    for row in earlier:
        require(int(row["n_shared_start_reference"]) == int(row["n_shared_start_treated"]) == 0,
                "Earlier query-specific exclusion table is inconsistent")
        require(int(row["n_shared_start_original_all_control"]) == int(row["n_start"]),
                "Earlier starts must belong to original control set")
    external = read("Table_S27_external_role_overlap.csv")
    require(len(external) == 30, "Unexpected external query count")
    for row in external:
        require(int(row["n_start"]) == 16 and int(row["n_start_guides"]) == 2,
                "External native budget differs")
        require(all(int(row[k]) == 0 for k in ("n_start_reference_overlap", "n_start_PCA_overlap",
                                               "n_start_treated_overlap")),
                "External disjoint-role table is inconsistent")
    galleries = read("Table_S28_microglia_gallery_versions.csv")
    require(len(galleries) == 14, "Unexpected gallery comparison count")
    for row in galleries:
        for label in ("old", "current"):
            names = [n for n in row[label + "_only"].split(";") if n]
            require(len(names) == len(set(names)), "Duplicate gallery-only target")
            require(int(row[label + "_K"]) == int(row["intersection"]) + len(names),
                    "Gallery partition counts disagree")

    return {"status": "passed", "scope": "Supplied derived-table aggregation and consistency only; no model inference or fresh cell-ID intersections.",
            "figure4_conditions": 6, "figure4_draws_per_condition": 100,
            "figure4_estimands": 42, "numeric_comparisons": len(differences),
            "max_absolute_difference": max(differences), "current_overlap_records": len(current),
            "earlier_query_rows": len(earlier), "external_query_rows": len(external),
            "gallery_version_rows": len(galleries), "input_hashes_verified": len(manifest["files"]),
            "not_recomputed": ["Original model outputs", "Original score matrices", "Original cell-ID intersections", "Unique native starting-cell union counts"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=BASE / "data")
    parser.add_argument("--output", type=Path, help="Optional new JSON report path; existing paths are rejected")
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error("Output already exists")
    report = verify(args.data)
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
