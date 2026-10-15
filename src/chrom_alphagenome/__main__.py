from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence

from .ism_attribution import ChromBPNetISMAttributor
from .pipeline import RegulatoryPrioritizer, VariantRecord
from .splice_sashimi import MaxEntSashimiAnalyzer


def run_benchmark(seed: int = 42) -> List[Dict[str, Any]]:
    """Run reproducible regulatory prioritization, MaxEntScan Sashimi splice, and ChromBPNet ISM benchmarks."""
    _ = seed

    # 1. Regulatory variant cohort prioritization
    prioritizer = RegulatoryPrioritizer(af_threshold=1e-4)
    ranked = prioritizer.rank_cohort(
        [
            VariantRecord(
                chrom="chr11",
                pos=5225910,
                ref="A",
                alt="G",
                gene="HBB",
                alphagenome_rna_q=0.310,
                alphagenome_splice_q=0.215,
                alphagenome_atac_q=0.290,
                gnomad_af=0.045,
                motif_delta_llr=-0.4,
            ),
            VariantRecord(
                chrom="chr11",
                pos=5226774,
                ref="C",
                alt="T",
                gene="HBB",
                alphagenome_rna_q=0.994,
                alphagenome_splice_q=0.998,
                alphagenome_atac_q=0.812,
                gnomad_af=1.2e-5,
                motif_delta_llr=-6.4,
            ),
        ]
    )
    top_var = ranked[0]

    # 2. MaxEntScan Sashimi splice disruption
    sashimi = MaxEntSashimiAnalyzer()
    splice_rep = sashimi.evaluate_splice_window(
        variant_id="chr11:5226774:G>A",
        ref_window="ACGTACAGGTAAGTACGTACGTTTAGT",
        alt_window="ACGTACAGATAAGTACAGGTAAGTAGT",
        canonical_offset=5,
        site_type="donor",
    )

    # 3. ChromBPNet Tn5-debiased ISM & TF motif disruption
    attributor = ChromBPNetISMAttributor()
    ism_res = attributor.compute_ism_attribution(
        "CCCAGATAAGGG",
        "CCCACATAAGGG",
        [10.0, 12.0, 15.0, 42.0, 55.0, 48.0, 50.0, 44.0, 38.0, 18.0, 14.0, 11.0],
    )
    motif_name = ism_res.top_motif_hit.motif_name if ism_res.top_motif_hit else "NONE"
    motif_llr = ism_res.top_motif_hit.delta_llr if ism_res.top_motif_hit else 0.0

    return [
        {
            "benchmark": "regulatory_cohort_prioritization",
            "locus": top_var.variant_id,
            "primary_metric": "posterior_prob",
            "score": round(top_var.posterior_prob, 4),
            "secondary_metric": "composite_avi",
            "secondary_value": round(top_var.composite_avi, 4),
            "mechanism": top_var.mechanism_class,
        },
        {
            "benchmark": "maxent_sashimi_donor_disruption",
            "locus": splice_rep.variant_id,
            "primary_metric": "canonical_delta_psi",
            "score": round(splice_rep.canonical_delta_psi, 4),
            "secondary_metric": "ref_maxent_score",
            "secondary_value": round(splice_rep.ref_maxent_score, 4),
            "mechanism": "CrypticDonorActivated" if splice_rep.cryptic_activated else "CanonicalIntact",
        },
        {
            "benchmark": "chrombpnet_tn5_ism_attribution",
            "locus": f"HBB_enhancer:{motif_name}",
            "primary_metric": "motif_delta_llr",
            "score": round(motif_llr, 4),
            "secondary_metric": "core_G4C_ism",
            "secondary_value": round(ism_res.ism_matrix[4]["C"], 4),
            "mechanism": "MotifDisrupted" if (ism_res.top_motif_hit and ism_res.top_motif_hit.disrupted) else "Intact",
        },
    ]


def format_markdown_table(rows: Sequence[Dict[str, Any]]) -> str:
    headers = [
        "benchmark",
        "locus",
        "primary_metric",
        "score",
        "secondary_metric",
        "secondary_value",
        "mechanism",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row[h]) for h in headers) + " |")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run reproducible chrom-alphagenome non-coding variant benchmarks."
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
        help="Output format for benchmark results.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible benchmark execution.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional file path to write benchmark output.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    rows = run_benchmark(seed=args.seed)
    rendered = (
        json.dumps({"seed": args.seed, "benchmarks": rows}, indent=2)
        if args.format == "json"
        else format_markdown_table(rows)
    )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
