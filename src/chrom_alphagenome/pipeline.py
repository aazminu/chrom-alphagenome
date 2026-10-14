from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import List


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    alphagenome_rna_q: float
    alphagenome_splice_q: float
    alphagenome_atac_q: float
    gnomad_af: float
    motif_delta_llr: float = 0.0
    clinvar_status: str = "VUS"


@dataclass(frozen=True)
class PrioritizationResult:
    variant_id: str
    gene: str
    composite_avi: float
    splice_delta_psi: float
    rarity_weight: float
    posterior_prob: float
    mechanism_class: str


class RegulatoryPrioritizer:
    """Calibrates multi-modal AlphaGenome & ChromBPNet tracks into posterior pathogenicity probabilities."""

    def __init__(self, af_threshold: float = 1e-4, prior_odds: float = 0.15) -> None:
        self.af_threshold = af_threshold
        self.prior_odds = prior_odds

    @staticmethod
    def _logit(p: float) -> float:
        clipped = min(max(p, 1e-6), 1.0 - 1e-6)
        return math.log(clipped / (1.0 - clipped))

    @staticmethod
    def _sigmoid(x: float) -> float:
        return 1.0 / (1.0 + math.exp(-x))

    def score_variant(self, v: VariantRecord) -> PrioritizationResult:
        vid = f"{v.chrom}:{v.pos}:{v.ref}>{v.alt}"
        # Weighted root-mean-power aggregation emphasizing splicing & expression extremes
        composite_avi = (
            0.40 * (v.alphagenome_splice_q ** 2)
            + 0.35 * (v.alphagenome_rna_q ** 2)
            + 0.25 * (v.alphagenome_atac_q ** 2)
        ) ** 0.5

        splice_delta_psi = max(0.0, (v.alphagenome_splice_q - 0.5) * 1.85)
        rarity_weight = 1.0 if v.gnomad_af <= 0 else math.exp(-v.gnomad_af / self.af_threshold)
        motif_bonus = min(max(-v.motif_delta_llr / 8.0, 0.0), 1.2)

        log_lr = (
            1.15 * self._logit(composite_avi)
            + 1.45 * rarity_weight
            + 0.95 * motif_bonus
            - 0.65
        )
        post_odds = self.prior_odds * math.exp(log_lr)
        posterior_prob = post_odds / (1.0 + post_odds)

        if v.alphagenome_splice_q >= 0.90:
            mechanism = "Cryptic Splice Disruption (Sashimi High-Impact)"
        elif v.alphagenome_atac_q >= 0.85 and abs(v.motif_delta_llr) >= 3.0:
            mechanism = "Cis-Regulatory TF Footprint Disruption (cCRE)"
        else:
            mechanism = "Distal Transcriptional Modulation"

        return PrioritizationResult(
            variant_id=vid,
            gene=v.gene,
            composite_avi=round(composite_avi, 5),
            splice_delta_psi=round(splice_delta_psi, 5),
            rarity_weight=round(rarity_weight, 5),
            posterior_prob=round(posterior_prob, 5),
            mechanism_class=mechanism,
        )

    def rank_cohort(self, variants: List[VariantRecord]) -> List[PrioritizationResult]:
        scored = [self.score_variant(v) for v in variants]
        return sorted(scored, key=lambda r: r.posterior_prob, reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="ChromAlphaGenome Variant Prioritizer")
    parser.add_argument("--locus", default="chr11:5225464-5227071")
    parser.add_argument("--gene", default="HBB")
    args = parser.parse_args()

    sample_cohort = [
        VariantRecord("chr11", 5226774, "C", "T", args.gene, 0.994, 0.998, 0.812, 1.2e-5, -6.4, "Pathogenic/VUS"),
        VariantRecord("chr11", 5227002, "G", "A", args.gene, 0.962, 0.420, 0.975, 0.0, -7.8, "VUS"),
        VariantRecord("chr11", 5225910, "A", "G", args.gene, 0.310, 0.215, 0.290, 4.5e-2, -0.4, "Benign"),
    ]
    prioritizer = RegulatoryPrioritizer()
    ranked = prioritizer.rank_cohort(sample_cohort)
    print(json.dumps({"locus": args.locus, "ranked_variants": [asdict(r) for r in ranked]}, indent=2))


if __name__ == "__main__":
    main()
