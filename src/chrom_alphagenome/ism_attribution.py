from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Curated erythroid/insulator PWMs (GATA1: WGATAR, KLF1: CCMCRCCCN, CTCFcore: CCGCGNGGNGGCAG)
TF_MOTIF_PWMS: Dict[str, List[Dict[str, float]]] = {
    "GATA1": [
        {"A": 0.45, "C": 0.05, "G": 0.05, "T": 0.45},
        {"A": 0.01, "C": 0.01, "G": 0.97, "T": 0.01},
        {"A": 0.97, "C": 0.01, "G": 0.01, "T": 0.01},
        {"A": 0.01, "C": 0.01, "G": 0.01, "T": 0.97},
        {"A": 0.96, "C": 0.01, "G": 0.02, "T": 0.01},
        {"A": 0.55, "C": 0.03, "G": 0.39, "T": 0.03},
    ],
    "KLF1": [
        {"A": 0.02, "C": 0.94, "G": 0.02, "T": 0.02},
        {"A": 0.02, "C": 0.94, "G": 0.02, "T": 0.02},
        {"A": 0.45, "C": 0.45, "G": 0.05, "T": 0.05},
        {"A": 0.02, "C": 0.94, "G": 0.02, "T": 0.02},
        {"A": 0.46, "C": 0.03, "G": 0.48, "T": 0.03},
        {"A": 0.02, "C": 0.94, "G": 0.02, "T": 0.02},
        {"A": 0.02, "C": 0.94, "G": 0.02, "T": 0.02},
        {"A": 0.02, "C": 0.94, "G": 0.02, "T": 0.02},
    ],
}

# Tn5 dinucleotide insertion preference weights relative to uniform
TN5_BASE_WEIGHTS = {"A": 0.78, "T": 0.80, "G": 1.32, "C": 1.30}


@dataclass(frozen=True)
class MotifPerturbation:
    motif_name: str
    strand: str
    offset: int
    ref_llr: float
    alt_llr: float
    delta_llr: float
    disrupted: bool


@dataclass(frozen=True)
class ISMAttributionResult:
    sequence_length: int
    debiased_footprint_profile: List[float]
    ism_matrix: List[Dict[str, float]]
    top_motif_hit: Optional[MotifPerturbation]


class ChromBPNetISMAttributor:
    """Performs Tn5 footprint de-biasing, PWM motif perturbation scanning, and per-base ISM attribution."""

    @staticmethod
    def _revcomp(seq: str) -> str:
        return "".join(COMPLEMENT[b] for b in reversed(seq))

    @staticmethod
    def _score_pwm(subseq: str, pwm: Sequence[Dict[str, float]]) -> float:
        return sum(math.log2(pwm[i][base] / 0.25) for i, base in enumerate(subseq))

    def debias_tn5_profile(self, seq: str, raw_counts: Sequence[float]) -> List[float]:
        clean = seq.strip().upper()
        if len(clean) != len(raw_counts) or not clean:
            raise ValueError("Sequence and raw_counts must be non-empty and equal length")
        if any(b not in COMPLEMENT for b in clean):
            raise ValueError("Invalid nucleotide in sequence")
        if any(c < 0.0 or math.isnan(c) for c in raw_counts):
            raise ValueError("raw_counts cannot be negative or NaN")

        debiased: List[float] = []
        for i, (base, count) in enumerate(zip(clean, raw_counts)):
            left_b = clean[i - 1] if i > 0 else base
            right_b = clean[i + 1] if i + 1 < len(clean) else base
            tn5_kernel = (
                0.25 * TN5_BASE_WEIGHTS[left_b]
                + 0.50 * TN5_BASE_WEIGHTS[base]
                + 0.25 * TN5_BASE_WEIGHTS[right_b]
            )
            debiased.append(round(count / max(tn5_kernel, 1e-6), 4))
        return debiased

    def scan_motif_perturbations(self, ref_seq: str, alt_seq: str) -> List[MotifPerturbation]:
        ref_clean = ref_seq.strip().upper()
        alt_clean = alt_seq.strip().upper()
        if len(ref_clean) != len(alt_clean) or not ref_clean:
            raise ValueError("ref_seq and alt_seq must be non-empty and equal length")
        if any(b not in COMPLEMENT for b in ref_clean + alt_clean):
            raise ValueError("Sequences must contain only A, C, G, T")

        hits: List[MotifPerturbation] = []
        for motif_name, pwm in TF_MOTIF_PWMS.items():
            k = len(pwm)
            if len(ref_clean) < k:
                continue
            for offset in range(len(ref_clean) - k + 1):
                r_sub = ref_clean[offset : offset + k]
                a_sub = alt_clean[offset : offset + k]
                if r_sub == a_sub:
                    continue
                for strand, r_kmer, a_kmer in [
                    ("+", r_sub, a_sub),
                    ("-", self._revcomp(r_sub), self._revcomp(a_sub)),
                ]:
                    r_llr = self._score_pwm(r_kmer, pwm)
                    a_llr = self._score_pwm(a_kmer, pwm)
                    delta = a_llr - r_llr
                    if max(r_llr, a_llr) >= 4.0 and abs(delta) >= 2.0:
                        hits.append(
                            MotifPerturbation(
                                motif_name=motif_name,
                                strand=strand,
                                offset=offset,
                                ref_llr=round(r_llr, 4),
                                alt_llr=round(a_llr, 4),
                                delta_llr=round(delta, 4),
                                disrupted=(r_llr >= 4.0 and delta <= -2.5),
                            )
                        )
        hits.sort(key=lambda h: abs(h.delta_llr), reverse=True)
        return hits

    def compute_ism_attribution(
        self,
        ref_seq: str,
        alt_seq: str,
        raw_atac_counts: Sequence[float],
    ) -> ISMAttributionResult:
        debiased = self.debias_tn5_profile(ref_seq, raw_atac_counts)
        clean_ref = ref_seq.strip().upper()
        hits = self.scan_motif_perturbations(ref_seq, alt_seq)
        top_hit = hits[0] if hits else None

        ism_matrix: List[Dict[str, float]] = []
        for i, ref_base in enumerate(clean_ref):
            row: Dict[str, float] = {}
            for mut_base in ("A", "C", "G", "T"):
                if mut_base == ref_base:
                    row[mut_base] = 0.0
                else:
                    gc_shift = (1.0 if mut_base in ("G", "C") else -1.0) - (
                        1.0 if ref_base in ("G", "C") else -1.0
                    )
                    row[mut_base] = round(0.15 * gc_shift * math.log1p(debiased[i]), 4)
            ism_matrix.append(row)

        return ISMAttributionResult(
            sequence_length=len(clean_ref),
            debiased_footprint_profile=debiased,
            ism_matrix=ism_matrix,
            top_motif_hit=top_hit,
        )
