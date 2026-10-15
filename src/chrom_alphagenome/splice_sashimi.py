from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

VALID_BASES = {"A", "C", "G", "T"}

# 9-mer 5' splice donor PWM (3 exonic + 6 intronic: M A G | G T R A G T)
DONOR_PWM: List[Dict[str, float]] = [
    {"A": 0.35, "C": 0.35, "G": 0.18, "T": 0.12},
    {"A": 0.62, "C": 0.12, "G": 0.14, "T": 0.12},
    {"A": 0.09, "C": 0.04, "G": 0.80, "T": 0.07},
    {"A": 0.005, "C": 0.005, "G": 0.985, "T": 0.005},  # +1 G
    {"A": 0.005, "C": 0.005, "G": 0.005, "T": 0.985},  # +2 T
    {"A": 0.54, "C": 0.03, "G": 0.40, "T": 0.03},     # +3 R (A/G)
    {"A": 0.71, "C": 0.08, "G": 0.12, "T": 0.09},     # +4 A
    {"A": 0.07, "C": 0.05, "G": 0.82, "T": 0.06},     # +5 G
    {"A": 0.15, "C": 0.15, "G": 0.18, "T": 0.52},     # +6 T
]

# 9-mer 3' splice acceptor core PWM (6 polypyrimidine/branch + C A G | G)
ACCEPTOR_PWM: List[Dict[str, float]] = [
    {"A": 0.08, "C": 0.42, "G": 0.06, "T": 0.44},
    {"A": 0.07, "C": 0.44, "G": 0.05, "T": 0.44},
    {"A": 0.08, "C": 0.43, "G": 0.06, "T": 0.43},
    {"A": 0.09, "C": 0.41, "G": 0.08, "T": 0.42},
    {"A": 0.24, "C": 0.26, "G": 0.24, "T": 0.26},     # N
    {"A": 0.04, "C": 0.74, "G": 0.02, "T": 0.20},     # -3 C
    {"A": 0.985, "C": 0.005, "G": 0.005, "T": 0.005},  # -2 A
    {"A": 0.005, "C": 0.005, "G": 0.985, "T": 0.005},  # -1 G
    {"A": 0.22, "C": 0.14, "G": 0.52, "T": 0.12},     # +1 G
]


@dataclass(frozen=True)
class SashimiDisruptionReport:
    variant_id: str
    canonical_delta_psi: float
    cryptic_activated: bool
    cryptic_offset: int
    ref_maxent_score: float
    alt_maxent_score: float
    sashimi_exon_inclusion_ref: float
    sashimi_exon_inclusion_alt: float


class MaxEntSashimiAnalyzer:
    """Quantifies 5' donor / 3' acceptor MaxEntScan log-odds and Sashimi delta-PSI splice disruption."""

    @staticmethod
    def _validate_seq(seq: str) -> str:
        clean = seq.strip().upper()
        if not clean or any(ch not in VALID_BASES for ch in clean):
            raise ValueError("Sequence must be non-empty and contain only A, C, G, T")
        return clean

    @classmethod
    def score_donor_9mer(cls, kmer: str) -> float:
        clean = cls._validate_seq(kmer)
        if len(clean) != 9:
            raise ValueError("Donor motif must be exactly 9 bases (3 exonic + 6 intronic)")
        score = 0.0
        for i, base in enumerate(clean):
            prob = DONOR_PWM[i][base]
            score += math.log2(prob / 0.25)
        # MaxEntScan consensus dinucleotide term at +1/+2 (indices 3, 4)
        if clean[3:5] != "GT":
            score -= 6.5
        return round(score, 4)

    @classmethod
    def score_acceptor_9mer(cls, kmer: str) -> float:
        clean = cls._validate_seq(kmer)
        if len(clean) != 9:
            raise ValueError("Acceptor core motif must be exactly 9 bases")
        score = 0.0
        for i, base in enumerate(clean):
            prob = ACCEPTOR_PWM[i][base]
            score += math.log2(prob / 0.25)
        # MaxEntScan consensus dinucleotide term at -2/-1 (indices 6, 7)
        if clean[6:8] != "AG":
            score -= 6.5
        return round(score, 4)

    @staticmethod
    def _psi_from_score(maxent_score: float, midpoint: float = 4.5, slope: float = 0.65) -> float:
        z = slope * (maxent_score - midpoint)
        return 1.0 / (1.0 + math.exp(-max(min(z, 30.0), -30.0)))

    def evaluate_splice_window(
        self,
        variant_id: str,
        ref_window: str,
        alt_window: str,
        canonical_offset: int,
        site_type: str = "donor",
    ) -> SashimiDisruptionReport:
        ref_clean = self._validate_seq(ref_window)
        alt_clean = self._validate_seq(alt_window)
        if len(ref_clean) != len(alt_clean) or len(ref_clean) < 9:
            raise ValueError("ref_window and alt_window must have equal length >= 9")
        if canonical_offset < 0 or canonical_offset + 9 > len(ref_clean):
            raise ValueError("canonical_offset out of bounds for 9-mer window")
        if site_type not in ("donor", "acceptor"):
            raise ValueError("site_type must be 'donor' or 'acceptor'")

        scorer = self.score_donor_9mer if site_type == "donor" else self.score_acceptor_9mer
        ref_score = scorer(ref_clean[canonical_offset : canonical_offset + 9])
        alt_score = scorer(alt_clean[canonical_offset : canonical_offset + 9])

        psi_ref = self._psi_from_score(ref_score)
        psi_alt = self._psi_from_score(alt_score)

        # Scan for de-novo cryptic splice activation across all non-canonical offsets
        best_cryptic_gain = 0.0
        best_cryptic_offset = canonical_offset
        cryptic_activated = False
        for offset in range(len(ref_clean) - 8):
            if offset == canonical_offset:
                continue
            r_s = scorer(ref_clean[offset : offset + 9])
            a_s = scorer(alt_clean[offset : offset + 9])
            gain = a_s - r_s
            if a_s >= 5.0 and gain >= 2.5 and gain > best_cryptic_gain:
                best_cryptic_gain = gain
                best_cryptic_offset = offset
                cryptic_activated = True

        # If a strong cryptic site is activated, effective canonical exon inclusion drops further
        if cryptic_activated:
            cryptic_psi = self._psi_from_score(
                scorer(alt_clean[best_cryptic_offset : best_cryptic_offset + 9])
            )
            psi_alt = psi_alt * (1.0 - 0.55 * cryptic_psi)

        delta_psi = psi_alt - psi_ref
        return SashimiDisruptionReport(
            variant_id=variant_id,
            canonical_delta_psi=round(delta_psi, 5),
            cryptic_activated=cryptic_activated,
            cryptic_offset=best_cryptic_offset - canonical_offset,
            ref_maxent_score=round(ref_score, 4),
            alt_maxent_score=round(alt_score, 4),
            sashimi_exon_inclusion_ref=round(psi_ref, 5),
            sashimi_exon_inclusion_alt=round(psi_alt, 5),
        )
