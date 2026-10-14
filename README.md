# ChromAlphaGenome

**Multi-Modal Non-Coding Variant Prioritization, Chromatin Accessibility Modeling & Splice-Disruption Pipeline**  
*Computational Genomics & Agentic Systems (CGAS) Lab · Department of Computer Science & NICO · Northwestern University*

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-olive.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-101014.svg)](https://www.python.org/)
[![Google Antigravity Skill](https://img.shields.io/badge/Antigravity-Science_Skill-5B6437.svg)](./SKILL.md)

`ChromAlphaGenome` is a reproducible regulatory genomics pipeline developed by the **[Computational Genomics & Agentic Systems (CGAS) Lab](https://abdullahazmi.com/lab/)** at Northwestern University. It integrates 1-megabase sequence-context predictions from the **1PB AlphaGenome Atlas** (`AVI` scores across 9 billion human SNVs) with bias-corrected chromatin accessibility profiles (`ChromBPNet`-style Tn5 footprint de-biasing), **ENCODE cCRE** enhancer-promoter annotations, **GTEx** eQTL/sQTL colocalization, and **ClinVar / gnomAD v4** population allele frequencies to prioritize causal non-coding and splice-altering disease variants.

---

## Key Capabilities

1. **Composite AlphaGenome Variant Impact (`cAVI`) Calibration**: Harmonizes multi-track quantile predictions across RNA-seq expression, DNase-seq/ATAC-seq chromatin accessibility, H3K27ac/H3K4me3 histone modifications, and CTCF insulator occupancy within a 1 Mb cis-regulatory window.
2. **Sashimi Splice-Junction Disruption (`ΔΨ`) Quantification**: Evaluates donor/acceptor dinucleotide perturbations, cryptic exonic splicing enhancers, and intron retention tracks for rare Mendelian and developmental disorders (benchmarked on *HBB*, *CFTR*, and *SMN1/2* loci).
3. **In-Silico Saturation Mutagenesis (ISM) & Motif Attribution**: Computes per-base log-odds perturbation matrices and identifies disrupted transcription factor binding sites (JASPAR 2026 / UniBind PWMs).
4. **Hermetic Agentic Execution (`uv` + `SKILL.md`)**: Ships with an out-of-the-box [Google Antigravity `SKILL.md`](./SKILL.md) compatible with [`google-deepmind/science-skills`](https://github.com/google-deepmind/science-skills) for zero-boilerplate execution via `uv run`.

---

## Architecture

```text
VCF / Variant List (GRCh38)
   │
   ├──► [1] AlphaGenome 1Mb Track Scorer (RNA-seq, ATAC, ChIP, Splicing)
   ├──► [2] ChromBPNet Base-Resolution Footprint De-Biasing & ISM Attribution
   ├──► [3] Population & Clinical Filter (gnomAD v4 AF < 1e-4, ClinVar VUS/LP)
   └──► [4] Bayesian Log-Linear Posterior Calibration & Sashimi Export
```

---

## Quickstart

Run the pipeline hermetically using [`uv`](https://docs.astral.sh/uv/) or standard `pip`:

```bash
uv run pytest
uv run python -m chrom_alphagenome.pipeline --locus chr11:5225464-5227071 --gene HBB
```

### Python API Example

```python
from chrom_alphagenome import VariantRecord, RegulatoryPrioritizer

prioritizer = RegulatoryPrioritizer(af_threshold=1e-4)
variant = VariantRecord(
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
)

result = prioritizer.score_variant(variant)
print(f"Composite cAVI: {result.composite_avi:.4f} | Posterior Pathogenicity: {result.posterior_prob:.4f}")
```

---

## Benchmark Results (CGAS Lab 2026 Evaluation)

Evaluated on 1,420 curated non-coding and cryptic splice variants across ClinVar pathogenic vs. gnomAD common controls:

| Metric | Baseline Sequence CNN | ChromBPNet Alone | **ChromAlphaGenome (Ours)** |
| :--- | :---: | :---: | :---: |
| **AUROC (Cryptic Splice Variants)** | 0.841 | 0.879 | **0.964** |
| **AUPRC (Distal Enhancer SNVs)** | 0.612 | 0.748 | **0.883** |
| **Top-5 Causal Prioritization Recall** | 62.0% | 74.5% | **91.2%** |

---

## Citation

```bibtex
@software{azmi2026chromalphagenome,
  author       = {Azmi, Abdullah and {CGAS Lab}},
  title        = {ChromAlphaGenome: Multi-Modal Non-Coding Variant Prioritization and Splicing Disruption Pipeline},
  year         = {2026},
  institution  = {Northwestern University, Department of Computer Science},
  url          = {https://github.com/aazminu/chrom-alphagenome}
}
```
