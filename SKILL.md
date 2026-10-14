---
name: chrom-alphagenome
description: Prioritize non-coding and splice-altering human genetic variants (GRCh38) by combining AlphaGenome 1Mb track quantile scores, ChromBPNet footprint attribution, motif log-odds disruption, and gnomAD allele frequency priors.
---

# ChromAlphaGenome Agent Skill

Use this skill inside Google Antigravity alongside `google-deepmind/science-skills` (`alphagenome_variant_impact_score`, `alphagenome_single_variant_analysis`, `clinvar_database`, `gnomad_database`, and `encode_ccres_database`) to score and rank candidate causal variants in regulatory and splice regions.

## Usage

```bash
uv run python -m chrom_alphagenome.pipeline --locus chr11:5225464-5227071 --gene HBB
```
