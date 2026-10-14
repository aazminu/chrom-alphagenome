from chrom_alphagenome import RegulatoryPrioritizer, VariantRecord


def test_cryptic_splice_prioritization() -> None:
    prioritizer = RegulatoryPrioritizer(af_threshold=1e-4)
    pathogenic = VariantRecord(
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
    benign = VariantRecord(
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
    )
    ranked = prioritizer.rank_cohort([benign, pathogenic])
    assert ranked[0].variant_id == "chr11:5226774:C>T"
    assert ranked[0].posterior_prob > 0.85
    assert ranked[1].posterior_prob < 0.20
    assert "Splice" in ranked[0].mechanism_class
