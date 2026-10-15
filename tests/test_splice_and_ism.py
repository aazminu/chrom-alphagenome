import pytest
from chrom_alphagenome import (
    ChromBPNetISMAttributor,
    MaxEntSashimiAnalyzer,
)


def test_maxent_donor_and_cryptic_sashimi_disruption() -> None:
    analyzer = MaxEntSashimiAnalyzer()
    # Canonical 5' donor: CAG|GTAAGT -> mutated to CAG|ATAAGT (+1 G>A) plus upstream cryptic activation
    ref_win = "ACGTACAGGTAAGTACGTACGTTTAGT"
    alt_win = "ACGTACAGATAAGTACAGGTAAGTAGT"
    report = analyzer.evaluate_splice_window(
        variant_id="chr11:5226774:G>A",
        ref_window=ref_win,
        alt_window=alt_win,
        canonical_offset=5,
        site_type="donor",
    )
    assert report.ref_maxent_score > 7.5
    assert report.alt_maxent_score < 2.0
    assert report.cryptic_activated is True
    assert report.canonical_delta_psi < -0.50


def test_maxent_acceptor_and_validation_errors() -> None:
    analyzer = MaxEntSashimiAnalyzer()
    strong_acc = "TTTTTCAGG"
    score = analyzer.score_acceptor_9mer(strong_acc)
    assert score > 5.0

    with pytest.raises(ValueError):
        analyzer.score_donor_9mer("ACGT")
    with pytest.raises(ValueError):
        analyzer.score_donor_9mer("ACGTNNGTA")


def test_chrombpnet_tn5_debiasing_and_gata1_motif_disruption() -> None:
    attributor = ChromBPNetISMAttributor()
    # Contains canonical GATA1 motif AGATAA at offset 3
    ref_seq = "CCCAGATAAGGG"
    alt_seq = "CCCACATAAGGG"  # G>C disrupts core GATA
    raw_counts = [10.0, 12.0, 15.0, 42.0, 55.0, 48.0, 50.0, 44.0, 38.0, 18.0, 14.0, 11.0]

    res = attributor.compute_ism_attribution(ref_seq, alt_seq, raw_counts)
    assert res.sequence_length == 12
    assert len(res.debiased_footprint_profile) == 12
    assert res.top_motif_hit is not None
    assert res.top_motif_hit.motif_name == "GATA1"
    assert res.top_motif_hit.disrupted is True
    assert res.top_motif_hit.delta_llr < -4.0
    # Position 4 is core G in AGATAA; G>C transversion must have strongly negative ISM attribution
    assert res.ism_matrix[4]["G"] == 0.0
    assert res.ism_matrix[4]["C"] < -5.0
    assert res.ism_matrix[4]["A"] < -5.0
    assert res.ism_matrix[4]["T"] < -5.0


def test_ism_edge_cases() -> None:
    attributor = ChromBPNetISMAttributor()
    with pytest.raises(ValueError):
        attributor.debias_tn5_profile("", [])
    with pytest.raises(ValueError):
        attributor.debias_tn5_profile("ACGT", [1.0, -2.0, 3.0, 4.0])

    analyzer = MaxEntSashimiAnalyzer()
    with pytest.raises(ValueError, match="site_type"):
        analyzer.evaluate_splice_window("v1", "ACGTACAGGTAAGT", "ACGTACAGGTAAGT", 4, site_type="invalid")


def test_cli_benchmark_markdown_json_and_reproducibility(
    capsys: pytest.CaptureFixture[str], tmp_path: pytest.TempPathFactory
) -> None:
    import json
    import subprocess
    import sys
    from pathlib import Path
    from chrom_alphagenome.__main__ import format_markdown_table, main, run_benchmark

    rows_1 = run_benchmark(seed=42)
    rows_2 = run_benchmark(seed=42)
    assert rows_1 == rows_2
    assert len(rows_1) == 3

    md = format_markdown_table(rows_1)
    assert "| benchmark | locus | primary_metric | score |" in md
    assert "maxent_sashimi_donor_disruption" in md
    assert "chrombpnet_tn5_ism_attribution" in md

    out_path = Path(str(tmp_path)) / "chrom_bench.json"
    rc = main(["--format", "json", "--seed", "42", "--output", str(out_path)])
    assert rc == 0
    data = json.loads(capsys.readouterr().out)
    assert data["seed"] == 42
    assert data["benchmarks"] == rows_1
    assert json.loads(out_path.read_text(encoding="utf-8")) == data

    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, str(repo_root / "scripts" / "benchmark.py"), "--format", "markdown"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "regulatory_cohort_prioritization" in proc.stdout

    with pytest.raises(SystemExit):
        main(["--format", "invalid_fmt"])


