import subprocess
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "rfcount_parse_summary.awk"


def _parse(tmp_path, log_text, *, sample, is_map, match_mode):
    log = tmp_path / "rfcount.log"
    log.write_text(log_text)
    result = subprocess.run(
        [
            "awk",
            "-f",
            str(SCRIPT),
            "-v",
            f"sample={sample}",
            "-v",
            f"is_map={is_map}",
            "-v",
            f"match_mode={match_mode}",
            str(log),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_parse_rt_stop_exact_match_uses_last_matching_row(tmp_path):
    out = _parse(
        tmp_path,
        "\n".join(
            [
                "other 9 1.0 2.0 3.0 4.0",
                "sample1 3 10.0 20.0 30.0 40.0",
                "sample1 4 11.0 21.0 31.0 41.0",
            ]
        ),
        sample="sample1",
        is_map="0",
        match_mode="exact",
    )

    assert out == "\n".join(
        [
            "sample\tcovered\tpct_a_muts\tpct_c_muts\tpct_g_muts\tpct_u_muts",
            "sample1\t4\t11.0\t21.0\t31.0\t41.0",
            "",
        ]
    )


def test_parse_map_exact_match_extracts_mutation_percentage(tmp_path):
    out = _parse(
        tmp_path,
        "sample1 8 3/10 (30.0%) 1.0 2.0 3.0 4.0\n",
        sample="sample1",
        is_map="1",
        match_mode="exact",
    )

    assert out == "\n".join(
        [
            "sample\tcovered\tmutated_alignments\tpct_mutated\tpct_a_muts\tpct_c_muts\tpct_g_muts\tpct_u_muts",
            "sample1\t8\t3/10\t30.0\t1.0\t2.0\t3.0\t4.0",
            "",
        ]
    )


def test_parse_genome_prefix_match_emits_requested_sample_id(tmp_path):
    out = _parse(
        tmp_path,
        "sample1.sorted.bam 0 12.5 22.5 32.5 42.5\n",
        sample="sample1",
        is_map="0",
        match_mode="prefix",
    )

    assert out == "\n".join(
        [
            "sample\tcovered\tpct_a_muts\tpct_c_muts\tpct_g_muts\tpct_u_muts",
            "sample1\t0\t12.5\t22.5\t32.5\t42.5",
            "",
        ]
    )
