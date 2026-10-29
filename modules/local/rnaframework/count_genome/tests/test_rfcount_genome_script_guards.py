from pathlib import Path


MODULE = Path(__file__).resolve().parents[1] / "main.nf"


def test_rfcount_genome_validates_rc_files_in_output_directory():
    module = MODULE.read_text()

    assert 'def outdir = "${prefix}_rfcount_genome"' in module
    assert 'rc_count=\\$(find "${outdir}" -type f -name' in module
    assert 'rc_count=\\$(find "${prefix}" -type f -name' not in module


def test_rfcount_genome_parses_mutated_alignment_percentage_field():
    module = MODULE.read_text()

    assert 'if (\\$3 ~ /\\\\// && \\$4 ~ /^\\(/)' in module
    assert '\\$4 ~ /^${prefix}_rfcount_genome/' not in module


def test_rfcount_genome_tolerates_completed_nonzero_exit_only_after_output_validation():
    module = MODULE.read_text()

    assert 'pipeline_statuses=( "\\${PIPESTATUS[@]}" )' in module
    assert 'grep -Fq \'[+] All done.\'' in module
    assert 'rfcount_completed_with_nonzero=1' in module
    assert 'continuing because RC files were produced' in module
