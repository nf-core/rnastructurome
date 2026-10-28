import subprocess
import sys
from pathlib import Path


def test_remap_wig_to_genome_ucsc_names_and_overlap_mean(tmp_path):
    repo_root = Path(__file__).resolve().parents[4]
    script = repo_root / "bin" / "remap_wig_to_genome.py"
    gtf = repo_root / "modules" / "local" / "wig_to_genome" / "tests" / "fixtures" / "annotation.gtf"
    wig = repo_root / "modules" / "local" / "wig_to_genome" / "tests" / "fixtures" / "transcript.wig"
    output_wig = tmp_path / "genomic.wig"
    chrom_sizes = tmp_path / "genomic.chrom.sizes"

    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--wig",
            str(wig),
            "--gtf",
            str(gtf),
            "--output-wig",
            str(output_wig),
            "--chrom-sizes",
            str(chrom_sizes),
            "--organism",
            "Homo sapiens",
            "--ucsc-common-chrom-names",
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert output_wig.read_text().splitlines() == [
        "track type=wiggle_0",
        "variableStep chrom=chr1",
        "101 3",
        "102 3",
        "variableStep chrom=chr2",
        "101 40",
        "201 20",
        "202 10",
    ]
    assert chrom_sizes.read_text().splitlines() == [
        "chr1\t202",
        "chr2\t202",
    ]
