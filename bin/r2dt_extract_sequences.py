#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""
Extract FASTA sequences from rf-fold dotbracket output for R2DT template search.

Usage
-----
    r2dt_extract_sequences.py --fold-dir <dir> --out <out.fa> \\
        [--gtf <ref.gtf>] [--allowed-biotypes rRNA,tRNA,...]

R2DT only has templates for structured ncRNA classes (rRNA, tRNA, snoRNA, snRNA,
SRP/RNase P, ...). Given a GTF and an allowlist of biotypes, only transcripts of
those biotypes are emitted, so R2DT never force-fits an mRNA/lncRNA to a wrong
template; everything else is left to ViennaRNA. With no --gtf, all sequences are
emitted (backward compatible). --fasta is accepted for backward compatibility.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

_TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
# Ensembl uses *_biotype; GENCODE uses *_type. Prefer transcript-level over gene-level.
_BIOTYPE_RES = {
    'transcript': re.compile(r'transcript_(?:biotype|type) "([^"]+)"'),
    'gene':       re.compile(r'gene_(?:biotype|type) "([^"]+)"'),
}
# Ensembl GTFs carry unversioned transcript_id (ENST…), but cDNA FASTA headers (and
# thus the folded .db IDs) are versioned (ENST….2). Match on both forms.
_VERSION_RE = re.compile(r'\.\d+$')


def _strip_version(tid):
    """Drop a trailing .<digits> Ensembl version suffix, if present."""
    return _VERSION_RE.sub('', tid)


def _open_text(path):
    """Open a plain or gzipped text file."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, encoding='utf-8')


def load_allowed_ids(gtf_path, allowed_biotypes):
    """Return the set of transcript IDs whose biotype is in the allowlist."""
    allowed = {b.strip().lower() for b in allowed_biotypes if b.strip()}
    allowed_ids = set()
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            tid_m = _TRANSCRIPT_ID_RE.search(line)
            if not tid_m:
                continue
            tbt = _BIOTYPE_RES['transcript'].search(line)
            gbt = _BIOTYPE_RES['gene'].search(line)
            biotype = (tbt.group(1) if tbt else (gbt.group(1) if gbt else None))
            if biotype and biotype.lower() in allowed:
                tid = tid_m.group(1)
                allowed_ids.add(tid)
                allowed_ids.add(_strip_version(tid))
    return allowed_ids


def main():
    """Extract FASTA sequences for R2DT input from dotbracket .db files."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fold-dir', required=True, type=Path,
                    help='Directory containing dotbracket/*.db files')
    ap.add_argument('--fasta',    required=False, type=Path,
                    help='Accepted for backward compatibility; not used')
    ap.add_argument('--gtf',      required=False, type=Path,
                    help='Reference GTF; restrict output to --allowed-biotypes transcripts')
    ap.add_argument('--allowed-biotypes', default='',
                    help='Comma-separated biotypes R2DT has templates for (with --gtf)')
    ap.add_argument('--out',      required=True, type=Path,
                    help='Output FASTA path')
    args = ap.parse_args()

    db_files = sorted((args.fold_dir / 'dotbracket').glob('*.db'))
    if not db_files:
        print(f'[R2DT] No .db files found in {args.fold_dir}/dotbracket', file=sys.stderr)
        args.out.write_text('')
        return

    # Restrict to transcripts R2DT can template (by GTF biotype), if a GTF is provided.
    allowed_ids = None
    if args.gtf and args.allowed_biotypes.strip():
        allowed_ids = load_allowed_ids(args.gtf, args.allowed_biotypes.split(','))
        print(
            f'[R2DT] {len(allowed_ids)} transcripts match allowed biotypes '
            f'({args.allowed_biotypes})',
            file=sys.stderr,
        )
        if not allowed_ids:
            print(
                '[R2DT] WARNING: no transcripts matched the biotype allowlist — '
                'R2DT will draw nothing; check the GTF has transcript/gene biotypes.',
                file=sys.stderr,
            )

    # Largest R2DT model is ~6000 nt; allow headroom above it and skip only clearly untemplatable
    # transcripts, so long LSU rRNAs still reach the template search.
    _max_len = 10000

    extracted     = []
    skipped_len   = 0
    skipped_type  = 0

    for db_path in db_files:
        lines = db_path.read_text(encoding='utf-8').splitlines()
        if len(lines) < 2 or not lines[0].startswith('>'):
            continue
        tid = lines[0][1:].split()[0]
        seq = lines[1].strip()
        if not seq:
            continue
        if allowed_ids is not None and tid not in allowed_ids \
                and _strip_version(tid) not in allowed_ids:
            skipped_type += 1
            continue
        if len(seq) <= _max_len:
            extracted.append((tid, seq))
        else:
            skipped_len += 1

    with open(args.out, 'w', encoding='utf-8') as fh:
        for h, s in extracted:
            fh.write(f'>{h}\n{s}\n')

    reasons = []
    if skipped_type:
        reasons.append(f'{skipped_type} not template-backed biotype')
    if skipped_len:
        reasons.append(f'{skipped_len} >{_max_len} nt')
    suffix = f' ({"; ".join(reasons)} skipped)' if reasons else ''
    print(
        f'[R2DT] Extracted {len(extracted)}/{len(db_files)} sequences for template search{suffix}',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
