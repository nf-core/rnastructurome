#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Download transcript FASTA sequences from NCBI Entrez by accession list or organism name."""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


EFETCH_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
MAX_RETRIES = 5
RETRY_BASE_DELAY = 10  # seconds; doubled on each retry


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download and concatenate FASTA sequences from NCBI Entrez."
    )
    parser.add_argument(
        "--accessions",
        default=None,
        help="Comma-separated NCBI nucleotide accession list (e.g. NC_002023.1,NC_002021.1). "
             "Takes precedence over --organism when both are provided.",
    )
    parser.add_argument(
        "--organism",
        required=True,
        help="Organism name (as in the samplesheet) used for logging and, when --accessions "
             "is absent, for an automatic NCBI esearch to find reference sequences.",
    )
    parser.add_argument("--output", required=True, help="Output gzip-compressed FASTA path")
    parser.add_argument(
        "--source-accessions",
        required=True,
        help="Output file listing one NCBI nuccore URL per accession",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API key (optional; increases rate limit from 3 to 10 req/s)",
    )
    parser.add_argument(
        "--efetch-url",
        default=EFETCH_URL,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--esearch-url",
        default=ESEARCH_URL,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--esummary-url",
        default=ESUMMARY_URL,
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def _get(url: str) -> bytes:
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504) and attempt < MAX_RETRIES - 1:
                delay = RETRY_BASE_DELAY * (2 ** attempt)
                print(
                    f"[NCBI_FASTA] HTTP {exc.code} on attempt {attempt + 1}/{MAX_RETRIES}. "
                    f"Retrying in {delay}s ...",
                    file=sys.stderr,
                )
                time.sleep(delay)
            else:
                raise
        except urllib.error.URLError as exc:
            if attempt < MAX_RETRIES - 1:
                delay = RETRY_BASE_DELAY * (2 ** attempt)
                print(
                    f"[NCBI_FASTA] URL error on attempt {attempt + 1}/{MAX_RETRIES}: {exc}. "
                    f"Retrying in {delay}s ...",
                    file=sys.stderr,
                )
                time.sleep(delay)
            else:
                raise
    raise RuntimeError("unreachable")  # pragma: no cover


def search_ncbi_accessions(
    organism: str, api_key: str | None,
    esearch_url: str = ESEARCH_URL, esummary_url: str = ESUMMARY_URL,
) -> list[str]:
    """Auto-search NCBI nuccore for RefSeq accessions matching the organism name.

    Tries an exact organism-name search first.  If that returns nothing it
    falls back to a broader term search so that minor name differences do not
    cause hard failures.  Prefers NC_ / NR_ RefSeq accessions when present.
    """
    def _esearch(term: str) -> list[str]:
        params: dict[str, str] = {
            "db": "nuccore",
            "term": term,
            "retmax": "100",
            "retmode": "json",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{esearch_url}?{urllib.parse.urlencode(params)}"
        data = json.loads(_get(url))
        return data["esearchresult"]["idlist"]

    # Exact organism name first, then relax to keyword search
    uid_list = _esearch(f'"{organism}"[Organism] AND RefSeq[Filter]')
    if not uid_list:
        uid_list = _esearch(f'{organism}[Organism] AND RefSeq[Filter]')
    if not uid_list:
        raise RuntimeError(
            f"[NCBI_FASTA] No RefSeq sequences found for organism '{organism}' on NCBI. "
            f"Check that the organism name matches NCBI taxonomy, or add explicit accessions "
            f"to params.ncbi_accessions_map in your pipeline config."
        )

    # Resolve UIDs → versioned accession strings via esummary
    id_str = ",".join(uid_list)
    params = {"db": "nuccore", "id": id_str, "retmode": "json"}
    if api_key:
        params["api_key"] = api_key  # type: ignore[assignment]
    summary = json.loads(_get(f"{esummary_url}?{urllib.parse.urlencode(params)}"))

    accessions: list[str] = []
    for uid in uid_list:
        acc = summary["result"].get(uid, {}).get("accessionversion", "")
        if acc:
            accessions.append(acc)

    if not accessions:
        raise RuntimeError(
            f"[NCBI_FASTA] Retrieved UIDs from NCBI but could not extract accession "
            f"numbers for organism '{organism}'."
        )

    # Prefer NC_ (genomic RefSeq) and NR_ (ribosomal RefSeq) when available
    canonical = [a for a in accessions if a.startswith(("NC_", "NR_"))]
    chosen = canonical if canonical else accessions

    print(
        f"[NCBI_FASTA] Auto-search for '{organism}' found {len(chosen)} accession(s): "
        f"{', '.join(chosen)}",
        file=sys.stderr,
    )
    return chosen


def fetch_fasta(accessions: list[str], api_key: str | None, efetch_url: str = EFETCH_URL) -> bytes:
    """Fetch all accessions as a single FASTA block from NCBI efetch."""
    id_str = ",".join(accessions)
    if efetch_url.startswith("file://"):
        url = efetch_url  # file:// handler ignores query params; used in tests
    else:
        params = f"db=nuccore&id={id_str}&rettype=fasta&retmode=text"
        if api_key:
            params += f"&api_key={api_key}"
        url = f"{efetch_url}?{params}"

    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                data = response.read()

            if not data.strip():
                raise RuntimeError("NCBI returned an empty response")

            # NCBI returns an HTML error page when accessions are not found
            if data.lstrip().startswith(b"<"):
                snippet = data[:300].decode("utf-8", errors="replace")
                raise RuntimeError(f"NCBI returned an HTML error page: {snippet!r}")

            return data

        except (urllib.error.URLError, RuntimeError) as exc:
            if attempt == MAX_RETRIES - 1:
                raise RuntimeError(
                    f"NCBI efetch failed after {MAX_RETRIES} attempts for "
                    f"{id_str!r}: {exc}"
                ) from exc
            delay = RETRY_BASE_DELAY * (2**attempt)
            print(
                f"[NCBI_FASTA] Attempt {attempt + 1}/{MAX_RETRIES} failed: {exc}. "
                f"Retrying in {delay}s ...",
                file=sys.stderr,
            )
            time.sleep(delay)

    raise RuntimeError("unreachable")  # pragma: no cover


def main() -> int:
    args = parse_args()

    if args.accessions:
        accessions = [a.strip() for a in args.accessions.split(",") if a.strip()]
        print(
            f"[NCBI_FASTA] Using {len(accessions)} pre-configured accession(s) "
            f"for '{args.organism}': {', '.join(accessions)}",
            file=sys.stderr,
        )
    else:
        print(
            f"[NCBI_FASTA] No accessions configured for '{args.organism}' — "
            f"searching NCBI automatically.",
            file=sys.stderr,
        )
        accessions = search_ncbi_accessions(
            args.organism, args.api_key,
            esearch_url=args.esearch_url, esummary_url=args.esummary_url,
        )

    if not accessions:
        print("[NCBI_FASTA] Error: no accessions resolved.", file=sys.stderr)
        return 1

    fasta_bytes = fetch_fasta(accessions, args.api_key, efetch_url=args.efetch_url)

    with gzip.open(args.output, "wb") as fh:
        fh.write(fasta_bytes)

    with open(args.source_accessions, "w", encoding="utf-8") as fh:
        for acc in accessions:
            fh.write(f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}\n")

    print(
        f"[NCBI_FASTA] Written {len(fasta_bytes):,} bytes compressed to {args.output!r}.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
