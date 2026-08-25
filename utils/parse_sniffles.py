#!/usr/bin/env python3
"""
parse_sniffles.py -> Parse a Sniffles VCF and generate text coord files.

Port of parse_sniffles.pl (Matthew Field, Feb 2015).

Writes two tab-separated files:
    <output>.intrachr.tsv   intra-chromosomal events (del, dup, inv, ins)
    <output>.interchr.tsv   inter-chromosomal events (bnd / translocations)

Example:
    parse_sniffles.py --sniffles_vcf in.vcf --output sniffles
"""

import argparse
import os
import re
import sys


def die(msg):
    """Replacement for the local modules::Exception->throw call."""
    sys.exit(f"ERROR: {msg}")


def parse_sniffles_result(output_vcf):
    """Parse the Sniffles VCF into a nested dict keyed by
    sv_type -> "chr1:start1" -> "chr2:start2" -> {fields}."""
    sv_data = {}

    try:
        fh = open(output_vcf)
    except OSError:
        die(f"Can't open output file {output_vcf}")

    with fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue

            cols = line.split("\t")
            # chr1 start1 id ref_base var_base . . gff_str . evidence_str
            chr1 = cols[0]
            start1 = cols[1]
            id_ = cols[2]
            ref_base = cols[3]      # noqa: F841 (kept for parity with original)
            var_base = cols[4]
            gff_str = cols[7]
            evidence_str = cols[9]

            chr2 = chr1

            m = re.search(r"END=(\d+);", gff_str)
            start2 = m.group(1) if m else None

            m = re.search(r"SVLEN=-?(\d+)", gff_str)
            length = m.group(1) if m else None

            # zyg : (ignored) : ref_reads : var_reads
            ev = evidence_str.split(":")
            zyg = ev[0] if len(ev) > 0 else None
            ref_reads = ev[2] if len(ev) > 2 else None
            var_reads = ev[3] if len(ev) > 3 else None

            m = re.search(r"SVTYPE=([A-Z]+)", gff_str)
            event_type = m.group(1).lower() if m else None

            m = re.search(r"VAF=([0-9.]+)", gff_str)
            vaf = m.group(1) if m else None

            if event_type == "bnd":
                m = re.search(r"[\[\]](.+):(\d+)[\[\]]", var_base)
                if m:
                    chr2, start2 = m.group(1), m.group(2)
                else:
                    die(f"with {var_base}")

                if not re.search(r"\S", chr2 or ""):
                    die(f"with {var_base}")

                # Don't record SVs twice
                if "bnd" in sv_data and f"{chr2}:{start2}" in sv_data["bnd"]:
                    continue

            key1 = f"{chr1}:{start1}"
            key2 = f"{chr2}:{start2}"

            entry = (
                sv_data
                .setdefault(event_type, {})
                .setdefault(key1, {})
                .setdefault(key2, {})
            )
            entry["ref_reads"] = ref_reads
            entry["var_reads"] = var_reads
            entry["length"] = length
            entry["gff"] = gff_str
            entry["id"] = id_
            entry["zyg"] = zyg
            entry["vaf"] = vaf

            if event_type == "ins":
                if len(var_base) > 1000:
                    var_base = ">1000bp"
                entry["seq"] = var_base

    return sv_data


def _chr_coord_sort_key(coord_str):
    """Sort by chromosome (string) then coordinate (numeric),
    matching the Perl 'a_chr cmp b_chr || a_coord <=> b_coord'."""
    m = re.search(r"(\S+):(\d+)", coord_str)
    if m:
        return (m.group(1), int(m.group(2)))
    return (coord_str, 0)


def print_files(results, file_base):
    intrachr_path = f"{file_base}.intrachr.tsv"
    interchr_path = f"{file_base}.interchr.tsv"

    try:
        intra = open(intrachr_path, "w")
    except OSError:
        die(f"Can't open file {intrachr_path}")
    try:
        tra = open(interchr_path, "w")
    except OSError:
        die(f"Can't open file {interchr_path}")

    with intra, tra:
        # intra-chromosomal header
        intra.write("\t".join([
            "CTG", "START_COORD", "END_COORD", "SV_TYPE", "ID",
            "VAR_READS", "REF_READS", "LEN", "ZYG", "COORD",
            "VAF", "INSERT_BP", "GFF_FULL",
        ]) + "\n")

        # inter-chromosomal header
        # (in the original Perl this was written to the intra file by mistake;
        #  here it correctly goes to the interchr file)
        tra.write("\t".join([
            "CTG1", "START_COORD", "CTG2", "END_COORD", "SV_TYPE", "ID",
            "VAR_READS", "REF_READS", "LEN", "ZYG", "COORD",
            "VAF", "INSERT_BP", "GFF_FULL",
        ]) + "\n")

        for sv_type in ("del", "dup", "inv", "ins", "bnd"):
            type_data = results.get(sv_type, {})
            for coord_str1 in sorted(type_data, key=_chr_coord_sort_key):
                for coord_str2 in type_data[coord_str1]:
                    rec = type_data[coord_str1][coord_str2]

                    chr1, coord1 = coord_str1.split(":", 1)
                    chr2, coord2 = coord_str2.split(":", 1)

                    ins = rec.get("seq", "N/A") if sv_type == "ins" else "N/A"
                    coord_combined = f"{coord_str1}+{coord_str2}"

                    if sv_type == "bnd":
                        tra.write("\t".join(str(x) for x in [
                            chr1, coord1, chr2, coord2, sv_type,
                            rec.get("id"), rec.get("var_reads"),
                            rec.get("ref_reads"), rec.get("length"),
                            rec.get("zyg"), coord_combined, rec.get("vaf"),
                            ins, rec.get("gff"),
                        ]) + "\n")
                    else:
                        intra.write("\t".join(str(x) for x in [
                            chr1, coord1, coord2, sv_type,
                            rec.get("id"), rec.get("var_reads"),
                            rec.get("ref_reads"), rec.get("length"),
                            rec.get("zyg"), coord_combined, rec.get("vaf"),
                            ins, rec.get("gff"),
                        ]) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Parse a Sniffles VCF and generate text coord files."
    )
    parser.add_argument("--sniffles_vcf", required=True,
                        help="Input Sniffles VCF")
    parser.add_argument("--output", default="sniffles",
                        help="Output file basename (default: sniffles)")
    args = parser.parse_args()

    if not os.path.exists(args.sniffles_vcf):
        die(f"File {args.sniffles_vcf} doesn't exist")

    sniffles_sv = parse_sniffles_result(args.sniffles_vcf)
    print_files(sniffles_sv, "./" + args.output)


if __name__ == "__main__":
    main()