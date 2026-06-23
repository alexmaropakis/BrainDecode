#!/usr/bin/env python3
"""
5_gnomad_parse.py

Annotate SAAP substitution sites against gnomAD v4.1 exomes (AF / SFS path).

For each SAAP (keyed on Ensembl transcript + protein position + encoded amino
acid), scan the gnomAD exome VCF over the gene region and read the embedded
Ensembl VEP CSQ. gnomAD's CSQ already provides, per transcript, the codon
(Codons), the amino-acid change (Amino_acids), the protein position, the strand,
and the transcript (Feature/ENST) -- so for any variant gnomAD observed there is
no need to reconstruct genomic coordinates. Allele count (AC) and frequency (AF)
are attached per matched alternate allele.

This is an ANNOTATION step: every input row is written with its gnomAD columns;
no row is filtered or dropped. Downstream subsetting (on AF, exact_mutation,
polymorphic, gnomad_status) is left to the analyst.

MATCH LOGIC
    A CSQ block matches a SAAP when, for the same alt allele:
        Feature (ENST, version-stripped) == SAAP transcript, AND
        Protein_position                 == SAAP protein_pos, AND
        encoded AA (Amino_acids ref side) == SAAP bp_aa
    Among matches, a record whose incorporated AA == mtp_aa (the gnomAD variant
    reproducing the SAME substitution as the SAAP) is preferred and reported as
    the polymorphism match. If no gnomAD variant reproduces the SAAP
    substitution but the site has some observed missense, the codon context is
    still reported (codon, strand, exact_mutation) WITHOUT borrowing the other
    variant's AF -- allele frequency is variant-specific.

INPUT
    --keep-csv   saap_mapped.csv (stage 4):
                   dataset, gene, transcript, protein, protein_pos,
                   bp_aa, mtp_aa, bp_seq, mtp_seq
    --gnomad-dir directory of gnomad.exomes.v4.1.sites.chr*.vcf.bgz (+ .tbi)
    --gene-bed   gene_coords.tsv (gene -> chrom/start/end), Ensembl-space contig
                   names (1,2,...,X,Y); chr-prefix is added at query time
    --out-dir    output directory
    --prefix     output filename prefix

OUTPUT
    {prefix}_gnomad_annotated.csv  (one row per input row)
"""

import argparse
import csv
import os
import pysam

# codon table to flag exact_mutation, not to call AAs
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

## Helper functions
def translate_codon(c):
    return CODON_TABLE.get(c.upper(), "X")

def snv_neighbors(codon):
    # yield (pos_in_codon[1..3], new_base, new_codon) for all 9 single-nt changes
    codon = codon.upper()
    for i in range(3):
        for b in "ACGT":
            if b == codon[i]:
                continue
            yield (i + 1, b, codon[:i] + b + codon[i + 1:])

def exact_mutation(ref_codon, bp_aa, mtp_aa):
    # True if a single-nt neighbour of ref_codon (which must encode bp_aa)
    # translates to mtp_aa; None if ref_codon does not encode bp_aa (guard).
    # Direction is BP -> MTP: ref codon encodes the base peptide residue,
    # neighbours that encode the incorporated residue are "exact".
    if translate_codon(ref_codon) != bp_aa:
        return None
    return any(translate_codon(nc) == mtp_aa for _, _, nc in snv_neighbors(ref_codon))

def parse_codons_field(codons_field):
    # VEP "Codons" e.g. "tCg/tGg" (changed base uppercased) ->
    #   (REF_CODON, ALT_CODON, changed_pos|None, ref_base|None, alt_base|None)
    # changed_pos is None for MNV/indel (more than one changed base).
    if "/" not in codons_field:
        return None
    rc, ac = codons_field.split("/")[0], codons_field.split("/")[1]
    if len(rc) != 3 or len(ac) != 3:
        return None
    changed = [i for i in range(3) if rc[i] != ac[i]]
    ref_codon, alt_codon = rc.upper(), ac.upper()
    if len(changed) != 1:
        return (ref_codon, alt_codon, None, None, None)
    i = changed[0]
    return (ref_codon, alt_codon, i + 1, rc[i].upper(), ac[i].upper())

def strip_ver(enst):
    return enst.split(".")[0] if enst else enst

def to_gnomad_contig(chrom):
    # gnomAD v4.1 exome VCFs use chr-prefixed contigs (verified: chr1).
    # Ensembl gene coords are bare (1); normalise to chr-prefixed for queries.
    return chrom if chrom.startswith("chr") else "chr" + chrom

def load_gene_coords(path):
    # gene_coords.tsv -> {GENE_UPPER: (chrom, start, end)}
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            n, c, s, e = parts[:4]
            d[n.upper()] = (c, int(s), int(e))
    return d

def vep_fields(vcf):
    # CSQ field order from the ##INFO=<ID=vep ...> header Description
    desc = vcf.header.info["vep"].description
    fmt = desc.split("Format:")[-1].strip()
    return fmt.split("|")

def passes(rec):
    # True for PASS or unfiltered records
    f = list(rec.filter.keys())
    return (not f) or ("PASS" in f)

def scan_gene(vcf, fields, chrom, start, end):
    """
    Scan a gene region and index observed missense variants.

    Returns (found, error) where found maps
        (ENST_noversion, protein_pos, encoded_aa, incorporated_aa)
    to a dict of {af, ac, varid, codons, strand, consequence,
                  chrom, pos_g, ref, alt}, keeping the max-AF record per key,
    and error is True if the tabix fetch failed (bad contig / index).
    """
    idx = {n: i for i, n in enumerate(fields)}
    found = {}

    try:
        recs = vcf.fetch(chrom, max(start - 1, 0), end)
    except (ValueError, OSError):
        return found, True

    for rec in recs:
        if not passes(rec):
            continue
        afs = rec.info.get("AF")
        acs = rec.info.get("AC")
        veps = rec.info.get("vep")
        if not veps:
            continue

        for ai, alt in enumerate(rec.alts):
            af = afs[ai] if afs and ai < len(afs) else None
            ac = acs[ai] if acs and ai < len(acs) else None

            for block in veps:
                p = block.split("|")

                def g(k):
                    j = idx.get(k)
                    return p[j] if j is not None and j < len(p) else ""

                if g("Allele") != alt:
                    continue
                if "missense_variant" not in g("Consequence"):
                    continue

                aa = g("Amino_acids")
                if "/" not in aa:
                    continue
                enc_aa, inc_aa = aa.split("/")[0], aa.split("/")[1]
                if len(enc_aa) != 1 or len(inc_aa) != 1:
                    continue

                pp = g("Protein_position").split("-")[0]
                if not pp.isdigit():
                    continue
                pos = int(pp)

                key = (strip_ver(g("Feature")), pos, enc_aa, inc_aa)
                cur = found.get(key)
                if cur is None or (af is not None and (cur["af"] is None or af > cur["af"])):
                    found[key] = {
                        "af": af, "ac": ac,
                        "varid": f"{rec.chrom}-{rec.pos}-{rec.ref}-{alt}",
                        "codons": g("Codons"), "strand": g("STRAND"),
                        "consequence": g("Consequence"),
                        "chrom": rec.chrom, "pos_g": rec.pos,
                        "ref": rec.ref, "alt": alt,
                    }
    return found, False

## Processing 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-csv", required=True, help="saap_mapped.csv")
    ap.add_argument("--gnomad-dir", required=True, help="dir of gnomad.exomes.v4.1.sites.chr*.vcf.bgq + .tbi files")
    ap.add_argument("--gene-bed", required=True, help="gene_coords.tsv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="Ping", help="likely dataset name, experimental dir")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    coords = load_gene_coords(a.gene_bed)

    with open(a.keep_csv, newline="") as f:
        rows = list(csv.DictReader(f))

    # group rows by gene for one VCF scan per gene
    by_gene = {}
    for r in rows:
        ok = (r.get("gene") and r.get("transcript")
              and r.get("protein_pos", "").isdigit()
              and len(r.get("bp_aa", "")) == 1 and len(r.get("mtp_aa", "")) == 1)
        if ok:
            by_gene.setdefault(r["gene"].upper(), []).append(r)
        else:
            r["_status"] = "unmappable_row"

    vcf_cache, fields_cache = {}, {}

    for gene, grp in by_gene.items():
        if gene not in coords:
            for r in grp:
                r["_status"] = "gene_not_found"
            continue

        chrom, start, end = coords[gene]
        gchrom = to_gnomad_contig(chrom)
        vcf_path = os.path.join(
            a.gnomad_dir, f"gnomad.exomes.v4.1.sites.{gchrom}.vcf.bgz")

        if vcf_path not in vcf_cache:
            if not os.path.exists(vcf_path):
                for r in grp:
                    r["_status"] = "vcf_missing"
                continue
            vcf_cache[vcf_path] = pysam.VariantFile(vcf_path)
            fields_cache[vcf_path] = vep_fields(vcf_cache[vcf_path])
        vcf, fields = vcf_cache[vcf_path], fields_cache[vcf_path]

        found, err = scan_gene(vcf, fields, gchrom, start, end)
        if err:
            for r in grp:
                r["_status"] = "fetch_error"
            continue

        for r in grp:
            enst = strip_ver(r["transcript"])
            pos = int(r["protein_pos"])
            bp, mtp = r["bp_aa"], r["mtp_aa"]

            # 1) preferred: gnomAD variant reproducing the SAAP substitution
            hit = found.get((enst, pos, bp, mtp))
            same_sub = hit is not None

            # 2) else: any observed missense at this site (encoded == bp), used
            #    only to recover the codon context -- NOT its allele frequency
            if hit is None:
                for (e, pp, ea, _ia), info in found.items():
                    if e == enst and pp == pos and ea == bp:
                        hit = info
                        break

            if hit is None:
                r["_status"] = "site_no_gnomad_missense"
                continue

            parsed = parse_codons_field(hit["codons"])
            if parsed:
                ref_codon, alt_codon, cpos, rb, ab = parsed
            else:
                ref_codon = alt_codon = ""
                cpos = rb = ab = None

            # codon context is valid for any observed missense at the site
            r["_codon"] = ref_codon
            r["_cpos"] = cpos
            r["_rb"] = rb
            r["_strand"] = hit["strand"]
            r["_exact"] = exact_mutation(ref_codon, bp, mtp) if ref_codon else None
            r["_same_sub"] = same_sub

            if same_sub:
                # variant-specific fields belong only to the matching variant
                r["_af"] = hit["af"]
                r["_ac"] = hit["ac"]
                r["_var"] = hit["varid"]
                r["_alt_codon"] = alt_codon
                r["_ab"] = ab
                r["_gchr"] = hit["chrom"]
                r["_gpos"] = hit["pos_g"]
                r["_gref"] = hit["ref"]
                r["_galt"] = hit["alt"]
                r["_conseq"] = hit["consequence"]
                r["_status"] = "reproduces_saap_sub"
            else:
                # codon context only; gnomAD has no AF for THIS substitution
                r["_af"] = None
                r["_ac"] = None
                r["_var"] = ""
                r["_alt_codon"] = ""
                r["_ab"] = None
                r["_gchr"] = r["_gpos"] = r["_gref"] = r["_galt"] = ""
                r["_conseq"] = hit["consequence"]
                r["_status"] = "site_variant_other_sub"

    # write files 
    base = [c for c in (rows[0].keys() if rows else []) if not c.startswith("_")]
    extra = ["gnomad_variant", "gnomad_chr", "gnomad_pos", "gnomad_ref",
             "gnomad_alt", "gnomad_ac", "gnomad_af", "polymorphic", "codon",
             "alt_codon", "codon_change_pos", "ref_base", "alt_base",
             "mutated_base", "strand", "consequence", "reproduces_saap_sub",
             "exact_mutation", "gnomad_status"]

    ann = os.path.join(a.out_dir, f"{a.prefix}_gnomad_annotated.csv")
    n_site = n_repro = n_poly = 0

    with open(ann, "w", newline="") as fa:
        wa = csv.DictWriter(fa, fieldnames=base + extra)
        wa.writeheader()

        for r in rows:
            st = r.get("_status", "site_no_gnomad_missense")
            has_site = st in ("reproduces_saap_sub", "site_variant_other_sub")
            repro = (st == "reproduces_saap_sub")
            af = r.get("_af") if repro else None
            ac = r.get("_ac") if repro else None

            if has_site:
                n_site += 1
                if repro:
                    n_repro += 1
                    poly = (af is not None and af > 0)
                    if poly:
                        n_poly += 1
                    poly_out = poly
                else:
                    poly_out = ""          # no variant-specific AF for this sub
            else:
                poly_out = ""

            out = {c: r.get(c, "") for c in base}
            out.update(
                gnomad_variant=r.get("_var", "") if repro else "",
                gnomad_chr=r.get("_gchr", "") if repro else "",
                gnomad_pos=r.get("_gpos", "") if repro else "",
                gnomad_ref=r.get("_gref", "") if repro else "",
                gnomad_alt=r.get("_galt", "") if repro else "",
                gnomad_ac=ac if (repro and ac is not None) else "",
                gnomad_af=af if (repro and af is not None) else "",
                polymorphic=poly_out,
                codon=r.get("_codon", "") if has_site else "",
                alt_codon=r.get("_alt_codon", "") if repro else "",
                codon_change_pos=(r.get("_cpos", "") if (repro and r.get("_cpos") is not None) else ""),
                ref_base=(r.get("_rb", "") if (has_site and r.get("_rb") is not None) else ""),
                alt_base=(r.get("_ab", "") if (repro and r.get("_ab") is not None) else ""),
                mutated_base=(r.get("_ab", "") if (repro and r.get("_ab") is not None) else ""),
                strand=r.get("_strand", "") if has_site else "",
                consequence=r.get("_conseq", "") if has_site else "",
                reproduces_saap_sub=r.get("_same_sub", "") if has_site else "",
                exact_mutation=("" if (not has_site or r.get("_exact") is None)
                                else r.get("_exact")),
                gnomad_status=st,
            )
            wa.writerow(out)

    print(f"{len(rows)} rows processed")
    print(f"{n_site} matched an observed gnomAD missense at the site")
    print(f"{n_repro} of those reproduce the SAAP substitution exactly")
    print(f"{n_poly} are polymorphic (AF>0)")
    print(f"wrote {ann}")


if __name__ == "__main__":
    main()
