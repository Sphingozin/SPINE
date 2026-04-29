"""
SPINE mutagenesis generator for In-Fusion cloning.

This is a separate file from the original SPINE.py and from the Golden Gate
alanine-scan copy. It does not use BsaI/BsmBI sites, barcode primers, or
Golden Gate overhangs.

Typical use:
    python SPINE_mutagenesis_infusion.py --fasta plasmid.fasta --gene-start 1800 --gene-end 3065 --mutation-regions 2050-2220,2500-2700 --output output --scan-mode alanine

Coordinates are 1-based nucleotide positions in the pasted/input plasmid.
Mutation regions are clipped to full codons inside the selected gene and split
into codon-aligned chunks that fit within the requested oligo length.
"""

import argparse
import csv
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt


SYNONYMOUS_CODONS = {
    'Cys': ['TGT', 'TGC'],
    'Asp': ['GAT', 'GAC'],
    'Ser': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'Gln': ['CAA', 'CAG'],
    'Met': ['ATG'],
    'Asn': ['AAC', 'AAT'],
    'Pro': ['CCT', 'CCG', 'CCA', 'CCC'],
    'Lys': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'Thr': ['ACC', 'ACA', 'ACG', 'ACT'],
    'Phe': ['TTT', 'TTC'],
    'Ala': ['GCA', 'GCC', 'GCG', 'GCT'],
    'Gly': ['GGT', 'GGG', 'GGA', 'GGC'],
    'Ile': ['ATC', 'ATA', 'ATT'],
    'Leu': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'His': ['CAT', 'CAC'],
    'Arg': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'Trp': ['TGG'],
    'Val': ['GTA', 'GTC', 'GTG', 'GTT'],
    'Glu': ['GAG', 'GAA'],
    'Tyr': ['TAT', 'TAC'],
}

HUMAN_USAGE = {
    'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.2,
    'CTT': 0.13, 'CTC': 0.2, 'CTA': 0.07, 'CTG': 0.41, 'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
    'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1, 'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
    'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
    'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52, 'TGG': 1,
    'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, 'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
    'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.2, 'AGG': 0.2,
    'GCT': 0.26, 'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25
}

ECOLI_USAGE = {
    'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41, 'TAA': 0.61, 'TAG': 0.09,
    'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
    'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
    'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
    'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54, 'TGA': 0.3, 'TGG': 1,
    'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
    'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04,
    'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15
}

MOUSE_USAGE = {
    'TTT': 0.44, 'TTC': 0.56, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.21,
    'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.08, 'CTG': 0.40, 'CAT': 0.42, 'CAC': 0.58, 'CAA': 0.26, 'CAG': 0.74,
    'ATT': 0.35, 'ATC': 0.50, 'ATA': 0.15, 'ATG': 1.0, 'AAT': 0.47, 'AAC': 0.53, 'AAA': 0.43, 'AAG': 0.57,
    'GTT': 0.18, 'GTC': 0.25, 'GTA': 0.11, 'GTG': 0.46, 'GAT': 0.47, 'GAC': 0.53, 'GAA': 0.43, 'GAG': 0.57,
    'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.46, 'TGC': 0.54, 'TGA': 0.51, 'TGG': 1.0,
    'CCT': 0.29, 'CCC': 0.31, 'CCA': 0.28, 'CCG': 0.12, 'CGT': 0.08, 'CGC': 0.18, 'CGA': 0.11, 'CGG': 0.20,
    'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.29, 'ACG': 0.11, 'AGT': 0.15, 'AGC': 0.25, 'AGA': 0.21, 'AGG': 0.22,
    'GCT': 0.27, 'GCC': 0.39, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.17, 'GGC': 0.33, 'GGA': 0.25, 'GGG': 0.25
}

AMINO_ACIDS = ['Cys', 'Asp', 'Ser', 'Gln', 'Met', 'Asn', 'Pro', 'Lys', 'Thr', 'Phe', 'Ala', 'Gly', 'Ile', 'Leu',
               'His', 'Arg', 'Trp', 'Val', 'Glu', 'Tyr']

CONSERVATIVE_MUTATIONS = {
    'Glu': ['Asp'],
    'Asp': ['Glu'],
    'Arg': ['Lys'],
    'Lys': ['Arg'],
    'Asn': ['Gln'],
    'Gln': ['Asn'],
    'Ser': ['Thr'],
    'Thr': ['Ser'],
    'Val': ['Ile', 'Leu'],
    'Ile': ['Val', 'Leu'],
    'Leu': ['Ile', 'Val'],
    'Phe': ['Tyr'],
    'Tyr': ['Phe'],
    'Cys': ['Ser'],
    'Met': ['Leu'],
    'His': ['Asn'],
    'Ala': ['Gly'],
    'Gly': ['Ala'],
    'Pro': ['Ala'],
    'Trp': ['Phe'],
}


def parse_ranges(text):
    ranges = []
    for item in [part.strip() for part in text.split(",") if part.strip()]:
        if "-" in item:
            start, end = [int(value.strip()) for value in item.split("-", 1)]
        else:
            start = int(item)
            end = start
        ranges.append((min(start, end), max(start, end)))
    return ranges


def circular_slice(sequence, start, length):
    start = start % len(sequence)
    doubled = sequence + sequence
    return doubled[start:start + length]


def circular_subsequence(sequence, start, end):
    start = start % len(sequence)
    end = end % len(sequence)
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]


def amino_acid_for_codon(codon):
    codon = str(codon).upper()
    for amino_acid, codons in SYNONYMOUS_CODONS.items():
        if codon in codons:
            return amino_acid
    raise ValueError("Unknown codon: " + codon)


def codon_usage_table(usage):
    if usage == "ecoli":
        return ECOLI_USAGE
    if usage == "mouse":
        return MOUSE_USAGE
    return HUMAN_USAGE


def preferred_codon(amino_acid, usage):
    usage_table = codon_usage_table(usage)
    return max(SYNONYMOUS_CODONS[amino_acid], key=lambda codon: usage_table.get(codon, 0))


def mutation_targets(wt_amino_acid, scan_mode):
    if scan_mode == "alanine":
        return [] if wt_amino_acid == 'Ala' else ['Ala']
    if scan_mode == "conservative":
        return CONSERVATIVE_MUTATIONS.get(wt_amino_acid, [])
    if scan_mode == "saturation":
        return [amino_acid for amino_acid in AMINO_ACIDS if amino_acid != wt_amino_acid]
    raise ValueError("Unknown scan mode: " + scan_mode)


def gc_content(seq):
    seq = str(seq).upper()
    if not seq:
        return 0.0
    return 100.0 * (seq.count('G') + seq.count('C')) / len(seq)


def primer_penalty(tm_value, gc_value, tm_min, tm_max, gc_min, gc_max):
    tm_target = (tm_min + tm_max) / 2
    tm_score = 0 if tm_min <= tm_value <= tm_max else abs(tm_value - tm_target)
    gc_score = 0 if gc_min <= gc_value <= gc_max else min(abs(gc_value - gc_min), abs(gc_value - gc_max))
    return tm_score + (gc_score / 10)


def primer_from_template(template, tm_min=56.0, tm_max=57.0, gc_min=40.0, gc_max=60.0, min_len=18, max_len=35):
    best = None
    for length in range(min_len, min(max_len, len(template)) + 1):
        primer = Seq(template[:length])
        tm_value = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
        gc_value = gc_content(primer)
        if tm_min <= tm_value <= tm_max and gc_min <= gc_value <= gc_max:
            return primer, round(tm_value, 1), round(gc_value, 1)
        penalty = primer_penalty(tm_value, gc_value, tm_min, tm_max, gc_min, gc_max)
        if best is None or penalty < best[3]:
            best = (primer, tm_value, gc_value, penalty)
    return best[0], round(best[1], 1), round(best[2], 1)


def codon_aligned_region(region_start, region_end, gene_start):
    region_start_zero = region_start - 1
    region_end_exclusive = region_end
    gene_start_zero = gene_start - 1
    codon_start = gene_start_zero + 3 * ((region_start_zero - gene_start_zero) // 3)
    codon_end = gene_start_zero + 3 * (((region_end_exclusive - gene_start_zero) + 2) // 3)
    return codon_start, codon_end


def split_region_into_chunks(region_start_zero, region_end_exclusive, max_insert_len, homology_len):
    max_mutable_len = max_insert_len - (2 * homology_len)
    if max_mutable_len < 3:
        raise ValueError("Oligo length is too short for the selected homology length.")
    max_mutable_len = 3 * (max_mutable_len // 3)
    chunks = []
    chunk_start = region_start_zero
    while chunk_start < region_end_exclusive:
        chunk_end = min(chunk_start + max_mutable_len, region_end_exclusive)
        chunk_end = chunk_start + 3 * ((chunk_end - chunk_start) // 3)
        if chunk_end <= chunk_start:
            raise ValueError("Could not create a codon-aligned In-Fusion chunk.")
        chunks.append((chunk_start, chunk_end))
        chunk_start = chunk_end
    return chunks


def design_region_vector_primers(record_id, sequence, region_start_zero, region_end_zero_exclusive, homology_len):
    downstream_template = circular_slice(sequence, region_end_zero_exclusive, 80)
    upstream_template = circular_slice(sequence, region_start_zero - 80, 80)
    forward, tm_forward, gc_forward = primer_from_template(downstream_template)
    reverse_template = Seq(upstream_template).reverse_complement()
    reverse, tm_reverse, gc_reverse = primer_from_template(reverse_template)
    region_name = str(region_start_zero + 1) + "-" + str(region_end_zero_exclusive)
    return [
        SeqRecord(
            forward,
            id=record_id + "_InFusion_Vector_" + region_name + "_F",
            description="Amplifies plasmid backbone after mutation region; Tm " + str(tm_forward) + "C; GC " + str(gc_forward) + "%"
        ),
        SeqRecord(
            reverse,
            id=record_id + "_InFusion_Vector_" + region_name + "_R",
            description="Amplifies plasmid backbone before mutation region; Tm " + str(tm_reverse) + "C; GC " + str(gc_reverse) + "%"
        ),
    ]


def generate_infusion_alanine_scan(
    fasta,
    gene_start,
    gene_end,
    mutation_regions,
    output,
    homology_len=15,
    oligo_len=230,
    usage="human",
    scan_mode="alanine",
):
    os.makedirs(output, exist_ok=True)
    records = list(SeqIO.parse(fasta, "fasta"))
    if not records:
        raise ValueError("No FASTA records found.")
    if (gene_end - gene_start + 1) % 3 != 0:
        raise ValueError("Gene start/end length must be divisible by 3.")

    inserts = []
    vector_primers = []
    summary_rows = []

    for record in records:
        plasmid = str(record.seq).upper()
        gene_start_zero = gene_start - 1
        gene_end_exclusive = gene_end
        if gene_start < 1 or gene_end > len(plasmid):
            raise ValueError(record.id + ": gene start/end is outside the plasmid sequence.")

        for input_start, input_end in mutation_regions:
            if input_start < gene_start or input_end > gene_end:
                raise ValueError("Mutation region " + str(input_start) + "-" + str(input_end) + " is outside the selected gene.")
            region_start_zero, region_end_exclusive = codon_aligned_region(input_start, input_end, gene_start)
            chunks = split_region_into_chunks(region_start_zero, region_end_exclusive, oligo_len, homology_len)
            print(
                record.id + " mutation region " + str(input_start) + "-" + str(input_end) +
                " split into " + str(len(chunks)) + " In-Fusion chunk(s)."
            )

            for chunk_start_zero, chunk_end_exclusive in chunks:
                region_seq = plasmid[chunk_start_zero:chunk_end_exclusive]
                left_homology = circular_slice(plasmid, chunk_start_zero - homology_len, homology_len)
                right_homology = circular_slice(plasmid, chunk_end_exclusive, homology_len)
                vector_primers.extend(
                    design_region_vector_primers(record.id, plasmid, chunk_start_zero, chunk_end_exclusive, homology_len)
                )

                for codon_start in range(chunk_start_zero, chunk_end_exclusive, 3):
                    codon = plasmid[codon_start:codon_start + 3]
                    wt_amino_acid = amino_acid_for_codon(codon)
                    aa_position = ((codon_start - gene_start_zero) // 3) + 1
                    for mutant_amino_acid in mutation_targets(wt_amino_acid, scan_mode):
                        mutant_codon = preferred_codon(mutant_amino_acid, usage)
                        mutation_offset = codon_start - chunk_start_zero
                        mutated_region = region_seq[:mutation_offset] + mutant_codon + region_seq[mutation_offset + 3:]
                        insert = left_homology + mutated_region + right_homology
                        insert_id = (
                            record.id + "_InFusion_" + scan_mode + "_nt" + str(chunk_start_zero + 1) + "-" +
                            str(chunk_end_exclusive) + "_" + wt_amino_acid + str(aa_position) + mutant_amino_acid
                        )
                        inserts.append(SeqRecord(Seq(insert), id=insert_id, description=""))
                        summary_rows.append([
                            record.id,
                            scan_mode,
                            input_start,
                            input_end,
                            region_start_zero + 1,
                            region_end_exclusive,
                            chunk_start_zero + 1,
                            chunk_end_exclusive,
                            aa_position,
                            wt_amino_acid,
                            codon,
                            mutant_amino_acid,
                            mutant_codon,
                            len(insert),
                        ])

    unique_primers = []
    seen_primers = set()
    for primer in vector_primers:
        if primer.id not in seen_primers:
            unique_primers.append(primer)
            seen_primers.add(primer.id)

    SeqIO.write(inserts, os.path.join(output, "InFusion_Mutagenesis_Inserts.fasta"), "fasta")
    SeqIO.write(unique_primers, os.path.join(output, "InFusion_Vector_Primers.fasta"), "fasta")
    with open(os.path.join(output, "InFusion_Mutagenesis_Summary.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "record_id",
            "scan_mode",
            "input_region_start_nt",
            "input_region_end_nt",
            "codon_aligned_region_start_nt",
            "codon_aligned_region_end_nt",
            "chunk_start_nt",
            "chunk_end_nt",
            "aa_position",
            "wt_amino_acid",
            "wt_codon",
            "mutant_amino_acid",
            "mutant_codon",
            "insert_length",
        ])
        writer.writerows(summary_rows)

    print("Generated " + str(len(inserts)) + " In-Fusion mutagenesis inserts.")
    print("Generated " + str(len(unique_primers)) + " vector primers.")
    print("Output folder: " + output)
    if not inserts:
        print("No inserts were generated. Check that mutation regions overlap mutable codons for the selected scan mode.")
    return inserts, unique_primers, summary_rows


def main():
    parser = argparse.ArgumentParser(description="Generate alanine-scan In-Fusion inserts from a full plasmid FASTA.")
    parser.add_argument("--fasta", required=True, help="Full plasmid FASTA file.")
    parser.add_argument("--gene-start", required=True, type=int, help="1-based nucleotide start of the gene in the plasmid.")
    parser.add_argument("--gene-end", required=True, type=int, help="1-based nucleotide end of the gene in the plasmid.")
    parser.add_argument("--mutation-regions", required=True, help="Comma-separated 1-based nucleotide regions, e.g. 2050-2220,2500-2700.")
    parser.add_argument("--output", required=True, help="Output folder.")
    parser.add_argument("--homology-len", default=15, type=int, help="In-Fusion homology length added to each insert end.")
    parser.add_argument("--oligo-len", default=230, type=int, help="Maximum total insert oligo length including homology arms.")
    parser.add_argument("--scan-mode", choices=["alanine", "conservative", "saturation"], default="alanine", help="Mutation scan type.")
    parser.add_argument("--usage", choices=["human", "ecoli", "mouse"], default="human", help="Codon usage for mutant codon choice.")
    args = parser.parse_args()

    generate_infusion_alanine_scan(
        fasta=args.fasta,
        gene_start=args.gene_start,
        gene_end=args.gene_end,
        mutation_regions=parse_ranges(args.mutation_regions),
        output=args.output,
        homology_len=args.homology_len,
        oligo_len=args.oligo_len,
        usage=args.usage,
        scan_mode=args.scan_mode,
    )


if __name__ == "__main__":
    main()
