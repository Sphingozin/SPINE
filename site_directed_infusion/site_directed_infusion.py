"""
Site-directed In-Fusion saturation mutagenesis primer designer.

Input:
    - full plasmid FASTA
    - gene start/end coordinates in the plasmid
    - exact amino-acid position to mutate
    - two restriction enzymes used to digest the receiving backbone

Output:
    - primer FASTA
    - primer CSV summary

Primer layout for each mutation:
    1. F1: restriction site 1 + gene-start binding sequence
    2. R1: reverse primer ending immediately before the mutation codon
    3. F2: overlap to fragment 1 + mutant codon + downstream binding sequence
    4. R2: restriction site 2 + gene-end reverse binding sequence

Two designs are supported:
    - full saturation: one F2 primer per target amino acid
    - degenerate: compact F2 primer set using NNT, VAA, ATG, and TGG
"""

import argparse
import csv
import os
import queue
import threading
import traceback
import tkinter as tk
from datetime import datetime
from tkinter import filedialog, messagebox, scrolledtext, ttk

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

AMINO_ACIDS = ['Cys', 'Asp', 'Ser', 'Gln', 'Met', 'Asn', 'Pro', 'Lys', 'Thr', 'Phe',
               'Ala', 'Gly', 'Ile', 'Leu', 'His', 'Arg', 'Trp', 'Val', 'Glu', 'Tyr']

HUMAN_USAGE = {
    'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57,
    'CTT': 0.13, 'CTC': 0.2, 'CTA': 0.07, 'CTG': 0.41, 'CAT': 0.41, 'CAC': 0.59,
    'CAA': 0.25, 'CAG': 0.75, 'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1,
    'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58, 'GTT': 0.18, 'GTC': 0.24,
    'GTA': 0.11, 'GTG': 0.47, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
    'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.45, 'TGC': 0.55,
    'TGG': 1, 'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, 'CGT': 0.08,
    'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28,
    'ACG': 0.12, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.2, 'AGG': 0.2, 'GCT': 0.26,
    'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25,
    'GGG': 0.25
}

MOUSE_USAGE = {
    **HUMAN_USAGE,
    'GCC': 0.39, 'GCT': 0.27, 'GCA': 0.23, 'GCG': 0.11,
    'GAA': 0.43, 'GAG': 0.57,
}

ECOLI_USAGE = {
    'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41,
    'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43,
    'CAA': 0.34, 'CAG': 0.66, 'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1,
    'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26, 'GTT': 0.28, 'GTC': 0.2,
    'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
    'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54,
    'TGG': 1, 'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49, 'CGT': 0.36,
    'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17,
    'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04, 'GCT': 0.18,
    'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13,
    'GGG': 0.15
}

ENZYMES = {
    'AgeI': 'ACCGGT',
    'BamHI': 'GGATCC',
    'BglII': 'AGATCT',
    'EcoRI': 'GAATTC',
    'HindIII': 'AAGCTT',
    'KpnI': 'GGTACC',
    'NcoI': 'CCATGG',
    'NdeI': 'CATATG',
    'NotI': 'GCGGCCGC',
    'NheI': 'GCTAGC',
    'SalI': 'GTCGAC',
    'SpeI': 'ACTAGT',
    'XbaI': 'TCTAGA',
    'XhoI': 'CTCGAG',
}

DEGENERATE_SETS = [
    ('NNT', 'NNT mixed codons'),
    ('VAA', 'VAA mixed codons'),
    ('ATG', 'Met codon'),
    ('TGG', 'Trp codon'),
]


def parse_fasta_text_or_file(fasta_path=None, fasta_text=None):
    if fasta_text:
        records = list(SeqIO.parse(fasta_text.splitlines(), "fasta"))
    else:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError("No FASTA sequence found.")
    return records[0]


def usage_table(name):
    if name == 'mouse':
        return MOUSE_USAGE
    if name == 'ecoli':
        return ECOLI_USAGE
    return HUMAN_USAGE


def aa_for_codon(codon):
    codon = str(codon).upper()
    for aa, codons in SYNONYMOUS_CODONS.items():
        if codon in codons:
            return aa
    raise ValueError("Unknown codon: " + codon)


def preferred_codon(aa, usage):
    table = usage_table(usage)
    return max(SYNONYMOUS_CODONS[aa], key=lambda codon: table.get(codon, 0))


def gc_content(seq):
    seq = str(seq).upper()
    return 100.0 * (seq.count('G') + seq.count('C')) / len(seq) if seq else 0.0


def choose_binding_primer(template, tm_min=56.0, tm_max=57.0, gc_min=40.0, gc_max=60.0, min_len=18, max_len=35):
    best = None
    target_tm = (tm_min + tm_max) / 2
    if not template:
        raise ValueError("No template sequence is available for primer design.")
    start_len = min(min_len, len(template))
    for length in range(start_len, min(max_len, len(template)) + 1):
        primer = Seq(str(template[:length]).upper())
        tm_value = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
        gc_value = gc_content(primer)
        tm_penalty = 0 if tm_min <= tm_value <= tm_max else abs(tm_value - target_tm)
        gc_penalty = 0 if gc_min <= gc_value <= gc_max else min(abs(gc_value - gc_min), abs(gc_value - gc_max)) / 10
        penalty = tm_penalty + gc_penalty
        if tm_min <= tm_value <= tm_max and gc_min <= gc_value <= gc_max:
            return primer, round(tm_value, 1), round(gc_value, 1)
        if best is None or penalty < best[3]:
            best = (primer, tm_value, gc_value, penalty)
    return best[0], round(best[1], 1), round(best[2], 1)


def enzyme_sequence(enzyme):
    enzyme = enzyme.strip()
    if enzyme in ENZYMES:
        return ENZYMES[enzyme]
    if all(base.upper() in "ACGT" for base in enzyme):
        return enzyme.upper()
    raise ValueError("Unknown enzyme or invalid sequence: " + enzyme)


def make_run_folder(output_parent):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(output_parent, "SiteDirected_InFusion_" + stamp)
    suffix = 1
    while os.path.exists(path):
        path = os.path.join(output_parent, "SiteDirected_InFusion_" + stamp + "_" + str(suffix))
        suffix += 1
    os.makedirs(path)
    return path


def primer_record(seq, primer_id, description):
    return SeqRecord(Seq(str(seq).upper()), id=primer_id, description=description)


def design_site_directed_infusion(
    fasta_path=None,
    fasta_text=None,
    gene_start=None,
    gene_end=None,
    aa_position=None,
    enzyme1='AgeI',
    enzyme2='NotI',
    design_mode='saturation',
    usage='human',
    overlap_len=18,
    output='.',
):
    record = parse_fasta_text_or_file(fasta_path=fasta_path, fasta_text=fasta_text)
    plasmid = str(record.seq).upper()
    gene_start_zero = gene_start - 1
    gene_end_exclusive = gene_end
    if gene_start < 1 or gene_end > len(plasmid):
        raise ValueError("Gene coordinates are outside the plasmid sequence.")
    if (gene_end - gene_start + 1) % 3 != 0:
        raise ValueError("Gene length must be divisible by 3.")

    codon_start = gene_start_zero + (aa_position - 1) * 3
    codon_end = codon_start + 3
    if codon_end > gene_end_exclusive:
        raise ValueError("Amino-acid position is outside the selected gene.")

    wt_codon = plasmid[codon_start:codon_end]
    wt_aa = aa_for_codon(wt_codon)
    enzyme1_seq = enzyme_sequence(enzyme1)
    enzyme2_seq = enzyme_sequence(enzyme2)

    run_folder = make_run_folder(output)
    primers = []
    rows = []

    gene_start_template = plasmid[gene_start_zero:gene_start_zero + 80]
    f1_binding, f1_tm, f1_gc = choose_binding_primer(gene_start_template)
    f1 = enzyme1_seq + str(f1_binding)
    primers.append(primer_record(
        f1,
        record.id + "_F1_" + enzyme1,
        "restriction_site=" + enzyme1_seq + "; binding_Tm=" + str(f1_tm) + "C; binding_GC=" + str(f1_gc) + "%"
    ))
    rows.append(["F1", "constant", "", "", f1, f1_tm, f1_gc, enzyme1_seq])

    r1_template_start = max(gene_start_zero, codon_start - 80)
    r1_template = plasmid[r1_template_start:codon_start]
    r1_binding, r1_tm, r1_gc = choose_binding_primer(str(Seq(r1_template).reverse_complement()))
    primers.append(primer_record(
        r1_binding,
        record.id + "_R1_before_" + wt_aa + str(aa_position),
        "ends_before_mutation; binding_Tm=" + str(r1_tm) + "C; binding_GC=" + str(r1_gc) + "%"
    ))
    rows.append(["R1", "constant", "", "", str(r1_binding), r1_tm, r1_gc, ""])

    gene_end_template = plasmid[max(gene_start_zero, gene_end_exclusive - 80):gene_end_exclusive]
    r2_binding, r2_tm, r2_gc = choose_binding_primer(str(Seq(gene_end_template).reverse_complement()))
    r2 = enzyme2_seq + str(r2_binding)
    primers.append(primer_record(
        r2,
        record.id + "_R2_" + enzyme2,
        "restriction_site=" + enzyme2_seq + "; binding_Tm=" + str(r2_tm) + "C; binding_GC=" + str(r2_gc) + "%"
    ))
    rows.append(["R2", "constant", "", "", r2, r2_tm, r2_gc, enzyme2_seq])

    overlap = plasmid[max(gene_start_zero, codon_start - overlap_len):codon_start]
    downstream_template = plasmid[codon_end:codon_end + 80]
    f2_binding, f2_tm, f2_gc = choose_binding_primer(downstream_template)

    if design_mode == "degenerate":
        targets = DEGENERATE_SETS
    else:
        targets = [(aa, preferred_codon(aa, usage)) for aa in AMINO_ACIDS if aa != wt_aa]

    for target, codon in targets:
        f2 = overlap + codon + str(f2_binding)
        primer_id = record.id + "_F2_" + wt_aa + str(aa_position) + target
        primers.append(primer_record(
            f2,
            primer_id,
            "mutant_codon=" + codon + "; downstream_binding_Tm=" + str(f2_tm) + "C; downstream_binding_GC=" + str(f2_gc) + "%"
        ))
        rows.append(["F2", design_mode, target, codon, f2, f2_tm, f2_gc, ""])

    SeqIO.write(primers, os.path.join(run_folder, "SiteDirected_InFusion_Primers.fasta"), "fasta")
    with open(os.path.join(run_folder, "SiteDirected_InFusion_Primers.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["primer", "mode", "target", "mutant_codon", "sequence", "binding_tm_c", "binding_gc_percent", "restriction_site"])
        writer.writerows(rows)

    with open(os.path.join(run_folder, "SiteDirected_InFusion_Design.txt"), "w") as handle:
        handle.write("Record: " + record.id + "\n")
        handle.write("Gene coordinates: " + str(gene_start) + "-" + str(gene_end) + "\n")
        handle.write("Mutation: " + wt_aa + str(aa_position) + "X; WT codon " + wt_codon + "\n")
        handle.write("Mode: " + design_mode + "\n")
        handle.write("Restriction enzymes: " + enzyme1 + "=" + enzyme1_seq + ", " + enzyme2 + "=" + enzyme2_seq + "\n")
        handle.write("Overlap length before mutation in F2 primers: " + str(len(overlap)) + " nt\n")
        handle.write("Output folder: " + run_folder + "\n")

    print("Designed " + str(len(primers)) + " primers.")
    print("WT mutation site: " + wt_aa + str(aa_position) + " (" + wt_codon + ")")
    print("Output folder: " + run_folder)
    return run_folder


class SiteDirectedGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Site-Directed In-Fusion Mutagenesis")
        self.geometry("920x700")
        self.output_queue = queue.Queue()
        self.worker = None
        self.build_ui()
        self.after(100, self.pump_output)

    def build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        form = ttk.Frame(self, padding=14)
        form.grid(row=0, column=0, sticky="ew")
        form.columnconfigure(1, weight=1)

        self.mode_var = tk.StringVar(value="file")
        self.fasta_var = tk.StringVar()
        self.output_var = tk.StringVar(value=os.getcwd())
        self.start_var = tk.StringVar()
        self.end_var = tk.StringVar()
        self.aa_var = tk.StringVar()
        self.enzyme1_var = tk.StringVar(value="AgeI")
        self.enzyme2_var = tk.StringVar(value="NotI")
        self.design_var = tk.StringVar(value="saturation")
        self.usage_var = tk.StringVar(value="human")
        self.overlap_var = tk.StringVar(value="18")

        mode_frame = ttk.Frame(form)
        mode_frame.grid(row=0, column=1, sticky="w", pady=4)
        ttk.Radiobutton(mode_frame, text="Use FASTA file", variable=self.mode_var, value="file", command=self.toggle_input).grid(row=0, column=0, padx=(0, 16))
        ttk.Radiobutton(mode_frame, text="Paste FASTA", variable=self.mode_var, value="paste", command=self.toggle_input).grid(row=0, column=1)

        ttk.Label(form, text="Full plasmid FASTA").grid(row=1, column=0, sticky="w", pady=4)
        self.fasta_entry = ttk.Entry(form, textvariable=self.fasta_var)
        self.fasta_entry.grid(row=1, column=1, sticky="ew", pady=4)
        self.fasta_button = ttk.Button(form, text="Browse", command=self.choose_fasta)
        self.fasta_button.grid(row=1, column=2, padx=(8, 0), pady=4)

        ttk.Label(form, text="Paste FASTA").grid(row=2, column=0, sticky="nw", pady=4)
        self.pasted_fasta = scrolledtext.ScrolledText(form, height=7, wrap="word")
        self.pasted_fasta.grid(row=2, column=1, columnspan=2, sticky="ew", pady=4)
        self.pasted_fasta.insert("1.0", ">my_plasmid\n")

        self.row(form, 3, "Output folder", self.output_var, browse=self.choose_output)
        self.row(form, 4, "Gene start", self.start_var, "1-based plasmid nt")
        self.row(form, 5, "Gene end", self.end_var, "1-based plasmid nt")
        self.row(form, 6, "AA position", self.aa_var, "exact residue number in selected gene")

        ttk.Label(form, text="Enzyme 1").grid(row=7, column=0, sticky="w", pady=4)
        ttk.Combobox(form, textvariable=self.enzyme1_var, values=tuple(ENZYMES.keys()), width=14).grid(row=7, column=1, sticky="w", pady=4)
        ttk.Label(form, text="Enzyme 2").grid(row=8, column=0, sticky="w", pady=4)
        ttk.Combobox(form, textvariable=self.enzyme2_var, values=tuple(ENZYMES.keys()), width=14).grid(row=8, column=1, sticky="w", pady=4)

        ttk.Label(form, text="Design mode").grid(row=9, column=0, sticky="w", pady=4)
        ttk.Combobox(form, textvariable=self.design_var, values=("saturation", "degenerate"), state="readonly", width=14).grid(row=9, column=1, sticky="w", pady=4)
        ttk.Label(form, text="Codon usage").grid(row=10, column=0, sticky="w", pady=4)
        ttk.Combobox(form, textvariable=self.usage_var, values=("human", "mouse", "ecoli"), state="readonly", width=14).grid(row=10, column=1, sticky="w", pady=4)
        self.row(form, 11, "F2 overlap", self.overlap_var, "nt upstream of mutated codon")

        run = ttk.Button(form, text="Design primers", command=self.start_run)
        run.grid(row=12, column=2, sticky="e", pady=(10, 0))
        self.run_button = run

        out = ttk.Frame(self, padding=(14, 0, 14, 14))
        out.grid(row=1, column=0, sticky="nsew")
        out.columnconfigure(0, weight=1)
        out.rowconfigure(0, weight=1)
        self.log = tk.Text(out, wrap="word")
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(out, orient="vertical", command=self.log.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log.configure(yscrollcommand=scroll.set)
        self.toggle_input()

    def row(self, parent, row, label, var, hint="", browse=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(parent, textvariable=var).grid(row=row, column=1, sticky="ew", pady=4)
        if browse:
            ttk.Button(parent, text="Browse", command=browse).grid(row=row, column=2, padx=(8, 0), pady=4)
        elif hint:
            ttk.Label(parent, text=hint).grid(row=row, column=2, sticky="w", padx=(8, 0), pady=4)

    def choose_fasta(self):
        path = filedialog.askopenfilename(filetypes=(("FASTA", "*.fasta *.fa *.fas"), ("All files", "*.*")))
        if path:
            self.fasta_var.set(path)

    def choose_output(self):
        path = filedialog.askdirectory()
        if path:
            self.output_var.set(path)

    def toggle_input(self):
        use_file = self.mode_var.get() == "file"
        self.fasta_entry.configure(state="normal" if use_file else "disabled")
        self.fasta_button.configure(state="normal" if use_file else "disabled")
        self.pasted_fasta.configure(state="disabled" if use_file else "normal")

    def read_config(self):
        fasta_text = None
        fasta_path = None
        if self.mode_var.get() == "paste":
            fasta_text = self.pasted_fasta.get("1.0", "end").strip()
            if not fasta_text.startswith(">"):
                raise ValueError("Pasted FASTA must start with a > header.")
        else:
            fasta_path = self.fasta_var.get().strip()
            if not os.path.isfile(fasta_path):
                raise ValueError("Choose a FASTA file.")
        output = self.output_var.get().strip()
        os.makedirs(output, exist_ok=True)
        return {
            "fasta_path": fasta_path,
            "fasta_text": fasta_text,
            "gene_start": int(self.start_var.get().strip()),
            "gene_end": int(self.end_var.get().strip()),
            "aa_position": int(self.aa_var.get().strip()),
            "enzyme1": self.enzyme1_var.get().strip(),
            "enzyme2": self.enzyme2_var.get().strip(),
            "design_mode": self.design_var.get(),
            "usage": self.usage_var.get(),
            "overlap_len": int(self.overlap_var.get().strip()),
            "output": output,
        }

    def start_run(self):
        try:
            config = self.read_config()
        except Exception as error:
            messagebox.showerror("Check inputs", str(error))
            return
        self.log.delete("1.0", "end")
        self.run_button.configure(state="disabled")
        self.worker = threading.Thread(target=self.run_design, args=(config,), daemon=True)
        self.worker.start()

    def run_design(self, config):
        try:
            import contextlib
            import io
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                design_site_directed_infusion(**config)
            self.output_queue.put(("done", buffer.getvalue()))
        except Exception:
            self.output_queue.put(("error", traceback.format_exc()))

    def pump_output(self):
        try:
            while True:
                kind, text = self.output_queue.get_nowait()
                self.log.insert("end", text)
                self.log.see("end")
                self.run_button.configure(state="normal")
                if kind == "error":
                    messagebox.showerror("Run failed", "See log for details.")
                else:
                    messagebox.showinfo("Finished", "Primer design complete.")
        except queue.Empty:
            pass
        self.after(100, self.pump_output)


def main():
    parser = argparse.ArgumentParser(description="Design site-directed In-Fusion mutagenesis primers.")
    parser.add_argument("--gui", action="store_true", help="Start the graphical interface.")
    parser.add_argument("--fasta", help="Full plasmid FASTA file.")
    parser.add_argument("--gene-start", type=int)
    parser.add_argument("--gene-end", type=int)
    parser.add_argument("--aa-position", type=int)
    parser.add_argument("--enzyme1", default="AgeI")
    parser.add_argument("--enzyme2", default="NotI")
    parser.add_argument("--design-mode", choices=["saturation", "degenerate"], default="saturation")
    parser.add_argument("--usage", choices=["human", "mouse", "ecoli"], default="human")
    parser.add_argument("--overlap-len", type=int, default=18)
    parser.add_argument("--output", default=".")
    args = parser.parse_args()

    if args.gui:
        SiteDirectedGui().mainloop()
        return

    design_site_directed_infusion(
        fasta_path=args.fasta,
        gene_start=args.gene_start,
        gene_end=args.gene_end,
        aa_position=args.aa_position,
        enzyme1=args.enzyme1,
        enzyme2=args.enzyme2,
        design_mode=args.design_mode,
        usage=args.usage,
        overlap_len=args.overlap_len,
        output=args.output,
    )


if __name__ == "__main__":
    main()
