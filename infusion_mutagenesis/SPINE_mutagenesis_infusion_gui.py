import contextlib
import io
import os
import queue
import sys
import threading
import traceback
import tkinter as tk
from datetime import datetime
from tkinter import filedialog, messagebox, scrolledtext, ttk

from SPINE_mutagenesis_infusion import generate_infusion_alanine_scan, parse_ranges


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))


def create_run_folder(parent_folder, prefix):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_folder = os.path.join(parent_folder, prefix + "_" + timestamp)
    suffix = 1
    while os.path.exists(run_folder):
        run_folder = os.path.join(parent_folder, prefix + "_" + timestamp + "_" + str(suffix))
        suffix += 1
    os.makedirs(run_folder)
    return run_folder


class QueueWriter(io.TextIOBase):
    def __init__(self, output_queue):
        self.output_queue = output_queue

    def write(self, text):
        if text:
            self.output_queue.put(("log", text))
        return len(text)

    def flush(self):
        return None


class SpineInfusionGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SPINE Mutagenesis - In-Fusion")
        self.geometry("920x680")
        self.minsize(820, 580)
        self.output_queue = queue.Queue()
        self.worker = None
        self._build_ui()
        self.after(100, self._pump_output)

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        form = ttk.Frame(self, padding=14)
        form.grid(row=0, column=0, sticky="ew")
        form.columnconfigure(1, weight=1)

        self.input_mode_var = tk.StringVar(value="file")
        self.fasta_var = tk.StringVar()
        self.output_var = tk.StringVar(value=SCRIPT_DIR)
        self.start_var = tk.StringVar()
        self.end_var = tk.StringVar()
        self.regions_var = tk.StringVar()
        self.homology_var = tk.StringVar(value="18")
        self.oligo_len_var = tk.StringVar(value="230")
        self.usage_var = tk.StringVar(value="human")
        self.scan_mode_var = tk.StringVar(value="alanine")

        input_mode = ttk.Frame(form)
        input_mode.grid(row=0, column=1, columnspan=2, sticky="w", pady=4)
        ttk.Radiobutton(input_mode, text="Use FASTA file", variable=self.input_mode_var, value="file", command=self._toggle_input_mode).grid(row=0, column=0, padx=(0, 18))
        ttk.Radiobutton(input_mode, text="Paste FASTA", variable=self.input_mode_var, value="paste", command=self._toggle_input_mode).grid(row=0, column=1)

        ttk.Label(form, text="Full plasmid FASTA").grid(row=1, column=0, sticky="w", pady=4)
        self.fasta_entry = ttk.Entry(form, textvariable=self.fasta_var)
        self.fasta_entry.grid(row=1, column=1, sticky="ew", pady=4)
        self.fasta_button = ttk.Button(form, text="Browse", command=self._choose_fasta)
        self.fasta_button.grid(row=1, column=2, padx=(8, 0), pady=4)

        ttk.Label(form, text="Paste FASTA").grid(row=2, column=0, sticky="nw", pady=4)
        self.pasted_fasta = scrolledtext.ScrolledText(form, height=8, wrap="word")
        self.pasted_fasta.grid(row=2, column=1, columnspan=2, sticky="ew", pady=4)
        self.pasted_fasta.insert("1.0", ">my_plasmid\n")

        self._file_row(form, 3, "Output folder", self.output_var, self._choose_output)
        self._entry_row(form, 4, "Gene start", self.start_var, "1-based nucleotide coordinate")
        self._entry_row(form, 5, "Gene end", self.end_var, "1-based nucleotide coordinate")
        self._entry_row(form, 6, "Mutation regions", self.regions_var, "Plasmid nt ranges, example: 2050-2220,2500-2700")
        self._entry_row(form, 7, "Homology length", self.homology_var, "Usually 18 bp for In-Fusion")
        self._entry_row(form, 8, "Oligo length", self.oligo_len_var, "Maximum total insert length")

        ttk.Label(form, text="Scan mode").grid(row=9, column=0, sticky="w", pady=4)
        ttk.Combobox(
            form,
            textvariable=self.scan_mode_var,
            values=("alanine", "conservative", "saturation"),
            state="readonly",
            width=16,
        ).grid(row=9, column=1, sticky="w", pady=4)

        ttk.Label(form, text="Codon usage").grid(row=10, column=0, sticky="w", pady=4)
        ttk.Combobox(form, textvariable=self.usage_var, values=("human", "mouse", "ecoli"), state="readonly", width=12).grid(
            row=10, column=1, sticky="w", pady=4
        )

        mode_notes = ttk.Label(form, text="alanine: X->A | conservative: similar amino acid | saturation: all non-WT amino acids")
        mode_notes.grid(row=11, column=1, columnspan=2, sticky="w", pady=4)

        actions = ttk.Frame(form)
        actions.grid(row=12, column=0, columnspan=3, sticky="ew", pady=(10, 0))
        actions.columnconfigure(0, weight=1)
        self.run_button = ttk.Button(actions, text="Run In-Fusion mutagenesis", command=self._start_run)
        self.run_button.grid(row=0, column=1, sticky="e")
        self._toggle_input_mode()

        output_frame = ttk.Frame(self, padding=(14, 0, 14, 14))
        output_frame.grid(row=1, column=0, sticky="nsew")
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)

        self.log_text = tk.Text(output_frame, wrap="word", height=18)
        self.log_text.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(output_frame, orient="vertical", command=self.log_text.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=scroll.set)

    def _file_row(self, parent, row, label, variable, command):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4)
        ttk.Button(parent, text="Browse", command=command).grid(row=row, column=2, padx=(8, 0), pady=4)

    def _entry_row(self, parent, row, label, variable, hint):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4)
        ttk.Label(parent, text=hint).grid(row=row, column=2, sticky="w", padx=(8, 0), pady=4)

    def _choose_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose full plasmid FASTA",
            filetypes=(("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")),
        )
        if path:
            self.fasta_var.set(path)
            if self.output_var.get() == SCRIPT_DIR:
                self.output_var.set(os.path.dirname(path))

    def _choose_output(self):
        path = filedialog.askdirectory(title="Choose output folder")
        if path:
            self.output_var.set(path)

    def _toggle_input_mode(self):
        use_file = self.input_mode_var.get() == "file"
        self.fasta_entry.configure(state="normal" if use_file else "disabled")
        self.fasta_button.configure(state="normal" if use_file else "disabled")
        self.pasted_fasta.configure(state="disabled" if use_file else "normal")

    def _start_run(self):
        if self.worker and self.worker.is_alive():
            return
        try:
            config = self._read_config()
        except Exception as error:
            messagebox.showerror("Check inputs", str(error))
            return
        self.log_text.delete("1.0", "end")
        self.run_button.configure(state="disabled")
        self._log("Starting In-Fusion mutagenesis...\n")
        self.worker = threading.Thread(target=self._run_infusion, args=(config,), daemon=True)
        self.worker.start()

    def _read_config(self):
        output_parent = self.output_var.get().strip()
        if not output_parent:
            raise ValueError("Choose an output folder.")
        os.makedirs(output_parent, exist_ok=True)
        output = create_run_folder(output_parent, "SPINE_InFusion_Mutagenesis")

        if self.input_mode_var.get() == "paste":
            fasta_text = self.pasted_fasta.get("1.0", "end").strip()
            if not fasta_text:
                raise ValueError("Paste a full plasmid sequence in FASTA format.")
            if not fasta_text.startswith(">"):
                raise ValueError("Pasted sequence must be FASTA format and start with a > header.")
            fasta = os.path.join(output, "SPINE_infusion_pasted_input.fasta")
            with open(fasta, "w") as handle:
                handle.write(fasta_text + "\n")
        else:
            fasta = self.fasta_var.get().strip()
            if not fasta:
                raise ValueError("Choose a full plasmid FASTA file.")
            if not os.path.isfile(fasta):
                raise ValueError("The FASTA file was not found.")

        return {
            "fasta": fasta,
            "output": output,
            "gene_start": int(self.start_var.get().strip()),
            "gene_end": int(self.end_var.get().strip()),
            "mutation_regions": parse_ranges(self.regions_var.get()),
            "homology_len": int(self.homology_var.get().strip()),
            "oligo_len": int(self.oligo_len_var.get().strip()),
            "usage": self.usage_var.get(),
            "scan_mode": self.scan_mode_var.get(),
        }

    def _run_infusion(self, config):
        writer = QueueWriter(self.output_queue)
        try:
            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                generate_infusion_alanine_scan(**config)
            self.output_queue.put(("done", None))
        except Exception:
            self.output_queue.put(("error", traceback.format_exc()))

    def _pump_output(self):
        try:
            while True:
                kind, payload = self.output_queue.get_nowait()
                if kind == "log":
                    self._log(payload)
                elif kind == "done":
                    self.run_button.configure(state="normal")
                    messagebox.showinfo("Finished", "In-Fusion mutagenesis files were generated.")
                elif kind == "error":
                    self._log("\n" + payload)
                    self.run_button.configure(state="normal")
                    messagebox.showerror("Run failed", "The run failed. See the log for details.")
        except queue.Empty:
            pass
        self.after(100, self._pump_output)

    def _log(self, text):
        self.log_text.insert("end", text)
        self.log_text.see("end")


if __name__ == "__main__":
    app = SpineInfusionGui()
    app.mainloop()
