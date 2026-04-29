# SPINE In-Fusion Mutagenesis GUI

This folder contains a separate In-Fusion cloning workflow for SPINE-style mutagenesis. It does not modify the original `SPINE/SPINE.py` code and does not use Golden Gate cloning, BsaI, BsmBI, barcode primers, or Type IIS overhangs.

The GUI supports three mutation modes:

- **Alanine scan**: every non-alanine residue in the selected regions is mutated to alanine.
- **Conservative scan**: residues are mutated to similar amino acids, such as `E->D`, `D->E`, `R->K`, `K->R`, `S->T`, and related substitutions.
- **Full saturation scan**: every selected residue is mutated to the other 19 standard amino acids.

## Start the GUI

Double-click:

```text
run_SPINE_mutagenesis_infusion_gui.bat
```

Or run from PowerShell:

```powershell
python SPINE_mutagenesis_infusion_gui.py
```

## Inputs

The GUI accepts either:

- a full plasmid FASTA file, or
- pasted FASTA text.

Required coordinates:

- **Gene start**: 1-based nucleotide coordinate in the full plasmid.
- **Gene end**: 1-based nucleotide coordinate in the full plasmid.
- **Mutation regions**: comma-separated 1-based nucleotide ranges in the full plasmid, for example:

```text
2050-2220,2500-2700
```

The mutation regions must fall inside the selected gene. Regions are automatically expanded to full codons when needed.

## Oligo Chunking

Long mutation regions are split into codon-aligned chunks so insert oligos fit the selected oligo length.

Default settings:

- homology length: `18 bp` on each side
- oligo length: `230 bp`
- mutable chunk size: up to about `192 bp`

Each chunk receives its own vector primer pair. The insert format is:

```text
left In-Fusion homology + mutated chunk + right In-Fusion homology
```

Only one codon is changed per insert. All other bases in that insert remain wild type.

## Codon Usage

The mutant codon is selected using the chosen codon usage table:

- `human`
- `mouse`
- `ecoli`

Mouse codon usage is available from the GUI dropdown and the command line option `--usage mouse`.

## Primer Design

Vector primers are designed around each chunk.

The complementary primer region is selected to prefer:

- Tm: `56-57 C`
- GC content: `40-60%`

If no candidate satisfies both windows exactly, the closest compromise is selected. Primer FASTA descriptions include both Tm and GC content.

## Outputs

Each GUI run creates a new timestamped folder:

```text
SPINE_InFusion_Mutagenesis_YYYYMMDD_HHMMSS
```

The output files are:

- `InFusion_Mutagenesis_Inserts.fasta`
- `InFusion_Vector_Primers.fasta`
- `InFusion_Mutagenesis_Summary.csv`

If FASTA was pasted into the GUI, the pasted input is also saved in the run folder.

## Command Line

The same pipeline can be run without the GUI:

```powershell
python SPINE_mutagenesis_infusion.py `
  --fasta plasmid.fasta `
  --gene-start 1800 `
  --gene-end 3065 `
  --mutation-regions 2050-2220,2500-2700 `
  --output output_folder `
  --scan-mode saturation `
  --usage mouse `
  --homology-len 18 `
  --oligo-len 230
```

Allowed scan modes:

```text
alanine
conservative
saturation
```

Allowed codon usage values:

```text
human
mouse
ecoli
```
