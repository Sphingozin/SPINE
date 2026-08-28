import csv
import tempfile
import unittest
from pathlib import Path

from SPINE_mutagenesis_infusion import (
    generate_infusion_alanine_scan,
    mutation_targets,
    parse_membrane_environments,
)


class MembraneConservativeTests(unittest.TestCase):
    def test_environment_parser(self):
        parsed = parse_membrane_environments(
            "tm_lipid:1-3;tm_packed:4;hydrated:5-6;functional:7,9"
        )
        self.assertEqual(parsed[1], "tm_lipid")
        self.assertEqual(parsed[6], "hydrated")
        self.assertEqual(parsed[9], "functional")

    def test_conflicting_environment_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "assigned to both"):
            parse_membrane_environments("tm_lipid:1-3;hydrated:3-5")

    def test_environment_changes_targets(self):
        self.assertEqual(mutation_targets("Ser", "membrane_conservative", "tm_lipid"), ["Thr", "Ala"])
        self.assertEqual(mutation_targets("Ser", "membrane_conservative", "tm_packed"), ["Thr"])
        self.assertEqual(mutation_targets("Arg", "membrane_conservative", "tm_lipid"), [])
        self.assertEqual(mutation_targets("Arg", "membrane_conservative", "functional"), ["Lys"])

    def test_generation_records_environment(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fasta"
            # S1, R2, L3, P4
            fasta.write_text(">test\nAGCCGCCTGCCC\n", encoding="ascii")
            output = root / "output"
            inserts, primers, rows = generate_infusion_alanine_scan(
                fasta=str(fasta), gene_start=1, gene_end=12,
                mutation_regions=[(1, 12)], output=str(output),
                homology_len=3, oligo_len=18, scan_mode="membrane_conservative",
                membrane_environments=parse_membrane_environments(
                    "tm_lipid:1;functional:2;tm_packed:3;hydrated:4"
                ),
            )
            self.assertTrue(inserts)
            self.assertTrue(primers)
            self.assertEqual({row[-2] for row in rows}, {"tm_lipid", "functional", "tm_packed", "hydrated"})
            with (output / "InFusion_Mutagenesis_Summary.csv").open(newline="") as handle:
                summary = list(csv.DictReader(handle))
            self.assertEqual(summary[0]["membrane_environment"], "tm_lipid")

    def test_missing_annotation_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fasta"
            fasta.write_text(">test\nAGCCGC\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "position 2"):
                generate_infusion_alanine_scan(
                    fasta=str(fasta), gene_start=1, gene_end=6,
                    mutation_regions=[(1, 6)], output=str(root / "output"),
                    homology_len=3, oligo_len=12, scan_mode="membrane_conservative",
                    membrane_environments={1: "tm_lipid"},
                )


if __name__ == "__main__":
    unittest.main()

