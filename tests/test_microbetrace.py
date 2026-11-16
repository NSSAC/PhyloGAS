"""Unit and integration tests for phylogas.microbetrace (`export-microbetrace`)."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from phylogas.microbetrace import (  # noqa: E402
    _snp_and_mutations,
    compute_snp_threshold_concordance,
    export_microbetrace_bundle,
    main as microbetrace_main,
)


def _mutate(base_seq: str, mutations: list[tuple[int, str]]) -> str:
    """Apply 1-indexed (pos, new_base) substitutions to base_seq."""
    chars = list(base_seq)
    for pos_1, nt in mutations:
        chars[pos_1 - 1] = nt
    return "".join(chars)


class TestMicrobeTraceExport(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tmpdir.name)
        self.base_seq = ("ACGT" * 250)  # 1,000 bp synthetic genome

        # Construct a ground-truth transmission dataset with two multi-step chains
        # and three background singletons:
        #
        # Chain 10 (AY.44, Fairfax VA):
        #   100.10 (sampled) -> 101.15 (sampled, +1 SNP)
        #                    -> 102.16 (UNSAMPLED, +1 SNP) -> 103.20 (sampled, +1 SNP, so 2 hops & 2 SNPs from 100.10)
        #   101.15 (sampled) -> 104.22 (sampled, +0 SNPs from 101.15)
        #
        # Chain 20 (AY.103, Accomack VA):
        #   200.12 (sampled, +8 SNPs from base) -> 201.17 (sampled, +1 SNP) -> 202.21 (sampled, +1 SNP)
        #
        # Singletons (30, 31, 32):
        #   300.30, 301.31, 302.32 (sampled, distinct chains)
        s_100 = self.base_seq
        s_101 = _mutate(s_100, [(10, "T")])
        s_102 = _mutate(s_100, [(20, "G")])
        s_103 = _mutate(s_102, [(30, "A")])
        s_104 = s_101

        s_200 = _mutate(self.base_seq, [(100 + k * 5, "T") for k in range(8)])
        s_201 = _mutate(s_200, [(200, "A")])
        s_202 = _mutate(s_201, [(210, "C")])

        s_300 = _mutate(self.base_seq, [(400 + k * 5, "G") for k in range(6)])
        s_301 = _mutate(self.base_seq, [(500 + k * 5, "A") for k in range(6)])
        s_302 = _mutate(self.base_seq, [(600 + k * 5, "C") for k in range(6)])

        self.all_records = [
            # (alias_pid, alias_contact, comp_id, county, fips, rucc, exp_date, col_date, variant, sampled, seq)
            ("100.10", "-1",     "10", "Fairfax VA",  "51059", "1", "2021-08-01", "2021-08-04", "AY.44",  True,  s_100),
            ("101.15", "100.10", "10", "Fairfax VA",  "51059", "1", "2021-08-06", "2021-08-09", "AY.44",  True,  s_101),
            ("102.16", "100.10", "10", "Fairfax VA",  "51059", "1", "2021-08-07", "2021-08-10", "AY.44",  False, s_102),
            ("103.20", "102.16", "10", "Fairfax VA",  "51059", "1", "2021-08-11", "2021-08-14", "AY.44",  True,  s_103),
            ("104.22", "101.15", "10", "Fairfax VA",  "51059", "1", "2021-08-13", "2021-08-16", "AY.44",  True,  s_104),
            ("200.12", "-1",     "20", "Accomack VA", "51001", "9", "2021-08-03", "2021-08-06", "AY.103", True,  s_200),
            ("201.17", "200.12", "20", "Accomack VA", "51001", "9", "2021-08-08", "2021-08-11", "AY.103", True,  s_201),
            ("202.21", "201.17", "20", "Accomack VA", "51001", "9", "2021-08-12", "2021-08-15", "AY.103", True,  s_202),
            ("300.30", "-1",     "30", "Roanoke VA",  "51161", "2", "2021-08-21", "2021-08-24", "AY.25",  True,  s_300),
            ("301.31", "-1",     "31", "Roanoke VA",  "51161", "2", "2021-08-22", "2021-08-25", "AY.25",  True,  s_301),
            ("302.32", "-1",     "32", "Roanoke VA",  "51161", "2", "2021-08-23", "2021-08-26", "AY.25",  True,  s_302),
        ]

        self.infections_csv = self.root / "allevents_variants.csv"
        self.samples_csv = self.root / "1S__surs_samples.csv"
        self.fasta_file = self.root / "sequences.fasta"

        fieldnames = [
            "strain", "alias_pid", "alias_contact", "sim_pid", "component_id",
            "county", "county_fips", "rucc_code", "exposure_date", "date",
            "exposure_tick", "division", "latitude", "longitude", "age_group",
            "sex", "race", "smh_race", "ses_category", "variant_benchmark",
        ]

        with open(self.infections_csv, "w", newline="") as f_inf, \
             open(self.samples_csv, "w", newline="") as f_sam, \
             open(self.fasta_file, "w") as f_fas:
            w_inf = csv.DictWriter(f_inf, fieldnames=fieldnames)
            w_sam = csv.DictWriter(f_sam, fieldnames=fieldnames)
            w_inf.writeheader()
            w_sam.writeheader()

            # Include reference sequence to verify context sequences are ignored
            f_fas.write(f">Wuhan/Hu-1/2019\n{self.base_seq}\n")

            for al, con, cid, cty, fips, rucc, edate, cdate, var, is_sampled, seq in self.all_records:
                pid, tick = al.split(".")
                strain = f"USA/VA-EHip-{al}/2021"
                row = {
                    "strain": strain,
                    "alias_pid": al,
                    "alias_contact": con,
                    "sim_pid": pid,
                    "component_id": cid,
                    "county": cty,
                    "county_fips": fips,
                    "rucc_code": rucc,
                    "exposure_date": edate,
                    "date": cdate,
                    "exposure_tick": tick,
                    "division": "Virginia",
                    "latitude": "38.8462",
                    "longitude": "-77.3064",
                    "age_group": "Senior (65+)" if cid == "20" else "Adult (18-49)",
                    "sex": "female",
                    "race": "White alone",
                    "smh_race": "White",
                    "ses_category": "medium",
                    "variant_benchmark": var,
                }
                w_inf.writerow(row)
                if is_sampled:
                    w_sam.writerow(row)
                f_fas.write(f">{strain} extra_header_comment\n{seq}\n")

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_snp_and_mutations(self):
        s1 = "ACGTACGTAC"
        s2 = "ACTTACGAAC"  # G3T (transversion) and T8A (transversion)
        snps, tn93, muts = _snp_and_mutations(s1, s2)
        self.assertEqual(snps, 2)
        self.assertGreater(tn93, 0.2)
        self.assertEqual(muts, "G3T,T8A")

    def test_cluster_enriched_export_and_ancestral_bridge(self):
        outdir = self.root / "out_enriched"
        summary = export_microbetrace_bundle(
            outdir=outdir,
            prefix="va_test",
            fasta_path=self.fasta_file,
            samples_path=self.samples_csv,
            infections_path=self.infections_csv,
            selection_mode="cluster-enriched",
            max_sequences=8,
            max_hops=5,
            background_fraction=0.125,  # 7 cluster nodes + 1 background singleton
            seed=43,
        )

        # 7 clustered sampled nodes (4 in chain 10, 3 in chain 20) + 1 singleton = 8 nodes
        self.assertEqual(summary["n_nodes"], 8)
        self.assertEqual(summary["n_sequences"], 8)
        self.assertEqual(summary["n_multicase_true_components"], 2)
        self.assertEqual(summary["n_truth_links_total"], 5)
        self.assertEqual(summary["n_truth_links_direct_1hop"], 4)
        self.assertEqual(summary["n_truth_links_ancestral_multihop"], 1)

        # Verify that 103.20 bridged across unsampled 102.16 to 100.10 with hops=2 and snp_distance=2
        links_csv = outdir / "va_test_truth_links.csv"
        with open(links_csv, newline="") as fh:
            links = list(csv.DictReader(fh))
        ancestral = [e for e in links if e["edge_type"] == "ancestral"]
        self.assertEqual(len(ancestral), 1)
        bridge = ancestral[0]
        self.assertEqual(bridge["source"], "USA/VA-EHip-100.10/2021")
        self.assertEqual(bridge["target"], "USA/VA-EHip-103.20/2021")
        self.assertEqual(int(bridge["transmission_hops"]), 2)
        self.assertEqual(int(bridge["unsampled_intermediates"]), 1)
        self.assertEqual(int(bridge["snp_distance"]), 2)
        self.assertEqual(int(bridge["delta_exposure_days"]), 10)

        # Verify FASTA headers match node `id` column 1:1 with no trailing spaces
        fasta_out = outdir / "va_test_sequences.fasta"
        fasta_ids = [
            line[1:].strip()
            for line in fasta_out.read_text().splitlines()
            if line.startswith(">")
        ]
        nodes_csv = outdir / "va_test_nodes.csv"
        with open(nodes_csv, newline="") as fh:
            nodes = list(csv.DictReader(fh))
        node_ids = [n["id"] for n in nodes]
        self.assertEqual(fasta_ids, node_ids)
        self.assertNotIn("Wuhan/Hu-1/2019", fasta_ids)

        # Verify SNP threshold clustering concordance at tau=2 recovers perfect clusters (F1=1.0, ARI=1.0)
        conc_by_tau = {c["snp_threshold"]: c for c in summary["snp_threshold_concordance"]}
        self.assertAlmostEqual(conc_by_tau[2]["pairwise_f1"], 1.0, places=4)
        self.assertAlmostEqual(conc_by_tau[2]["adjusted_rand_index"], 1.0, places=4)

    def test_complete_chains_includes_unsampled_intermediates(self):
        outdir = self.root / "out_complete"
        summary = export_microbetrace_bundle(
            outdir=outdir,
            prefix="va_chains",
            fasta_path=self.fasta_file,
            samples_path=self.samples_csv,
            infections_path=self.infections_csv,
            selection_mode="complete-chains",
            min_chain_size=3,
            max_chain_size=10,
            max_sequences=8,
        )
        # Chain 10 (5 total nodes including unsampled 102.16) + Chain 20 (3 nodes) = 8 nodes
        self.assertEqual(summary["n_nodes"], 8)
        self.assertEqual(summary["n_truth_links_direct_1hop"], 6)
        self.assertEqual(summary["n_truth_links_ancestral_multihop"], 0)

        nodes_csv = outdir / "va_chains_nodes.csv"
        with open(nodes_csv, newline="") as fh:
            nodes = {n["alias_pid"]: n for n in csv.DictReader(fh)}
        self.assertIn("102.16", nodes)
        self.assertEqual(nodes["102.16"]["sampled_in_surveillance"], "False")
        self.assertEqual(nodes["103.20"]["sampled_in_surveillance"], "True")

    def test_cli_entrypoint_and_filters(self):
        outdir = self.root / "out_cli"
        rc = microbetrace_main([
            "--fasta", str(self.fasta_file),
            "--samples", str(self.samples_csv),
            "--infections", str(self.infections_csv),
            "--outdir", str(outdir),
            "--prefix", "accomack_only",
            "--county", "Accomack",
            "--id-column", "alias_pid",
        ])
        self.assertEqual(rc, 0)
        summary = json.loads((outdir / "accomack_only_summary.json").read_text())
        self.assertEqual(summary["n_nodes"], 3)
        self.assertEqual(summary["urbanicity_breakdown"], {"Non-metro / Rural (RUCC 4-9)": 3})

    def test_household_coincidental_exposures_and_secondary_cases(self):
        # Enrich infections_csv with explicit household IDs (`hid`):
        # - 100.10 and 101.15 share household HH_A (direct household transmission)
        # - 300.30 (chain 30) and 301.31 (chain 31) share household HH_COINC
        #   (coincidental household co-exposure from two independent introductions!)
        with open(self.infections_csv, newline="") as fh:
            rows = list(csv.DictReader(fh))
        fieldnames = list(rows[0].keys()) + ["hid", "rlid", "admin3", "occupation_socp"]
        for r in rows:
            al = r["alias_pid"]
            if al in ("100.10", "101.15"):
                r["hid"] = "HH_A"
            elif al in ("300.30", "301.31"):
                r["hid"] = "HH_COINC"
            else:
                r["hid"] = f"HH_{al}"
            r["rlid"] = f"RL_{r['hid']}"
            r["admin3"] = "90100"
            r["occupation_socp"] = "291141"
        with open(self.infections_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

        outdir = self.root / "out_hh"
        summary = export_microbetrace_bundle(
            outdir=outdir,
            prefix="va_hh",
            fasta_path=self.fasta_file,
            samples_path=self.samples_csv,
            infections_path=self.infections_csv,
            selection_mode="all",
            max_sequences=20,
        )
        # 1 coincidental household link between 300.30 (chain 30) and 301.31 (chain 31)
        self.assertEqual(summary["n_coincidental_exposure_links"], 1)
        self.assertEqual(summary["transmission_exposure_settings"].get("Household"), 1)

        with open(outdir / "va_hh_truth_links.csv", newline="") as fh:
            links = list(csv.DictReader(fh))
        coinc = [e for e in links if e["edge_type"] == "coincidental_exposure"]
        self.assertEqual(len(coinc), 1)
        self.assertEqual(coinc[0]["exposure_setting"], "Household")
        self.assertEqual(coinc[0]["same_true_chain"], "False")
        # 300.30 and 301.31 differ by 9 SNPs -> MicrobeTrace 0-2 SNP filter rules it out
        self.assertEqual(int(coinc[0]["snp_distance"]), 9)

        with open(outdir / "va_hh_nodes.csv", newline="") as fh:
            nodes = {n["alias_pid"]: n for n in csv.DictReader(fh)}
        # 100.10 infected both 101.15 and 102.16 in the full graph -> true_secondary_cases == 2
        self.assertEqual(int(nodes["100.10"]["true_secondary_cases"]), 2)
        self.assertEqual(nodes["100.10"]["household_id"], "HH_A")


if __name__ == "__main__":
    unittest.main()

