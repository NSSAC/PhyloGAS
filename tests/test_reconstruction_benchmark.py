"""Unit and integration tests for phylogenetic & transmission reconstruction benchmarks."""

from __future__ import annotations

import csv
import datetime as _dt
import json
import tempfile
import unittest
from pathlib import Path

from phylogas.benchmark import outbreak_prep as op
from phylogas.benchmark import reconstruction as recon


class TestReconstructionBenchmark(unittest.TestCase):
    def test_id_normalization(self) -> None:
        lkp = {"S000001": "USA/VA-EHip-101.150/2021", "S000002": "202.160"}
        self.assertEqual(recon.alias_from_any_id("hCoV-19/USA/VA-EHip-12.34/2021"), "12.34")
        self.assertEqual(recon.alias_from_any_id("USA/MA-EHip-99.200/2021"), "99.200")
        self.assertEqual(recon.alias_from_any_id("S000001", lkp), "101.150")
        self.assertEqual(recon.alias_from_any_id("S000002", lkp), "202.160")
        self.assertIsNone(recon.alias_from_any_id("Wuhan/Hu-1/2019", lkp))
        self.assertIsNone(recon.alias_from_any_id("ref_genome", lkp))
        self.assertIsNone(recon.alias_from_any_id("HOST_000001", lkp))

    def _make_synthetic_truth_csv(self, tmp: Path) -> Path:
        """Build a 2-chain synthetic ABM truth table:
          Chain 1:
            1.10 (2021-06-01, unsampled root)
              |-- 2.15 (2021-06-06, sampled) -> infects 3.20 (2021-06-11, sampled)
              `-- 4.16 (2021-06-07, unsampled) -> infects 5.22 (2021-06-13, sampled) & 6.23 (2021-06-14, sampled)
          Chain 2 (independent introduction):
            7.12 (2021-06-03, sampled) -> infects 8.18 (2021-06-09, sampled)
        Sampled set S = {2.15, 3.20, 5.22, 6.23, 7.12, 8.18} (6 tips).
        """
        p = tmp / "allevents_variants.csv"
        rows = [
            ("1.10", "-1", "2021-06-01", "1", "AY.44", "51059"),
            ("2.15", "1.10", "2021-06-06", "1", "AY.44", "51059"),
            ("3.20", "2.15", "2021-06-11", "1", "AY.44", "51059"),
            ("4.16", "1.10", "2021-06-07", "1", "AY.44", "51059"),
            ("5.22", "4.16", "2021-06-13", "1", "AY.44", "51059"),
            ("6.23", "4.16", "2021-06-14", "1", "AY.44", "51059"),
            ("7.12", "-1", "2021-06-03", "2", "AY.103", "51059"),
            ("8.18", "7.12", "2021-06-09", "2", "AY.103", "51059"),
        ]
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["alias_pid", "alias_contact", "exposure_date", "component_id", "variant_benchmark", "county_fips"])
            w.writerows(rows)
        return p

    def test_exact_phylogeny_and_polytomy_collapse(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            inf_csv = self._make_synthetic_truth_csv(tmp)
            tips = {"2.15", "3.20", "5.22", "6.23", "7.12", "8.18"}
            truth = recon.load_truth_forest(inf_csv, target_aliases=tips)
            true_tree = recon.build_pruned_truth_genealogy(truth, tips)

            # Scoring true_tree against itself must achieve perfect scores
            perfect = recon.score_phylogeny(true_tree, truth, min_mutations_collapse=0.5)
            self.assertEqual(perfect["n_tips_scored"], 6)
            self.assertEqual(perfect["strict_nrf"], 0.0)
            self.assertEqual(perfect["collapsed_nrf"], 0.0)
            self.assertEqual(perfect["collapsed_clade_f1"], 1.0)
            self.assertEqual(perfect["info_rf_similarity"], 1.0)
            self.assertEqual(perfect["quartet_concordance"], 1.0)
            self.assertEqual(perfect["kc_spearman_r"], 1.0)
            self.assertEqual(perfect["cross_chain_merge_pct"], 0.0)
            self.assertEqual(perfect["clade_tmrca_mae_days"], 0.0)

            # Now test a forced-binary IQ-TREE style Newick that introduces a
            # zero-mutation (0.0 SNP) binary split at the root merging Chain 1 & Chain 2
            # only via a 0-SNP branch while keeping >=1 SNP on real clades:
            nwk = (
                "(((('USA/VA-EHip-2.15/2021':0.0,'USA/VA-EHip-3.20/2021':0.0001)N_215:0.0001,"
                "('USA/VA-EHip-5.22/2021':0.0001,'USA/VA-EHip-6.23/2021':0.0001)N_416:0.0001)N_110:0.0002,"
                "('USA/VA-EHip-7.12/2021':0.0,'USA/VA-EHip-8.18/2021':0.0001)N_712:0.0002)N_SUB:0.0,"
                "'Wuhan/Hu-1/2019':0.001)ROOT;"
            )
            node_data = {
                "clock": {"rate": 0.0007},
                "nodes": {
                    "N_110": {
                        "mutation_length": 5.0,
                        "num_date": recon.date_to_decimal_year(_dt.date(2021, 6, 1)),
                        "num_date_confidence": [
                            recon.date_to_decimal_year(_dt.date(2021, 5, 29)),
                            recon.date_to_decimal_year(_dt.date(2021, 6, 4)),
                        ],
                    },
                    "N_215": {
                        "mutation_length": 1.0,
                        "num_date": recon.date_to_decimal_year(_dt.date(2021, 6, 6)),
                        "num_date_confidence": [
                            recon.date_to_decimal_year(_dt.date(2021, 6, 5)),
                            recon.date_to_decimal_year(_dt.date(2021, 6, 8)),
                        ],
                    },
                    "N_416": {
                        "mutation_length": 2.0,
                        "num_date": recon.date_to_decimal_year(_dt.date(2021, 6, 7)),
                        "num_date_confidence": [
                            recon.date_to_decimal_year(_dt.date(2021, 6, 5)),
                            recon.date_to_decimal_year(_dt.date(2021, 6, 9)),
                        ],
                    },
                    "N_712": {
                        "mutation_length": 6.0,
                        "num_date": recon.date_to_decimal_year(_dt.date(2021, 6, 3)),
                        "num_date_confidence": [
                            recon.date_to_decimal_year(_dt.date(2021, 6, 1)),
                            recon.date_to_decimal_year(_dt.date(2021, 6, 5)),
                        ],
                    },
                    "N_SUB": {"mutation_length": 0.0, "num_date": 2021.30},
                },
            }
            inf_tree = recon.parse_newick(nwk, node_data=node_data)
            sc = recon.score_phylogeny(inf_tree, truth, min_mutations_collapse=0.5)
            self.assertEqual(sc["collapsed_clade_f1"], 1.0)
            self.assertEqual(sc["collapsed_nrf"], 0.0)
            self.assertEqual(sc["clade_tmrca_mae_days"], 0.0)
            self.assertEqual(sc["clade_tmrca_95ci_cov"], 1.0)
            self.assertEqual(sc["tree_chain_ari"], 1.0)

            # Also verify TransPhylo decimal-year Newick export drops Wuhan outgroup
            nwk_p = tmp / "tree.nwk"
            nd_p = tmp / "branch_lengths.json"
            out_nwk = tmp / "dated_years.nwk"
            out_meta = tmp / "dated_meta.json"
            # Populate tip dates in node_data
            for st_id, d_iso in (
                ("USA/VA-EHip-2.15/2021", "2021-06-10"),
                ("USA/VA-EHip-3.20/2021", "2021-06-15"),
                ("USA/VA-EHip-5.22/2021", "2021-06-17"),
                ("USA/VA-EHip-6.23/2021", "2021-06-18"),
                ("USA/VA-EHip-7.12/2021", "2021-06-07"),
                ("USA/VA-EHip-8.18/2021", "2021-06-13"),
            ):
                node_data["nodes"][st_id] = {
                    "num_date": recon.date_to_decimal_year(_dt.date.fromisoformat(d_iso)),
                    "mutation_length": 1.0,
                }
            nwk_p.write_text(nwk, encoding="utf-8")
            nd_p.write_text(json.dumps(node_data), encoding="utf-8")
            meta_tp = recon.export_dated_newick_for_transphylo(nwk_p, nd_p, out_nwk, out_meta)
            self.assertEqual(meta_tp["n_tips"], 6)
            exported_txt = out_nwk.read_text(encoding="utf-8")
            self.assertNotIn("Wuhan/Hu-1/2019", exported_txt)

    def test_direct_and_indirect_transmission_scoring(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            inf_csv = self._make_synthetic_truth_csv(tmp)
            tips = {"2.15", "3.20", "5.22", "6.23", "7.12", "8.18"}
            truth = recon.load_truth_forest(inf_csv, target_aliases=tips)
            t_dir, t_ind, hops, _ = recon.build_truth_donors(truth, tips)

            # Direct donors in S: 3.20 <- 2.15, 8.18 <- 7.12; all others have unsampled/root parent!
            self.assertEqual(t_dir["3.20"], "2.15")
            self.assertEqual(t_dir["8.18"], "7.12")
            self.assertEqual(t_dir["5.22"], recon.UNOBSERVED)
            self.assertEqual(t_ind["5.22"], recon.UNOBSERVED)  # 4.16 and 1.10 are both unsampled!
            self.assertEqual(hops["3.20"], 1)

            # Add 1.10 into tips: now 5.22 has UNOBSERVED direct donor (4.16 unsampled),
            # but 1.10 as 2-hop INDIRECT sampled donor!
            tips_with_root = tips | {"1.10"}
            truth2 = recon.load_truth_forest(inf_csv, target_aliases=tips_with_root)
            t_dir2, t_ind2, hops2, _ = recon.build_truth_donors(truth2, tips_with_root)
            self.assertEqual(t_dir2["5.22"], recon.UNOBSERVED)
            self.assertEqual(t_ind2["5.22"], "1.10")
            self.assertEqual(hops2["5.22"], 2)

            # Perfect posterior distribution on tips_with_root
            perfect_probs = {
                "direct": {r: {d: 1.0} for r, d in t_dir2.items()},
                "indirect": {r: {d: 1.0} for r, d in t_ind2.items()},
            }
            sc = recon.score_transmission(perfect_probs, truth2, tips_with_root)
            self.assertEqual(sc["direct_donor_prob_mean"], 1.0)
            self.assertEqual(sc["direct_top1_donor_acc"], 1.0)
            self.assertEqual(sc["direct_obs_donor_recall"], 1.0)
            self.assertEqual(sc["direct_obs_link_precision"], 1.0)
            self.assertEqual(sc["direct_donor_95set_cov"], 1.0)
            self.assertEqual(sc["indirect_top1_donor_acc"], 1.0)
            self.assertEqual(sc["indirect_obs_donor_recall"], 1.0)

    def test_outbreak_prep_and_scoring_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            # Construct a synthetic county outbreak of 180 cases:
            #   - 150 cases in focal chain C1 (0-2 SNPs from founder)
            #   - 15 cases in a 1-SNP cousin re-introduction (same statewide chain C1, outside parent)
            #   - 15 cases in a distant 12-SNP unrelated lineage C99 (should be dropped by <=3 SNP pre-isolation!)
            ref_seq = "A" * 200
            founder_c1 = "C" + "A" * 199
            distant_c99 = ("G" * 12) + ("A" * 188)

            ref_fa = tmp / "ref.fasta"
            ref_fa.write_text(f">Wuhan/Hu-1/2019\n{ref_seq}\n", encoding="utf-8")

            inf_csv = tmp / "allevents_variants.csv"
            painted_fa = tmp / "painted.fasta"
            base_date = _dt.date(2021, 7, 1)

            inf_rows = []
            fa_lines = []

            # Focal local chain (150 cases: 100.1 -> 101.2 -> ... branching)
            for idx in range(1, 151):
                al = f"{100 + idx}.{10 + (idx // 3)}"
                parent = "-1" if idx == 1 else f"{100 + max(1, idx // 2)}.{10 + (max(1, idx // 2) // 3)}"
                d_str = (base_date + _dt.timedelta(days=idx // 4)).isoformat()
                st = f"USA/VA-EHip-{al}/2021"
                inf_rows.append((al, parent, d_str, "1", "AY.44", "51690", "36.68", "-79.87", st))
                seq_list = list(founder_c1)
                if idx > 50:
                    seq_list[5] = "T"  # 1 SNP sub-branch
                fa_lines.append(f">{st}\n{''.join(seq_list)}\n")

            # Cousin 1-SNP re-introduction into same county (15 cases, parent outside county "9999.1")
            for idx in range(1, 16):
                al = f"{500 + idx}.{20 + idx}"
                parent = "9999.1" if idx == 1 else f"{500 + idx - 1}.{20 + idx - 1}"
                d_str = (base_date + _dt.timedelta(days=10 + idx)).isoformat()
                st = f"USA/VA-EHip-{al}/2021"
                inf_rows.append((al, parent, d_str, "1", "AY.44", "51690", "36.68", "-79.87", st))
                seq_list = list(founder_c1)
                seq_list[15] = "G"  # 1 SNP cousin
                fa_lines.append(f">{st}\n{''.join(seq_list)}\n")

            # Outside parent row in another county so truth forest knows 501.21 came from 9999.1
            inf_rows.append(("9999.1", "-1", "2021-06-25", "1", "AY.44", "51089", "36.68", "-79.87", "USA/VA-EHip-9999.1/2021"))

            # Distant 12-SNP unrelated lineage in same county (15 cases)
            for idx in range(1, 16):
                al = f"{800 + idx}.{25 + idx}"
                parent = "-1" if idx == 1 else f"{800 + idx - 1}.{25 + idx - 1}"
                d_str = (base_date + _dt.timedelta(days=5 + idx)).isoformat()
                st = f"USA/VA-EHip-{al}/2021"
                inf_rows.append((al, parent, d_str, "99", "AY.3", "51690", "36.68", "-79.87", st))
                fa_lines.append(f">{st}\n{distant_c99}\n")

            with open(inf_csv, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow([
                    "alias_pid", "alias_contact", "exposure_date", "component_id",
                    "variant_benchmark", "county_fips", "latitude", "longitude", "strain"
                ])
                w.writerows(inf_rows)
            painted_fa.write_text("".join(fa_lines), encoding="utf-8")

            ob_out = tmp / "07_outbreak_bench"
            manifests = op.prepare_state_outbreaks(
                state="va",
                infections_path=inf_csv,
                painted_fasta=painted_fa,
                reference_fasta=ref_fa,
                outdir=ob_out,
                n_clusters=1,
                sampling_rates=[0.50],
                max_cluster_snps=3,
            )
            self.assertEqual(len(manifests), 1)
            mf = manifests[0]
            rdir = Path(mf["run_dir"])

            # Verify JUNIPER safety contract: ONLY `aligned.fasta` and `metadata.csv` inside `input_data/`
            in_files = sorted(p.name for p in (rdir / "input_data").iterdir())
            self.assertEqual(in_files, ["aligned.fasta", "metadata.csv"])

            # Verify distant 12-SNP C99 cases were dropped while cousin C1 cases remained compatible
            self.assertGreater(mf["n_dropped_distant_snps"], 0)
            self.assertEqual(mf["n_distinct_statewide_comps"], 1)
            self.assertGreaterEqual(mf["n_compatible_local_intros"], 2)

            # Simulate a JUNIPER output (`juniper_nodes.csv` + `juniper_donor_probs.csv` + `juniper_params.json`)
            # and verify `score_outbreak_directory` scores both MicrobeTrace_1SNP and JUNIPER!
            with open(rdir / "sample_lookup.csv", encoding="utf-8") as fh_lkp:
                lkp_rows = list(csv.DictReader(fh_lkp))
            with open(rdir / "juniper_nodes.csv", "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["idx", "infector_idx", "infection_date", "observed", "sample"])
                w.writerow([1, "", "2021-06-28", True, "ref_genome"])
                for idx_r, lr in enumerate(lkp_rows, start=2):
                    w.writerow([idx_r, max(1, idx_r - 1), lr["exposure_date"], True, lr["sample"]])

            with open(rdir / "juniper_donor_probs.csv", "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["recipient_id", "donor_id", "direct_prob", "indirect_prob"])
                for idx_r, lr in enumerate(lkp_rows):
                    don = "UNOBSERVED" if idx_r == 0 else lkp_rows[idx_r - 1]["sample"]
                    w.writerow([lr["sample"], don, 0.8, 0.9])
                    w.writerow([lr["sample"], "UNOBSERVED", 0.2, 0.1])

            (rdir / "juniper_params.json").write_text(
                json.dumps({"runtime_secs": 12.5, "mu_est_yr": 0.0007, "pi_est": 0.48, "R_est": 1.4, "ess_min": 95.0}),
                encoding="utf-8",
            )

            score_rows = op.score_outbreak_directory(rdir, inf_csv)
            methods = [r["method"] for r in score_rows]
            self.assertIn("MicrobeTrace_1SNP", methods)
            self.assertIn("JUNIPER", methods)
            jun_row = next(r for r in score_rows if r["method"] == "JUNIPER")
            self.assertIsNotNone(jun_row["direct_donor_prob_mean"])
            self.assertIsNotNone(jun_row["collapsed_clade_f1"])


if __name__ == "__main__":
    unittest.main()
