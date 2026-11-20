"""Unit tests for canonical public-health specimen-pool & subsampling algorithms."""

from __future__ import annotations

import csv
import lzma
import random
import unittest
from pathlib import Path
import tempfile

from phylogas.sampling_canonical import (
    CANONICAL_STRATEGIES,
    PoolConfig,
    allocate_waterfill_quotas,
    haversine_km_scalar,
    run_canonical_from_csv,
    select_regional_hubs,
)


class TestCanonicalSampling(unittest.TestCase):
    def test_waterfill_geocap_redistributes_capped_counties(self) -> None:
        caps = {"King": 200, "Spokane": 100, "Clark": 4, "SanJuan": 2}
        rng = random.Random(42)
        q = allocate_waterfill_quotas(caps, weights=None, budget=100, rng=rng)
        self.assertEqual(sum(q.values()), 100)
        # Small counties below 100/4 = 25 must get 100% of their available capacity
        self.assertEqual(q["Clark"], 4)
        self.assertEqual(q["SanJuan"], 2)
        # Remaining 94 split evenly (47 and 47) across King and Spokane
        self.assertEqual(q["King"], 47)
        self.assertEqual(q["Spokane"], 47)

    def test_waterfill_caseprop_respects_caps_and_budget(self) -> None:
        # Clark has 30% of statewide cases (weight 300) but only 5 shipped swabs in sentinel pool
        caps = {"King": 300, "Spokane": 150, "Clark": 5}
        weights = {"King": 500.0, "Spokane": 200.0, "Clark": 300.0}
        rng = random.Random(7)
        q = allocate_waterfill_quotas(caps, weights=weights, budget=100, rng=rng)
        self.assertEqual(sum(q.values()), 100)
        # Ideal Clark = 30 > cap 5 -> Clark must lock at all 5 available swabs
        self.assertEqual(q["Clark"], 5)
        # Remaining 95 split 5:2 between King (~68) and Spokane (~27)
        self.assertEqual(q["King"], 68)
        self.assertEqual(q["Spokane"], 27)

    def test_select_regional_hubs_disperses_geographically(self) -> None:
        cases = {
            "Seattle_A": 10000,
            "Seattle_Suburb": 9000,  # right next to Seattle_A
            "Spokane_East": 4000,    # 360 km east
            "Yakima_Central": 2500,  # 180 km southeast
        }
        coords = {
            "Seattle_A": (47.60, -122.33),
            "Seattle_Suburb": (47.62, -122.30),
            "Spokane_East": (47.65, -117.42),
            "Yakima_Central": (46.60, -120.50),
        }
        hubs = select_regional_hubs(cases, coords, k_hubs=3, pop_alpha=0.35, min_hub_share=0.01)
        self.assertEqual(hubs[0], "Seattle_A")
        self.assertIn("Spokane_East", hubs)
        self.assertIn("Yakima_Central", hubs)
        self.assertNotIn("Seattle_Suburb", hubs)
        self.assertGreater(haversine_km_scalar(*coords["Seattle_A"], *coords["Spokane_East"]), 300.0)

    def test_end_to_end_canonical_csv_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tpath = Path(td)
            ll_csv = tpath / "synthetic_ll.csv"
            counties = [
                ("King WA", 47.60, -122.33, "1", 300),
                ("Spokane WA", 47.65, -117.42, "2", 120),
                ("Yakima WA", 46.60, -120.50, "3", 80),
                ("Clark WA", 45.63, -122.67, "2", 90),
                ("Ferry WA", 48.50, -118.50, "8", 30),
            ]
            with open(ll_csv, "w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["alias_pid", "alias_contact", "date", "county", "latitude", "longitude", "rucc_code", "exit_state", "component_id"])
                pid = 0
                for wk in range(4):
                    dt = f"2021-07-{2 + wk * 7:02d}"
                    for cty, la, lo, rucc, n_c in counties:
                        for j in range(n_c):
                            pid += 1
                            ex = "hM2_a" if j % 10 == 0 else ("A2_a" if j % 7 == 0 else "I2_a")
                            w.writerow([f"p{pid}", f"p{max(1, pid-1)}", dt, cty, str(la), str(lo), rucc, ex, "10"])

            out_dir = tpath / "03_sampled_canonical"
            diag_dir = tpath / "diag"
            written, diag = run_canonical_from_csv(
                linelist_path=ll_csv,
                outdir=out_dir,
                start_date_iso="2021-07-08",
                min_pool=50,
                weekly_budget=100,
                strides=(1, 4),
                strategies=CANONICAL_STRATEGIES,
                pool_cfg=PoolConfig(sentinel_k_hubs=2, convenience_k_hubs=1),
                seed=43,
                diagnostics_dir=diag_dir,
                state_label="wa_test",
            )
            self.assertEqual(len(written), 2 * len(CANONICAL_STRATEGIES))
            # Verify sentinel under-represents inpatient (< 1.0) while convenience over-represents inpatient (> 1.0)
            self.assertLess(diag["clinical_tier_OE"]["inpatient"]["sentinel_OE"], 0.85)
            self.assertGreater(diag["clinical_tier_OE"]["inpatient"]["convenience_OE"], 1.30)

            # Verify sample CSV is valid xz with exact weekly budget (4 weeks * 100 = 400)
            s_geocap = out_dir / "1S__statewide_geocap_samples.csv.xz"
            with lzma.open(s_geocap, "rt", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 400)


if __name__ == "__main__":
    unittest.main()
