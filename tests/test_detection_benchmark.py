"""Unit tests for the T95 / UT95 / ΔT95 detection lead-time benchmark."""

from __future__ import annotations

import csv
import math
import tempfile
import unittest
from pathlib import Path

from phylogas.benchmark.detection import (
    _hypergeom_miss_prob,
    _strategy_label_from_path,
    benchmark_detection,
    compute_stratified_t95_weeks,
    compute_ut95_weeks,
)


class TestDetectionBenchmark(unittest.TestCase):
    def test_strategy_label_strips_scenario_prefix(self) -> None:
        self.assertEqual(
            _strategy_label_from_path(Path("scenario1_1S__surs_samples.csv.xz")),
            "1S__surs",
        )
        self.assertEqual(
            _strategy_label_from_path(Path("1S__sentinel_geocap_samples.csv.xz")),
            "1S__sentinel_geocap",
        )
    def test_hypergeom_miss_prob_exact_values(self) -> None:
        # Pool N=10, k=2 targets, n=3 draws: C(8,3)/C(10,3) = 56/120 = 7/15
        self.assertAlmostEqual(_hypergeom_miss_prob(10, 2, 3), 56.0 / 120.0, places=12)
        # Edge cases
        self.assertEqual(_hypergeom_miss_prob(10, 0, 3), 1.0)
        self.assertEqual(_hypergeom_miss_prob(10, 2, 0), 1.0)
        # Guaranteed hit when n > N - k
        self.assertEqual(_hypergeom_miss_prob(10, 3, 8), 0.0)
        self.assertEqual(_hypergeom_miss_prob(10, 10, 1), 0.0)

    def test_ut95_weeks_matches_closed_form(self) -> None:
        # Constant prevalence p = 0.05, weekly budget B = 20 -> per-week miss = 0.95^20 = 0.358486
        # After 1 wk (lag 0): 1 - 0.3585 = 0.6415
        # After 2 wks (lag 1): 1 - 0.3585^2 = 0.8715
        # After 3 wks (lag 2): 1 - 0.3585^3 = 0.9539 >= 0.95 -> lag = 2 weeks
        target = {w: 5 for w in range(10)}
        total = {w: 100 for w in range(10)}
        budget = {w: 20 for w in range(10)}
        self.assertEqual(compute_ut95_weeks(target, total, budget, origin_week=1, n_weeks=10), 2)

        # 4S stride budget (only week 0, 4, 8 sampled with B=60): 1 - 0.95^60 = 0.9539 at week 4 (lag 3)
        stride_budget = {0: 60, 4: 60, 8: 60}
        self.assertEqual(
            compute_ut95_weeks(target, total, stride_budget, origin_week=1, n_weeks=10), 3
        )

    def test_stratified_t95_geocap_beats_convenience_on_rural_chain(self) -> None:
        # Rural Montgomery (N=50, k=10) vs Metro Fairfax (N=1000, k=0) each week
        weekly_N = {w: {"Montgomery": 50, "Fairfax": 1000} for w in range(6)}
        target_k = {w: {"Montgomery": 10} for w in range(6)}

        # GeoCap: n=10 in Montgomery, n=10 in Fairfax -> miss prob per week = C(40,10)/C(50,10) ≈ 0.0825
        # 2 weeks -> 1 - 0.0825^2 = 0.993 >= 0.95 -> T95 = 1 week
        geocap_n = {w: {"Montgomery": 10, "Fairfax": 10} for w in range(6)}
        t95_geo = compute_stratified_t95_weeks(target_k, weekly_N, geocap_n, origin_week=0, n_weeks=6)
        self.assertEqual(t95_geo, 1)

        # Convenience: n=0 in Montgomery, n=20 in Fairfax -> never detects (None)
        conv_n = {w: {"Montgomery": 0, "Fairfax": 20} for w in range(6)}
        t95_conv = compute_stratified_t95_weeks(target_k, weekly_N, conv_n, origin_week=0, n_weeks=6)
        self.assertIsNone(t95_conv)

        # Compare against UT95 (statewide p = 10/1050, B=20 -> miss ≈ 0.826/wk -> needs 16 wks > 6)
        ut95 = compute_ut95_weeks(
            {w: 10 for w in range(6)},
            {w: 1050 for w in range(6)},
            {w: 20 for w in range(6)},
            origin_week=0,
            n_weeks=6,
        )
        self.assertIsNone(ut95)  # GeoCap detected in 1 week when homogeneous UT95 couldn't in 6

    def test_end_to_end_benchmark_detection(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tpath = Path(td)
            inf_csv = tpath / "allevents.csv"
            ll_csv = tpath / "linelist.csv"
            geo_csv = tpath / "1S__geocap_samples.csv"
            conv_csv = tpath / "1S__convenience_samples.csv"
            out_sum = tpath / "Detection_Summary.csv"
            out_ev = tpath / "Detection_Events.csv"

            # 6 weeks starting 2021-07-01 (week 0 = 2021-07-01..2021-07-07, start_date = 2021-07-08)
            # Background component "1" (1000/wk in Fairfax metro, 50/wk in Montgomery rural)
            # Metro target "10" (starts week 1, 25/wk in Fairfax -> 125 total >= 100)
            # Rural target "20" (starts week 1, 25/wk in Montgomery -> 125 total >= 100)
            # Die-out component "99" (starts week 1, 5 total < 100 -> must be filtered out)
            inf_rows = []
            ll_rows = []
            geo_rows = []
            conv_rows = []
            pid = 0

            for wk in range(6):
                dt = f"2021-07-{2 + wk * 7:02d}" if (2 + wk * 7) <= 31 else f"2021-08-{(2 + wk * 7) - 31:02d}"
                # Background in Fairfax (RUCC 1)
                ffx_wk = []
                for j in range(1000):
                    pid += 1
                    row = {
                        "alias_pid": f"p{pid}",
                        "date": dt,
                        "sim_tick": str(100 + wk * 7),
                        "county": "Fairfax VA",
                        "rucc_code": "1",
                        "component_id": "1",
                    }
                    inf_rows.append(row)
                    ll_rows.append(row)
                    ffx_wk.append(row)
                # Duplicate alias_pid row to verify deduplication
                inf_rows.append(dict(ffx_wk[0], sim_tick=str(102 + wk * 7)))

                # Background in Montgomery (RUCC 6)
                mo_wk = []
                for j in range(40):
                    pid += 1
                    row = {
                        "alias_pid": f"p{pid}",
                        "date": dt,
                        "sim_tick": str(100 + wk * 7),
                        "county": "Montgomery VA",
                        "rucc_code": "6",
                        "component_id": "1",
                    }
                    inf_rows.append(row)
                    ll_rows.append(row)
                    mo_wk.append(row)

                if wk >= 1:
                    # Metro target "10" in Fairfax
                    for j in range(25):
                        pid += 1
                        row = {
                            "alias_pid": f"p{pid}",
                            "date": dt,
                            "sim_tick": str(100 + wk * 7),
                            "county": "Fairfax VA",
                            "rucc_code": "1",
                            "component_id": "10",
                        }
                        inf_rows.append(row)
                        ll_rows.append(row)
                        ffx_wk.append(row)

                    # Rural target "20" in Montgomery
                    for j in range(25):
                        pid += 1
                        row = {
                            "alias_pid": f"p{pid}",
                            "date": dt,
                            "sim_tick": str(100 + wk * 7),
                            "county": "Montgomery VA",
                            "rucc_code": "6",
                            "component_id": "20",
                        }
                        inf_rows.append(row)
                        ll_rows.append(row)
                        mo_wk.append(row)

                if wk == 1:
                    # Die-out component "99" (size 5 < 100)
                    for j in range(5):
                        pid += 1
                        row = {
                            "alias_pid": f"p{pid}",
                            "date": dt,
                            "sim_tick": str(100 + wk * 7),
                            "county": "Montgomery VA",
                            "rucc_code": "6",
                            "component_id": "99",
                        }
                        inf_rows.append(row)
                        ll_rows.append(row)

                # GeoCap draws 20 from Fairfax + 20 from Montgomery (40/wk total)
                geo_rows.extend(ffx_wk[-20:])
                geo_rows.extend(mo_wk[-20:])
                # Convenience draws 40 from Fairfax + 0 from Montgomery (40/wk total)
                conv_rows.extend(ffx_wk[-40:])

            cols = ["alias_pid", "date", "sim_tick", "county", "rucc_code", "component_id"]
            for path, rows in [
                (inf_csv, inf_rows),
                (ll_csv, ll_rows),
                (geo_csv, geo_rows),
                (conv_csv, conv_rows),
            ]:
                with open(path, "w", encoding="utf-8", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=cols)
                    w.writeheader()
                    w.writerows(rows)

            summary, events = benchmark_detection(
                samples_globs=[str(geo_csv), str(conv_csv)],
                infections=inf_csv,
                linelist=ll_csv,
                target_col="component_id",
                start_date="2021-07-08",
                min_pool=50,
                min_chain_size=100,
                stratum_cols=("county",),
                out_csv=out_sum,
                out_events_csv=out_ev,
            )

            self.assertTrue(out_sum.is_file())
            self.assertTrue(out_ev.is_file())

            # Only components "1", "10", "20" have >= 100 cases; "99" is excluded
            target_ids = {e["target_id"] for e in events}
            self.assertEqual(target_ids, {"1", "10", "20"})

            # Check rural component "20": GeoCap beats UT95 (delta_t95 < 0), Convenience fails to reach 95%
            ev_map = {(e["strategy"], e["target_id"]): e for e in events}
            geo_20 = ev_map[("1S__geocap", "20")]
            conv_20 = ev_map[("1S__convenience", "20")]
            self.assertEqual(geo_20["geo_group"], "rural")
            self.assertEqual(geo_20["t95_weeks"], 0)  # 20 draws out of 65 with 25 targets -> P(det) > 0.9999 in week 0!
            self.assertLess(float(geo_20["delta_t95_weeks"]), 0.0)
            self.assertEqual(conv_20["t95_weeks"], "")  # never detected
            # Right-censored delta_t95 penalizes missed rural chain instead of dropping it
            self.assertGreater(float(conv_20["delta_t95_censored_weeks"]), 0.0)

            # Summary contains all / metro / rural rows for each strategy
            sum_keys = {(r["strategy"], r["geo_group"]) for r in summary}
            self.assertIn(("1S__geocap", "all"), sum_keys)
            self.assertIn(("1S__geocap", "metro"), sum_keys)
            self.assertIn(("1S__geocap", "rural"), sum_keys)
            self.assertIn(("1S__convenience", "rural"), sum_keys)


if __name__ == "__main__":
    unittest.main()
