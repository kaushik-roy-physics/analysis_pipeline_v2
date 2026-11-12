"""
The 2026-11 additions to analysis_pipeline_v2, end to end on a synthetic
confocal experiment with known answers (tests/synthetic_fields.py): the
replicate unit, permutation tests with their p floor, signaling readouts,
saturated nuclei, markers left out by name, skipped analyses, the 2-D / 3-D
folder name, numeric condition order, stale checkpoints, the comment-keeping GUI
writers and the HPC one-off settings.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
for p in (ROOT / "src", ROOT / "tests", ROOT / "gui", ROOT):
    sys.path.insert(0, str(p))

import synthetic_fields as sf  # noqa: E402

_CACHE: dict = {}


def _run(cfg: Path, *extra) -> int:
    import main as pipeline_main
    args = ["--config", str(cfg), "--log-level", "WARNING"]
    assert pipeline_main.main(["manifest"] + args) == 0
    code = pipeline_main.main(["downstream"] + args + list(extra))
    logging.getLogger().handlers.clear()
    return code


def _experiment():
    """One full run shared by the end-to-end checks (the signaling step is slow)."""
    if "root" not in _CACHE:
        root = Path(tempfile.mkdtemp())
        cfg = sf.make_experiment(root)
        assert _run(cfg) == 0
        _CACHE.update(root=root, cfg=cfg)
    return _CACHE["root"], _CACHE["cfg"]


def _table(name):
    root, _ = _experiment()
    return pd.read_csv(root / "analysis_2d" / "tables" / name)


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------

def test_known_answers_and_well_replicates():
    root, _ = _experiment()
    assert (root / "analysis_2d" / "report.html").exists()        # results_dir: auto
    unit = json.loads((root / "analysis_2d" / "tables" / "replicate_unit.json").read_text())
    assert unit["unit"] == "well"                                 # 3 wells per dose
    c = _table("conditions.csv").set_index(["marker", "condition"])
    for dose, truth in sf.POS_FRAC.items():
        got = c.loc[("SOX2", dose), "positive_fraction"]
        assert abs(got - truth) < 0.08, (dose, got, truth)
        assert c.loc[("SOX2", dose), "n_replicates"] == 3
    comp = _table("comparisons_positive_fraction.csv")
    sox = comp[comp["marker"] == "SOX2"]
    assert set(sox["unit"]) == {"well"} and sox["exact"].all()
    assert (sox["smallest_possible_p"].round(6) == 0.1).all()     # 3 vs 3 wells: 2/20
    assert (sox["p_value"] >= sox["smallest_possible_p"] - 1e-9).all()
    print("ok  % positive matches the truth; n = wells; exact permutation p with its floor")


def test_signaling_relative_to_reference_and_never_gated():
    c = _table("conditions.csv").set_index(["marker", "condition"])
    for dose, truth in sf.BCAT.items():
        got = c.loc[("b-CAT", dose), "median"]
        assert abs(got - truth) / truth < 0.1, (dose, got, truth)    # N:C ÷ control
    thr = _table("thresholds.csv")
    assert "b-CAT" not in set(thr["marker"])                       # not gated
    ref = _table("signaling_reference.csv")
    assert ref.loc[ref["marker"] == "b-CAT", "relative"].all()
    root, _ = _experiment()
    assert list((root / "analysis_2d" / "figures" / "signaling").glob("ring_qc_*.png"))
    print("ok  b-CAT N:C ÷ control = 1 / 2 / 4 (truth), no gate, ring QC pictures")


def test_saturated_nucleus_and_excluded_marker():
    sat = _table("qc_saturated_nuclei.csv")
    assert len(sat) == 1 and sat["channel_name"].iloc[0] == "SOX2"
    root, _ = _experiment()
    inten = pd.read_parquet(root / "analysis_2d" / "intensities.parquet")
    assert set(inten["marker"].astype(str)) == {"SOX2", "b-CAT"}  # ECAD and DAPI out
    html = (root / "analysis_2d" / "report.html").read_text()
    assert "Replicate = well" in html and "detector ceiling" in html
    print("ok  the clipped nucleus is listed and left out; ECAD excluded by name")


def test_effects_find_the_dose_trend():
    tr = _table("effects_trends.csv")
    row = tr[(tr["subject"] == "SOX2") & (tr["metric"] == "positive_fraction")].iloc[0]
    assert row["slope"] > 0 and row["slope_ci_lo"] > 0, row
    var = _table("effects_variance.csv")
    assert {"sd_field", "sd_well"} <= set(var.columns)
    print(f"ok  effects: SOX2 % positive rises with dose (slope {row['slope']:.2f}/decade)")


def test_skip_steps_set_aside():
    root = Path(tempfile.mkdtemp())
    src, _ = _experiment()
    shutil.copytree(src, root, dirs_exist_ok=True)
    cfg = root / "config.yaml"
    d = yaml.safe_load(cfg.read_text())
    d["downstream"]["skip_steps"] = ["joint", "plate"]
    cfg.write_text(yaml.safe_dump(d, sort_keys=False))
    assert _run(cfg) == 0
    out = root / "analysis_2d"
    assert not (out / "figures" / "joint").exists() and not (out / "figures" / "plate").exists()
    assert (out / "_skipped" / "joint" / "figures" / "joint").exists()
    assert "Joint marker views" not in (out / "report.html").read_text()
    from downstream import steps
    try:
        steps.resolve_skips(["thresholds"])
        raise AssertionError("skipping a core step must fail")
    except ValueError:
        pass
    print("ok  skip_steps: joint and plate left out, earlier outputs set aside")


# ---------------------------------------------------------------------------
# Units of replication without running the pipeline
# ---------------------------------------------------------------------------

def test_one_well_per_condition_falls_back_to_fields_with_a_caveat():
    from downstream import aggregate
    fields = pd.DataFrame({"field_id": range(6), "well": [1, 1, 1, 2, 2, 2],
                           "plate_id": "p", "condition": ["a"] * 3 + ["b"] * 3})
    unit, note = aggregate.choose_unit(fields, {})
    assert unit == "field" and "single well" in note
    fields["exp"] = ["d1", "d2", "d3"] * 2
    unit, _ = aggregate.choose_unit(fields, {"stats": {"experiment_column": "exp"}})
    assert unit == "experiment"
    print("ok  replicate: one well per condition → fields (said so); differentiations → experiment")


def test_numeric_condition_order():
    from downstream.plots.base import ordered_conditions
    got = ordered_conditions(["0.1", "100", "5", "10", "0.5"], {})
    assert got == ["0.1", "0.5", "5", "10", "100"], got
    got = ordered_conditions(["Activin | 100", "Activin | 5"], {})
    assert got == ["Activin | 5", "Activin | 100"]
    print("ok  doses sort as numbers on every axis")


# ---------------------------------------------------------------------------
# Shared with the micropattern pipeline (same code)
# ---------------------------------------------------------------------------

def test_stale_extraction_after_new_masks():
    from image_processing import intensity_extractor as ie
    root, _ = _experiment()
    res = next(p for p in (root / "zmax").iterdir() if p.is_dir())
    tmp = Path(tempfile.mkdtemp()) / res.name
    shutil.copytree(res, tmp)
    assert ie._extraction_current(tmp)
    time.sleep(1.1)
    os.utime(tmp / "nuclear_masks.tif")
    assert not ie._extraction_current(tmp)
    print("ok  new masks make the intensities stale")


def test_gui_writers_keep_comments():
    import core
    t = core.ConfigTemplate()
    existing = t.render(t.merged({})).replace(
        "  skip_steps: []", "  skip_steps: []  # MY NOTE", 1)
    vals = t.merged(yaml.safe_load(existing))
    assert t.render_into(vals, existing) == existing
    vals["downstream"]["stats"]["replicate"] = "well"
    out = t.render_into(vals, existing)
    assert "MY NOTE" in out and yaml.safe_load(out)["downstream"]["stats"]["replicate"] == "well"
    text = ("plate: 18well\nplate_id: p\nmode: predefined\npoints_per_well: 2\n"
            "stains:\n  default: {0: DAPI, 1: SOX2}\nwells:\n  # re-seeded\n"
            "  1-2: {dose: 0}   # ctrl\n  3: {dose: 5}\n")
    st = core.plate_map_from_yaml(yaml.safe_load(text))
    assert core.plate_map_text(st, "x", text) == text
    st["wells"][-1]["dose"] = "10"
    out = core.plate_map_text(st, "x", text)
    assert "# re-seeded" in out and "# ctrl" in out and yaml.safe_load(out)["wells"][3]["dose"] == 10
    print("ok  GUI saves into the files on disk, comments kept")


def test_hpc_one_off_overrides_and_dry_run():
    exp = Path(tempfile.mkdtemp())
    (exp / "config.yaml").write_text("datasets: []\n")
    env = (ROOT / "templates" / "hpc_job.env").read_text()
    env = env.replace('PIPELINE_DIR="/data/Warmflash/KR/analysis_pipeline_v2"',
                      f'PIPELINE_DIR="{ROOT}"')
    env = env.replace('VENV_DIR="/scratch/$USER/confocal_v2_venv"', f'VENV_DIR="{exp}/novenv"')
    (exp / "hpc_job.env").write_text(env)
    base = {k: v for k, v in os.environ.items() if not k.startswith("MP_")}
    r = subprocess.run(["bash", str(ROOT / "hpc" / "run.sh"), "test"], cwd=exp,
                       capture_output=True, text=True, env={**base, "SEQ": "1-3"}, timeout=120)
    assert r.returncode == 0, r.stderr
    assert r.stdout.count("--seq 1-3") >= 2, r.stdout
    print("ok  HPC: run.sh test passes one-off settings to every job")


def test_power_rows_flag_the_permutation_floor():
    """MDD rows say when a permutation test cannot reach alpha (E vs E: floor 2 / C(2E, E))."""
    from downstream import effects as ef_mod
    assert abs(ef_mod.permutation_floor(3) - 0.1) < 1e-12
    assert abs(ef_mod.permutation_floor(4) - 2 / 70) < 1e-12
    prm = ef_mod.settings({})
    prm["power"]["experiments"] = [3, 4, 6]
    vc = {"sd_experiment": 8.0, "sd_well": 5.0, "sd_field": 10.0}
    rows = {r["experiments"]: r for r in ef_mod.power_table(vc, prm, 3, 5, 50.0)}
    assert not rows[3]["permutation_can_reach_alpha"]
    assert rows[4]["permutation_can_reach_alpha"] and rows[6]["permutation_can_reach_alpha"]
    assert abs(rows[3]["mdd"] - 27.0) < 0.1 and abs(rows[6]["mdd"] - 15.9) < 0.1
    print("ok  power: 3 vs 3 flagged (p floor 0.1); MDD 27.0 / 15.9 as in the compendium")


if __name__ == "__main__":
    failures = 0
    for name, func in list(globals().items()):
        if name.startswith("test_") and callable(func):
            try:
                func()
            except Exception as exc:                          # noqa: BLE001
                import traceback
                traceback.print_exc()
                failures += 1
                print(f"FAIL {name}: {exc!r}")
    print("\n" + ("all tests passed" if not failures else f"{failures} test(s) failed"))
    sys.exit(1 if failures else 0)
