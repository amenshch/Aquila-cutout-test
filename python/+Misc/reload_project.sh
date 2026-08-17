#!/usr/bin/env bash
# Rebuild the working state for the RT mass-correction project in a fresh chat.
# The container filesystem is wiped between chats; this refetches everything from
# GitHub (amenshch/Aquila-cutout-test, branch main).  Run from /home/claude/proj.
#
#   bash reload_project.sh
#
# Prefer any newer files the user uploads to /mnt/user-data/uploads over these.

set -e
BASE="https://raw.githubusercontent.com/amenshch/Aquila-cutout-test/main"
mkdir -p cats python latex

echo "== catalogues =="
for f in \
  cats/bes_grid_final2.txt \
  cats/bes_model_grid_final2_catalog \
  cats/bes_model_grid_final2_recovery_tables_v2 \
  cats/bes_model_grid_final2_catalog_rec \
  cats/bes_model_params_subdivision_catalog \
  cats/bes_model_params_subdivision_catalog_rec ; do
  curl -s -o "$f" -w "[%{http_code}] $f\n" -L "$BASE/$f" || true
done

echo "== code =="
for f in \
  python/fit_mass.py \
  python/tabulate_recovery.py \
  python/getsf_catalog.py \
  python/add_recoverable_mass.py \
  python/calibrate_floors.py ; do
  curl -s -o "$f" -w "[%{http_code}] $f\n" -L "$BASE/$f" || true
done

echo "== draft =="
for f in latex/rt_mass_correction.tex latex/rt_mass_correction.bib ; do
  curl -s -o "$f" -w "[%{http_code}] $f\n" -L "$BASE/$f" || true
done

echo "== per-cloud getsf catalogues (Guoyin) =="
# name pattern: <Cloud>-Guoyin/<Cloud>.sw.sources.ok.cat=<Cloud>.sw.sources.ok.add.cat=thin.<Cloud>.sw.sources.ok.00.cat
# (W3/W4/W5 uses prefix "TriRegion").  Fetch on demand — large; uncomment as needed.
# for C in Aquila OrionA California Ophiuchus Scorpius CygnusX ; do
#   d="${C}-Guoyin"; mkdir -p "$d"
#   n="${C}.sw.sources.ok.cat=${C}.sw.sources.ok.add.cat=thin.${C}.sw.sources.ok.00.cat"
#   curl -s -o "$d/$n" -w "[%{http_code}] $d\n" -L "$BASE/$d/$n" || true
# done

echo
echo "Reload done.  Distances: Aquila 260, OrionA 432, California 470,"
echo "Ophiuchus 144, Scorpius 130, CygnusX 1150, W3/W4/W5 (TriRegion) 1706 pc."
echo "Primary corrector: CMInvariantCorrector (fit_mass.py) on Σ, conc_peakmean, R_foot."
