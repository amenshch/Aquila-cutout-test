#!/usr/bin/env python
"""inject_config_v3_flat_iso.py -- isolated control, background flattened.

Run with:
    python -B run_inject_v3.py inject_config_v3_iso

Identical to inject_config_v3_flat.py in every setting, including RANDOM_SEED,
except that CONTROL_ISOLATED is on.  The placement is carried out in full and
unchanged; only afterwards are the cores that overlap a kept one removed,
leaving a maximal non-overlapping subset in placement order.

Every core that survives sits at the same position, with the same model, on the
same background, as in the companion run.  The comparison of the two therefore
isolates the effect of blending between injected cores, core by core, in the
same way the flattened pair isolates the effect of the cloud fluctuations.

A control produced instead by re-running with a stricter separation rule would
place different models in different places, and the comparison would confound
blending with sampling.  That is not a hypothetical risk: comparing the
flattened and unflattened runs sample against sample suggested the cloud
fluctuations biased the recovered mass high by 39%, and the effect disappeared
once the same cores were compared in both.

Output files are named cmf<slope>_n<cores>_s<seed>_iso, so they cannot
overwrite the companion's.  Verify the pairing with

    diff <(grep -v '^#' inj_..._iso_truth.txt | awk '{print $2,$3,$4}') \\
         <(grep -v '^#' inj_..._truth.txt     | awk '{print $2,$3,$4}')

which should list only the removed cores, never a changed position.

Four fields can be made from one placement, by combining this switch with
FLATTEN_BACKGROUND: as placed, as placed with the background flattened,
isolated, and isolated with the background flattened.  Every core appears in
all four wherever it survives the filter, so the blending and fluctuation
contributions can be separated simultaneously.
"""
import os

_here = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(_here, 'inject_config_v3_flat.py')).read())

CONTROL_ISOLATED = True
