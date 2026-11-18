#!/usr/bin/env python3
"""
Check that a uniform hollow cylinder (rin < r < rout) is loaded with the
correct total weight, and only inside the plasma, for every value of
radial_numpercell_power.
"""

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

n0 = 1.0e24
rin = 0.5e-6
rout = 1.0e-6
Lz = 1.0e-6
expected = n0 * np.pi * (rout**2 - rin**2) * Lz
tolerance = 0.01

ts = OpenPMDTimeSeries("./diags/diag1")
for species in ("power0", "power1", "power2"):
    x, y, w = ts.get_particle(["x", "y", "w"], species=species, iteration=0)
    r = np.sqrt(x**2 + y**2)
    rel_err = abs(w.sum() - expected) / expected
    print(
        f"{species}: total weight = {w.sum():.6e}, expected = {expected:.6e}, "
        f"relative error = {rel_err:.2e}, min r/rin = {r.min() / rin:.4f}"
    )
    assert rel_err < tolerance
    assert r.min() >= rin * (1.0 - 1.0e-12)
