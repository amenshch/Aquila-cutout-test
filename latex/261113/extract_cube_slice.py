#!/usr/bin/env python3
"""extract_cube_slice.py -- a vertical slice of the RADMC-3D diffuse cloud, for the figure of the paper.

The slice is the plane that contains one row of the field and the line of sight: the horizontal coordinate x across the field and
the depth z along the line of sight. For every cell of that plane it saves the H2 number density and the dust temperature. The
margins of the cube are left out. The cubes are read through memory maps, so the script needs little memory even for a cube of
12 arcsec cells.

Run it in the directory of the RADMC-3D run (with cirrus3d.json, amr_grid.inp, dust_density.binp and dust_temperature.bdat or
dust_temperature.dat):

    python3 extract_cube_slice.py [row]
    python3 extract_cube_slice.py --column [column]

The first form cuts along a row, the second along a column of the field. Rows and columns are counted from the bottom and from the
left of the field, margins excluded; the defaults are the middle row and the last (rightmost) column. The output is the compressed
file cirrus3d_slice_row<row>.npz or cirrus3d_slice_col<column>.npz, about 0.5 MB, with the arrays n (cm^-3) and T (K), both of
shape (nz, number of cells along the cut), the surface density along the cut (cm^-2), the cell size (pc), and the parameters of the
cube.
"""

import json
import os
import sys

import numpy as np

AMU = 1.66053907e-24


def open_binary(name, shape):
    with open(name, "rb") as f:
        hdr = np.fromfile(f, dtype=np.int64, count=4)
    return np.memmap(name, dtype=np.float64 if hdr[1] == 8 else np.float32, mode="r", offset=32, shape=shape)


def main():
    par = json.load(open("cirrus3d.json"))
    nx, ny, nz, m = par["nx"], par["ny"], par["nz"], par["margin_cells"]
    fy, fx = ny - 2 * m, nx - 2 * m
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    bycol = "--column" in sys.argv
    if bycol:
        idx = int(args[0]) if args else fx - 1
        if not 0 <= idx < fx:
            sys.exit(" the column must be between 0 and %d" % (fx - 1))
        cut = (slice(None), slice(m, ny - m), m + idx)
    else:
        idx = int(args[0]) if args else fy // 2
        if not 0 <= idx < fy:
            sys.exit(" the row must be between 0 and %d" % (fy - 1))
        cut = (slice(None), m + idx, slice(m, nx - m))
    rho = open_binary("dust_density.binp", (nz, ny, nx))
    if os.path.exists("dust_temperature.bdat"):
        temp = np.array(open_binary("dust_temperature.bdat", (nz, ny, nx))[cut])
    else:
        with open("dust_temperature.dat") as f:
            f.readline(), f.readline(), f.readline()
            temp = np.array(f.read().split(), float)[:nx * ny * nz].reshape(nz, ny, nx)[cut]
    n = np.array(rho[cut]) / (par["dust2gas"] * par["muh2"] * AMU)
    dz = par["cell_pc"] * 3.0856775814913673e18
    sigma = n.sum(axis=0) * dz
    kind = "col" if bycol else "row"
    out = "cirrus3d_slice_%s%d.npz" % (kind, idx)
    np.savez_compressed(out, n=n.astype(np.float32), T=temp.astype(np.float32), sigma=sigma.astype(np.float32),
                        cell_pc=par["cell_pc"], row=(-1 if bycol else idx), col=(idx if bycol else -1), params=json.dumps(par))
    print(" %s %d: slice of %d x %d cells, surface density %.2e to %.2e cm-2, temperature %.2f to %.2f K"
          % ("column" if bycol else "row", idx, n.shape[0], n.shape[1], sigma.min(), sigma.max(), temp[temp > 0].min(),
             temp.max()))
    print(" written %s (%.2f MB)" % (out, os.path.getsize(out) / 1e6))


if __name__ == "__main__":
    main()
