#!/usr/bin/env python
import phx_files as phx
import numpy as np
import argparse
import json
import h5py
import os

desc = """
Pack a set of .TSn files and associated metadata, to an .HDF5 file.
The output .HDF5 will have the same base file name as the base name for
the input .TBL and .TSn files. This program is expecting a .TS2, .TS3,
and .TS4 file with the same base name as tbl_fn to exist in ts_dir.

This program was originally written to package data produced by an MTU-5A.

WARNING: If the HDF5 file already exists, this program will OVERWRITE IT!

If you have access to Phoenix's proprietary SysCal.exe program, and have
created an appropriate .PFC file for running it, then this program can
also use SysCal.exe to generate a system response for all channels on
instrument (meaning the box and sensor responses, including dipole
lengths), based on the frequencies provided in the .PFC file, and then
will pack the calibration into the .HDF file as well.

Data packed to HDF5:
    tbl: .JSON formatted copy of .TBL file contents
    TS2: 24 kHz samples (containing tags and recs)
    TS3: 2.4 kHz samples
    TS4: 150 Hz samples (recorded CONTINUOUSLY)
Optional:
    cal: N_freqs x (5 + 2) complex array of system response functions.
       First column is frequency (Hz), second column is the associated
       sampling rate level (i.e. the n in TSn). Channels 1 and 2 have units
       of 1/(V/m), channels 3, 4, and 5 units of 1/T. Timeseries data must
       be transformed to unit scale prior to application (divide by 2**23).

Author: Taylor Viti
Date Created: 08/05/2019
"""

parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse
                                 .RawDescriptionHelpFormatter)
parser.add_argument("tbl_fn", type=str,
                    help="Path to the table file for this deployment")
parser.add_argument("ts_dir", type=str,
                    help="""Path to the directory containing the
                    associated .TSn files""")
parser.add_argument("HDF5_out_dir", type=str,
                    help="Path to dump the resultant .HDF5 file")
parser.add_argument("--syscal_fn", type=str, default=None,
                    help="Path to SysCal.exe")
parser.add_argument("--cts_out_dir", type=str, default="/tmp/",
                    help="""Location to dump ascii cal data
                    (defaults to /tmp)""")
parser.add_argument("--pfc_fn", type=str, default=None,
                    help="Location of SysCal parameter file")
parser.add_argument("--clb_dir", type=str, default=None,
                    help="Path to dir containing box cal files")
parser.add_argument("--clc_dir", type=str, default=None,
                    help="Path to dir containing mag cal files")

args = parser.parse_args()
tbl_fn = args.tbl_fn
ts_dir = args.ts_dir
cts_dir = args.cts_out_dir
out_dir = args.HDF5_out_dir
pfc_fn = args.pfc_fn
clb_dir = args.clb_dir
clc_dir = args.clc_dir
syscal_fn = args.syscal_fn

base_fn, ext = os.path.splitext(tbl_fn)
base_fn = base_fn.split("/")[-1]

# Read in the table file
print("Reading .TBL file %s" % tbl_fn)
tbl, grps, smphs, types = phx.read_TBL(tbl_fn)

# Generate a calibration, if asked
if args.syscal_fn is not None:
    # Read in the freqs from the PFC file
    pfc = np.loadtxt(pfc_fn,
                     delimiter=",",
                     skiprows=1,
                     usecols=(1, 2))
    levels = pfc[:, 0]
    freqs = pfc[:, 1]

    # Generate the cal data, and save the resultant cts to a file too
    cts_fn = "%s/%s.cts" % (cts_dir, base_fn)
    print("Creating calibration data: %s" % cts_fn)
    cts = phx.get_syscal(freqs, tbl_fn, clb_dir, clc_dir,
                         cts_fname=cts_fn,
                         SYSCAL_fname=syscal_fn,
                         levels=levels)

# Read in the .TS2, .TS3, and .TS4 data
TS2_fn = "%s/%s.TS2" % (ts_dir, base_fn)
TS3_fn = "%s/%s.TS3" % (ts_dir, base_fn)
TS4_fn = "%s/%s.TS4" % (ts_dir, base_fn)

print("Reading in TSn files...")
print("TS2: %s" % TS2_fn)
print("TS3: %s" % TS3_fn)
print("TS4: %s" % TS4_fn)
tags2, recs2 = phx.read_TSn(TS2_fn)
tags3, recs3 = phx.read_TSn(TS3_fn)
tags4, recs4 = phx.read_TSn(TS4_fn)

# Pack everything into an HDF5 file
out_fn = "%s/%s.HDF5" % (out_dir, base_fn)

# Delete the file if it exists, so we can make a new one
if os.path.exists(out_fn):
    os.remove(out_fn)

print("Writing to %s" % out_fn)
with h5py.File(out_fn, "w") as f:
    # Pack the deployment metadata
    f.create_dataset("tbl", data=json.dumps(tbl))

    # Pack the cal, if asked
    if args.syscal_fn is not None:
        f.create_dataset("cal", data=cts)

    # Pack the timeseries data
    TS2 = f.create_group("TS2")
    TS3 = f.create_group("TS3")
    TS4 = f.create_group("TS4")

    TS2.create_dataset("tags", data=tags2)
    TS2.create_dataset("recs", data=recs2)

    TS3.create_dataset("tags", data=tags3)
    TS3.create_dataset("recs", data=recs3)

    TS4.create_dataset("tags", data=tags4)
    TS4.create_dataset("recs", data=recs4)
