# Dump the contents of a Phoenix .TBL file to a json file
#
# Author: Taylor Viti
# Date Created: 07/14/2019

import phx_files as phx
import argparse
import json
import os

parser = argparse.ArgumentParser("Dump the contents of a .TBL file to a .json file")
parser.add_argument("tbl_fn", type=str)
parser.add_argument("json_fn", type=str)

args = parser.parse_args()
tbl_fn = args.tbl_fn
json_fn = args.json_fn

print("Input file: %s" % tbl_fn)

# If the output is a dir, then use the input filename as the output
if os.path.isdir(json_fn):
    # Remove any trailing /'s
    if json_fn[-1] == "/":
        json_fn = json_fn[0:-1]

    base_fn = tbl_fn.split("/")[-1]
    base_fn = base_fn.split(".")[-2]

    json_fn = "%s/%s.json" % (json_fn, base_fn)

print("Output file: %s" % json_fn)

tbl = phx.read_TBL(tbl_fn)
for k in tbl.keys():
    print(k, tbl[k])

with open(json_fn, "w") as f:
    json.dump(tbl, f)
