"""Read scores already written by AutoDock Vina or Amber MMPBSA.py."""

import argparse
import json
import os
import re
from pathlib import Path


def read_vina_result(out):
    if os.path.isfile(out):
        for l in open(out):
            m=re.match(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)",l)
            if m: return float(m.group(1))
    return None


def parse_delta(dat):
    txt=open(dat).read(); i=txt.find("Differences (Complex - Receptor - Ligand)")
    if i<0: return {}
    seg=txt[i:]
    def g(l):
        m=re.search(rf"\n{re.escape(l)}\s+(-?\d+\.\d+)", seg); return float(m.group(1)) if m else None
    return dict(deltaG_total=g("DELTA TOTAL"), dVDWAALS=g("VDWAALS"), dEEL=g("EEL"),
                dEGB=g("EGB"), dESURF=g("ESURF"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("format", choices=("vina", "mmgbsa"))
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    if not args.path.is_file():
        parser.error(f"File does not exist: {args.path}")
    if args.format == "vina":
        result = read_vina_result(args.path)
        if result is None:
            parser.error("No Vina score record found")
    else:
        result = parse_delta(args.path)
        if not result or all(value is None for value in result.values()):
            parser.error("No energy-difference values found")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
