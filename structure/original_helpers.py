"""Original study functions. See README.md for source line ranges and verification."""

import glob
import json
import os
import re
import time
import urllib.request

import numpy as np


def safe_name(s):
    return re.sub(r"[^a-z0-9_]","", s.lower().replace("-","_").replace(".","_"))[:60]


def job(name, seq, count):
    return {"name": safe_name(name), "modelSeeds": [],
            "sequences": [{"proteinChain": {"sequence": seq, "count": count}}]}


def esmfold(sequence, outpdb):
    if os.path.isfile(outpdb) and os.path.getsize(outpdb)>1000: return True
    for attempt in range(4):
        try:
            req=urllib.request.Request("https://api.esmatlas.com/foldSequence/v1/pdb/",
                                       data=sequence.encode(), method="POST")
            txt=urllib.request.urlopen(req, timeout=180).read().decode()
            if txt.startswith(("HEADER","ATOM","MODEL")):
                open(outpdb,"w").write(txt); return True
        except Exception as e:
            print(f"      esmfold retry {attempt+1}: {e}"); time.sleep(5)
    return False


def mean_plddt(pdb):
    b=[float(l[60:66]) for l in open(pdb) if l.startswith("ATOM") and l[12:16].strip()=="CA"]
    if not b: return 0.0
    m=float(np.mean(b))
    return round(m*100 if m<=1.5 else m, 1)   # ESMFold API returns pLDDT on 0-1; normalize to 0-100


def rank_of(path):
    """AF3 files are like fold_<job>_model_0.cif ; lower model index = higher rank."""
    m = re.search(r"model_(\d+)", os.path.basename(path).lower())
    return int(m.group(1)) if m else 99


def cif_to_pdb(cif, pdb):
    # 1) gemmi (best for mmCIF)
    try:
        import gemmi
        st = gemmi.read_structure(cif)
        st.setup_entities()
        st.write_pdb(pdb)
        if os.path.getsize(pdb) > 1000: return "gemmi"
    except Exception:
        pass
    # 2) Biopython
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        s = MMCIFParser(QUIET=True).get_structure("m", cif)
        io = PDBIO(); io.set_structure(s); io.save(pdb)
        if os.path.getsize(pdb) > 1000: return "biopython"
    except Exception:
        pass
    # 3) obabel
    try:
        import subprocess
        subprocess.run(["obabel", cif, "-O", pdb], capture_output=True, timeout=120)
        if os.path.isfile(pdb) and os.path.getsize(pdb) > 1000: return "obabel"
    except Exception:
        pass
    return None


def best_model_metrics(job_dir):
    """Return metrics of the best (highest ranking_score) of the 5 models."""
    summ = sorted(glob.glob(os.path.join(job_dir, "*summary_confidences_*.json")))
    best = None
    for f in summ:
        d = json.load(open(f))
        rs = d.get("ranking_score", -1)
        if best is None or rs > best["ranking_score"]:
            cp = d.get("chain_pair_iptm", [[None, None], [None, None]])
            iface = cp[0][1] if len(cp) > 1 and len(cp[0]) > 1 else None
            best = {"iptm": d.get("iptm"), "ptm": d.get("ptm"),
                    "ranking_score": rs, "interface_iptm": iface,
                    "n_models": len(summ)}
    return best
