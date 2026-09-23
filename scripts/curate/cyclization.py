from __future__ import annotations
import sys
from pathlib import Path
import gemmi
sys.path.insert(0, str(Path(__file__).resolve().parent))
import struct_utils

HEADTAIL_MAX = 1.9
DISULF_MAX = 2.5

def _ordered_aa_residues(chain):
    return [r for r in chain if struct_utils.aa_info(r)[1]]

def _is_cys(res):
    return struct_utils.aa_info(res)[3] == "CYS"

def _atom(res, name):
    for a in res:
        if a.name == name:
            return a
    return None

def _dist(a, b):
    return a.pos.dist(b.pos) if (a and b) else None

def detect_from_connections(st, chain_id):
    model = st[0]
    chain = model.find_chain(chain_id)
    order = {str(r.seqid): i for i, r in enumerate(chain)} if chain else {}
    n = len(order)
    out = {"disulf": False, "headtail": False, "sidechain": False,
           "bonds": [], "sidechain_bonds": []}
    for con in st.connections:
        p1, p2 = con.partner1, con.partner2
        if not (p1.chain_name == chain_id and p2.chain_name == chain_id):
            continue
        kind = str(con.type)
        i1, i2 = order.get(str(p1.res_id.seqid)), order.get(str(p2.res_id.seqid))
        atoms = {p1.atom_name, p2.atom_name}
        bond = {"atom1": p1.atom_name, "res1": str(p1.res_id.seqid),
                "atom2": p2.atom_name, "res2": str(p2.res_id.seqid), "kind": kind, "src": "conn"}
        out["bonds"].append(bond)
        if "Disulf" in kind or atoms == {"SG"}:
            out["disulf"] = True
        elif atoms == {"C", "N"} and i1 is not None and i2 is not None:
            if {i1, i2} == {0, n - 1} and n > 2:
                out["headtail"] = True
            elif abs(i1 - i2) > 1:
                out["sidechain"] = True; out["sidechain_bonds"].append(bond)
        elif i1 is not None and i2 is not None and abs(i1 - i2) > 1 and "Covale" in kind:
            out["sidechain"] = True; out["sidechain_bonds"].append(bond)
    return out

def subtype_sidechain(bonds):
    def elem(a):
        return a[0] if a else "?"
    seen = set()
    for b in bonds:
        e = {elem(b["atom1"]), elem(b["atom2"])}
        if "S" in e and "C" in e:
            seen.add("thioether")
        elif e == {"O", "C"}:
            seen.add("ester")
        elif e == {"N", "C"}:
            seen.add("lactam")
        elif e == {"C"}:
            seen.add("staple")
    for t in ("thioether", "ester", "lactam", "staple"):
        if t in seen:
            return t
    return "other"

def detect_geometric(st, chain_id):
    model = st[0]
    chain = model.find_chain(chain_id)
    out = {"headtail": False, "disulf": False, "n_c_dist": None, "bonds": []}
    if chain is None:
        return out
    aa = _ordered_aa_residues(chain)
    if len(aa) >= 3:
        n_first = _atom(aa[0], "N")
        c_last = _atom(aa[-1], "C")
        d = _dist(n_first, c_last)
        out["n_c_dist"] = round(d, 2) if d is not None else None
        if d is not None and d < HEADTAIL_MAX:
            out["headtail"] = True
            out["bonds"].append({"atom1": "N", "res1": str(aa[0].seqid),
                                 "atom2": "C", "res2": str(aa[-1].seqid),
                                 "kind": "head_to_tail", "src": "geom", "dist": round(d, 2)})
    sgs = [(r, _atom(r, "SG")) for r in aa if _is_cys(r) and _atom(r, "SG")]
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            d = _dist(sgs[i][1], sgs[j][1])
            if d is not None and d < DISULF_MAX:
                out["disulf"] = True
                out["bonds"].append({"atom1": "SG", "res1": str(sgs[i][0].seqid),
                                     "atom2": "SG", "res2": str(sgs[j][0].seqid),
                                     "kind": "disulfide", "src": "geom", "dist": round(d, 2)})
    return out

def classify(st, chain_id):
    conn = detect_from_connections(st, chain_id)
    geom = detect_geometric(st, chain_id)
    headtail = conn["headtail"] or geom["headtail"]
    disulf = conn["disulf"] or geom["disulf"]
    sidechain = conn["sidechain"]
    bonds = conn["bonds"] + geom["bonds"]
    if headtail and disulf:
        ctype = "bicyclic"
    elif headtail:
        ctype = "head_to_tail"
    elif disulf:
        ctype = "disulfide"
    elif sidechain:
        ctype = subtype_sidechain(conn["sidechain_bonds"])
    else:
        ctype = "linear"
    ev = f"conn(ht={conn['headtail']},ss={conn['disulf']},sc={conn['sidechain']})|" \
         f"geom(ht={geom['headtail']},ss={geom['disulf']},nc={geom['n_c_dist']})"
    return ctype, bonds, ev

if __name__ == "__main__":
    print("cyclization module ok; HEADTAIL_MAX=%.1f DISULF_MAX=%.1f" % (HEADTAIL_MAX, DISULF_MAX))
