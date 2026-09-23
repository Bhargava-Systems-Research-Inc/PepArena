import sys, math
inp, out = sys.argv[1], sys.argv[2]
atoms = []
for l in open(inp):
    if not l.startswith("ATOM"): continue
    el = (l[76:78].strip() or l[12:16].strip()[:1])
    if el == "H": continue
    atoms.append(l.rstrip("\n"))
sg = [(float(a[30:38]), float(a[38:46]), float(a[46:54]), a[22:27].strip())
      for a in atoms if a[12:16].strip() == "SG"]
cyx = set()
for i in range(len(sg)):
    for j in range(i + 1, len(sg)):
        if math.dist(sg[i][:3], sg[j][:3]) < 2.5:
            cyx.add(sg[i][3]); cyx.add(sg[j][3])
with open(out, "w") as w:
    for a in atoms:
        if a[17:20].strip() == "CYS" and a[22:27].strip() in cyx:
            a = a[:17] + "CYX" + a[20:]
        w.write(a + "\n")
    w.write("TER\n")
print(f"{out}: atoms={len(atoms)} cyx={len(cyx)}")
