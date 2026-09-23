#!/usr/bin/env python
import sys, os, collections, statistics

def parse(path):
    lines = open(path).read().splitlines()
    try:
        a0 = lines.index('@<TRIPOS>ATOM') + 1
    except ValueError:
        return lines, [], 0
    a1 = next((i for i in range(a0, len(lines)) if lines[i].startswith('@<TRIPOS>')), len(lines))
    return lines, list(range(a0, a1)), a0

def main(mol2dir):
    files = sorted(f for f in os.listdir(mol2dir) if f.endswith('.mol2'))
    types, charges = collections.defaultdict(list), collections.defaultdict(list)
    for f in files:
        lines, idxs, _ = parse(os.path.join(mol2dir, f))
        for i in idxs:
            p = lines[i].split()
            if len(p) < 9 or p[5] == 'Hev':
                continue
            key = (p[1], p[6])
            types[key].append(p[5])
            charges[key].append(float(p[8]))
    repaired = unrepairable = 0
    for f in files:
        path = os.path.join(mol2dir, f)
        lines, idxs, _ = parse(path)
        changed = False
        for i in idxs:
            p = lines[i].split()
            if len(p) < 9 or p[5] != 'Hev':
                continue
            key = (p[1], p[6])
            if not types[key]:
                unrepairable += 1
                continue
            t = collections.Counter(types[key]).most_common(1)[0][0]
            q = statistics.median(charges[key])
            lines[i] = '%7s %-8s %9s %9s %9s %-9s %s %-9s %9.4f' % (
                p[0], p[1], p[2], p[3], p[4], t, p[6], p[7], q)
            repaired += 1
            changed = True
        if changed:
            open(path, 'w').write('\n'.join(lines) + '\n')
    print(f'repaired {repaired} Hev atoms, {unrepairable} with no consensus')
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
