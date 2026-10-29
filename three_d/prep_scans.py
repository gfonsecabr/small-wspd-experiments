"""Converts the Stanford scans (Wavefront .obj files of the mirror github.com/alecjacobson/common-3d-test-models) into
integer point files for exp3d: data/scans/NAME.xyz, one "x y z" line per distinct vertex, scaled so that the
largest extent of the bounding box is 10^6 and rounded. Run from three_d/ (downloads the .obj files if missing)."""
import os
import urllib.request

URL = 'https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/master/data/{}.obj'
MODELS = {'bunny': 'stanford-bunny', 'armadillo': 'armadillo', 'dragon': 'xyzrgb_dragon', 'buddha': 'happy'}
D = 'data/scans'

os.makedirs(D, exist_ok=True)
for name, src in MODELS.items():
    obj = f'{D}/{src}.obj'
    if not os.path.exists(obj):
        urllib.request.urlretrieve(URL.format(src), obj)
    V = [tuple(map(float, l.split()[1:4])) for l in open(obj) if l.startswith('v ')]
    lo = [min(v[k] for v in V) for k in range(3)]
    ext = max(max(v[k] for v in V) - lo[k] for k in range(3))
    P = []
    seen = set()
    for v in V:
        p = tuple(round((v[k] - lo[k]) / ext * 1e6) for k in range(3))
        if p not in seen:
            seen.add(p)
            P.append(p)
    with open(f'{D}/{name}.xyz', 'w') as f:
        f.writelines(f'{x} {y} {z}\n' for x, y, z in P)
    print(f'{name}: {len(V)} vertices, {len(P)} distinct points')
