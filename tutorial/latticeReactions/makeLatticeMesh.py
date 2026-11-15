#!/usr/bin/env python3
"""Stair-step hexahedral mesh of a Fischer-Koch S lattice in a channel, periodic in x and y.

Writes constant/polyMesh (points, faces, owner, neighbour, boundary) directly:
uniform cubes of size h, the cells inside the solid are left out. The side patches
left/right and bottom/top are conformal cyclic patches.

    python3 makeLatticeMesh.py [h in mm, default 0.2]

Geometry (mm): cross-section 10 x 10 (two unit cells of 5 mm), lattice from z = 0 to 10,
inlet section from z = -4, outlet section to z = 16. Solid where
    cos(2X)sin(Y)cos(Z) + cos(X)cos(2Y)sin(Z) + sin(X)cos(Y)cos(2Z) > T,   X = 2 pi x / L
"""
import sys, os
import numpy as np

h_mm = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
L, T = 5.0e-3, 0.35                     # unit cell [m], level of the surface
x0, x1, z0, z1 = -5.0e-3, 5.0e-3, -4.0e-3, 16.0e-3
zl0, zl1 = 0.0, 10.0e-3                 # lattice
h = h_mm*1e-3
nx = ny = int(round((x1 - x0)/h)); nz = int(round((z1 - z0)/h))
assert abs(nx*h - (x1 - x0)) < 1e-12 and abs(nz*h - (z1 - z0)) < 1e-12, "h must divide 10 mm and 20 mm"

xc = x0 + (np.arange(nx) + 0.5)*h; zc = z0 + (np.arange(nz) + 0.5)*h
X, Y, Z = np.meshgrid(xc, xc, zc, indexing='ij')            # [i, j, k]
a = 2*np.pi/L
F = np.cos(2*a*X)*np.sin(a*Y)*np.cos(a*Z) + np.cos(a*X)*np.cos(2*a*Y)*np.sin(a*Z) + np.sin(a*X)*np.cos(a*Y)*np.cos(2*a*Z)
fluid = ~((F > T) & (Z > zl0) & (Z < zl1))

# keep the fluid that is connected to the inlet (periodic in x and y)
reach = np.zeros_like(fluid); reach[:, :, 0] = fluid[:, :, 0]
while True:
    g = reach.copy()
    g |= np.roll(reach, 1, 0) | np.roll(reach, -1, 0) | np.roll(reach, 1, 1) | np.roll(reach, -1, 1)
    g[:, :, 1:] |= reach[:, :, :-1]; g[:, :, :-1] |= reach[:, :, 1:]
    g &= fluid
    if (g == reach).all(): break
    reach = g
print("fluid cells not connected to the inlet, taken as solid:", int(fluid.sum() - reach.sum()))
fluid = reach

# cell numbers, k-major
cid = -np.ones((nx, ny, nz), dtype=np.int64)
order = np.argwhere(fluid.transpose(2, 1, 0))               # rows (k, j, i) sorted
kk, jj, ii = order[:, 0], order[:, 1], order[:, 2]
ncell = len(order); cid[ii, jj, kk] = np.arange(ncell)

def pid(i, j, k): return (k*(ny + 1) + j)*(nx + 1) + i     # grid point number

def quad(i, j, k, side):
    if side == 'x-': return np.stack([pid(i, j, k), pid(i, j, k+1), pid(i, j+1, k+1), pid(i, j+1, k)], 1)
    if side == 'x+': return np.stack([pid(i+1, j, k), pid(i+1, j+1, k), pid(i+1, j+1, k+1), pid(i+1, j, k+1)], 1)
    if side == 'y-': return np.stack([pid(i, j, k), pid(i+1, j, k), pid(i+1, j, k+1), pid(i, j, k+1)], 1)
    if side == 'y+': return np.stack([pid(i, j+1, k), pid(i, j+1, k+1), pid(i+1, j+1, k+1), pid(i+1, j+1, k)], 1)
    if side == 'z-': return np.stack([pid(i, j, k), pid(i, j+1, k), pid(i+1, j+1, k), pid(i+1, j, k)], 1)
    if side == 'z+': return np.stack([pid(i, j, k+1), pid(i+1, j, k+1), pid(i+1, j+1, k+1), pid(i, j+1, k+1)], 1)

def nb(di, dj, dk):
    """cell number of the neighbour, -1 solid, -2 outside in z, periodic wrap flag in x/y"""
    i2, j2, k2 = ii + di, jj + dj, kk + dk
    wrap = (i2 < 0) | (i2 >= nx) | (j2 < 0) | (j2 >= ny)
    out = (k2 < 0) | (k2 >= nz)
    n = cid[i2 % nx, j2 % ny, np.clip(k2, 0, nz - 1)]
    n = np.where(out, -2, n)
    return n, wrap

c = np.arange(ncell)
# internal faces: owner c, neighbour with the higher number, in the order +x, +y, +z
fo, fn, fq, rank = [], [], [], []
for r, (d, side) in enumerate((((1, 0, 0), 'x+'), ((0, 1, 0), 'y+'), ((0, 0, 1), 'z+'))):
    n, wrap = nb(*d); m = (n >= 0) & ~wrap
    fo.append(c[m]); fn.append(n[m]); fq.append(quad(ii[m], jj[m], kk[m], side)); rank.append(np.full(m.sum(), r))
fo, fn, fq, rank = map(np.concatenate, (fo, fn, fq, rank))
s = np.lexsort((fn, fo)); fo, fn, fq = fo[s], fn[s], fq[s]
nint = len(fo)

patches = []          # (name, type, owner cells, quads, extra entries)
def collect(d, side, sel):
    n, wrap = nb(*d); m = sel(n, wrap)
    return c[m], quad(ii[m], jj[m], kk[m], side)

patches.append(('inlet', 'patch') + collect((0, 0, -1), 'z-', lambda n, w: n == -2) + ('',))
patches.append(('outlet', 'patch') + collect((0, 0, 1), 'z+', lambda n, w: n == -2) + ('',))
cyc = lambda n, w: (n >= 0) & w
sep = x1 - x0
for name, d, side, nbr, sv in (('left', (-1, 0, 0), 'x-', 'right', (-sep, 0, 0)), ('right', (1, 0, 0), 'x+', 'left', (sep, 0, 0)),
                               ('bottom', (0, -1, 0), 'y-', 'top', (0, -sep, 0)), ('top', (0, 1, 0), 'y+', 'bottom', (0, sep, 0))):
    extra = "        neighbourPatch  %s;\n        transformType   translational;\n        separation      (%g %g %g);\n" % ((nbr,) + sv)
    patches.append((name, 'cyclic') + collect(d, side, cyc) + (extra,))
wo, wq = [], []
for d, side in (((-1, 0, 0), 'x-'), ((1, 0, 0), 'x+'), ((0, -1, 0), 'y-'), ((0, 1, 0), 'y+'), ((0, 0, -1), 'z-'), ((0, 0, 1), 'z+')):
    o, q = collect(d, side, lambda n, w: n == -1); wo.append(o); wq.append(q)
wo, wq = np.concatenate(wo), np.concatenate(wq); s = np.argsort(wo, kind='stable')
patches.append(('reactingWall', 'wall', wo[s], wq[s], ''))

allq = np.concatenate([fq] + [p[3] for p in patches]); allo = np.concatenate([fo] + [p[2] for p in patches])
used = np.unique(allq); newid = -np.ones((nx + 1)*(ny + 1)*(nz + 1), dtype=np.int64); newid[used] = np.arange(len(used))
allq = newid[allq]
pi_ = used % (nx + 1); pj = (used // (nx + 1)) % (ny + 1); pk = used // ((nx + 1)*(ny + 1))
pts = np.stack([x0 + pi_*h, x0 + pj*h, z0 + pk*h], 1)

os.makedirs('constant/polyMesh', exist_ok=True)
def head(cls, obj, note=''):
    return ("FoamFile\n{\n    format      ascii;\n    class       %s;\n%s    location    \"constant/polyMesh\";\n    object      %s;\n}\n\n" % (cls, note, obj))
note = '    note        "nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d";\n' % (len(pts), ncell, len(allq), nint)
with open('constant/polyMesh/points', 'w') as f:
    f.write(head('vectorField', 'points') + "%d\n(\n" % len(pts)); np.savetxt(f, pts, fmt='(%.10g %.10g %.10g)'); f.write(")\n")
with open('constant/polyMesh/faces', 'w') as f:
    f.write(head('faceList', 'faces') + "%d\n(\n" % len(allq)); np.savetxt(f, allq, fmt='4(%d %d %d %d)'); f.write(")\n")
with open('constant/polyMesh/owner', 'w') as f:
    f.write(head('labelList', 'owner', note) + "%d\n(\n" % len(allo)); np.savetxt(f, allo, fmt='%d'); f.write(")\n")
with open('constant/polyMesh/neighbour', 'w') as f:
    f.write(head('labelList', 'neighbour', note) + "%d\n(\n" % nint); np.savetxt(f, fn, fmt='%d'); f.write(")\n")
with open('constant/polyMesh/boundary', 'w') as f:
    f.write(head('polyBoundaryMesh', 'boundary') + "%d\n(\n" % len(patches)); start = nint
    for name, typ, o, q, extra in patches:
        f.write("    %s\n    {\n        type            %s;\n%s%s        nFaces          %d;\n        startFace       %d;\n    }\n"
                % (name, typ, "        inGroups        List<word> 1(wall);\n" if typ == 'wall' else '', extra, len(o), start)); start += len(o)
    f.write(")\n")

vol_l = (x1 - x0)**2*(zl1 - zl0); nl = int(fluid[:, :, (zc > zl0) & (zc < zl1)].sum())
print("h = %g mm, grid %d x %d x %d, cells %d, internal faces %d" % (h_mm, nx, ny, nz, ncell, nint))
print("porosity of the lattice %.3f, wall area %.1f mm2 (stair-step), per lattice volume %.0f m2/m3"
      % (nl*h**3/vol_l, len(wo)*h*h*1e6, len(wo)*h*h/vol_l))
for name, typ, o, q, extra in patches: print("  %-13s %-7s %6d faces" % (name, typ, len(o)))
