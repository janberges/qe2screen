#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize
import storylines

comm = elphmod.MPI.comm
info = elphmod.MPI.info

info('Load tight-binding, mass-spring, and coupling models')

el = elphmod.el.Model('TaS2')
ph = elphmod.ph.Model('dyn', lr=False)
elph = elphmod.elph.Model('work/TaS2.epmatwp', 'work/TaS2.wigner', el, ph)
elph.sample_orig()

info('Set up q-point path')

q, x, GM = elphmod.bravais.path([(0.0, 0.5 / 6e3, 0.0), (0.0, 0.5, 0.0)],
    ibrav=4, N=198)

q0, x0, w0 = elphmod.el.read_bands('ref.freq')

q0 = 2 * np.pi * np.dot(q0, ph.a.T) / np.linalg.norm(ph.a[0])
x0 += x[-1] - x0[-1]

info('Load reference data')

sqrtM = np.sqrt(np.repeat(ph.M, 3))

D0 = np.empty((len(q0), ph.size, ph.size), dtype=complex)

for iq in range(len(q0)):
    D0[iq] = elphmod.ph.read_flfrc('ref%d' % (iq + 1))[0][1][0]

D0 /= sqrtM[np.newaxis, np.newaxis, :]
D0 /= sqrtM[np.newaxis, :, np.newaxis]

d0 = np.empty((len(q0), ph.size))

iq = 0

with open('ph0.out') as lines:
    for line in lines:
        if 'Printing the electron-phonon matrix elements' in line:
            next(lines)

            for i in range(ph.size):
                d0[iq, i] = float(next(lines))

            iq += 1

d0 /= sqrtM[np.newaxis, :]

info('Optimize long-range separation parameter')

ph.lr = True

def objective(L):
    ph.L, = L
    ph.update_short_range()

    return ph.sum_force_constants()

scipy.optimize.minimize(objective, [1.0], tol=0.1)

elph.update_short_range()

info('Interpolate dynamical matrix and coupling for Q = 0')

def sample(q):
    D = elphmod.dispersion.sample(ph.D, q)

    d = elphmod.dispersion.sample(elph.g, q, elbnd=True,
        comm=elphmod.MPI.I)[:, :, 0, 0]

    return D, d

Di, di = sample(q)

info('Optimze quadrupole tensors')

def error():
    D, d = sample(q0)

    dD = (abs(D - D0) ** 2).sum()
    dd = (abs(abs(d) - d0) ** 2).sum()

    return dD, dd

dD0, dd0 = error()

def objective(Q):
    ph.Q = np.zeros((ph.nat, 3, 3, 3))

    ph.Q[1, 1, 1, 1] = Q[0] # Ta y y y
    ph.Q[2, 1, 1, 1] = Q[1] # S  y y y
    ph.Q[2, 2, 1, 1] = Q[2] # S  z y y

    ph.Q[:, 0, 0, 1] = ph.Q[:, 0, 1, 0] = ph.Q[:, 1, 0, 0] = -ph.Q[:, 1, 1, 1]
    ph.Q[:, 2, 0, 0] = ph.Q[:, 2, 1, 1]

    ph.Q[0, :2] = ph.Q[2, :2]
    ph.Q[0, 2] = -ph.Q[2, 2]

    ph.Q[Q == 0.0] = 0.0 # avoid negative zeros

    ph.update_short_range()
    elph.update_short_range()

    dD, dd = error()

    dD /= dD0
    dd /= dd0

    info('error(D) = %.10g%%' % (100 * dD))
    info('error(d) = %.10g%%' % (100 * dd))

    return np.sqrt(dD ** 2 + dd ** 2)

scipy.optimize.minimize(objective, np.ones(3), tol=1e-3)

info('Interpolate dynamical matrix and coupling for optimal Q')

D, d = sample(q)

info('Plot results')

if comm.rank != 0:
    raise SystemExit

elphmod.ph.write_quadrupole_fmt('quadrupole.fmt', ph.Q)

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

margin = 0.12

marks = dict(mark='*', mark_size='0.8pt', only_marks=True)

plot = storylines.Plot(
    style='APS',
    preamble=r'\usepackage{bm}',

    margin=margin,
    left=1.0,
    bottom=0.4,

    grid=True,

    xticks=[(0, r'$\Gamma$'), (x[-1] / 2, None), (x[-1], 'M')],

    lpos='rt',
    lopt='below left=1mm',
    lbox=True,
    )

plot.width = plot.single / 2

plot.ylabel = r'Diagonal of dynamical matrix (eV\textsuperscript2)'
plot.ymin = 0.08
plot.ymax = 0.1225
plot.ystep = 0.01
plot.axes()

for i in range(ph.size):
    plot.line(x, Di[:, i, i].real * elphmod.misc.Ry ** 2,
        color=mauve, thick=True, label=r'$\bm Z$')

for i in range(ph.size):
    plot.line(x, D[:, i, i].real * elphmod.misc.Ry ** 2,
        color=darkmauve, thick=True, label=r'$\bm Z^*, Q$')

for i in range(ph.size):
    plot.line(x0, D0[:, i, i].real * elphmod.misc.Ry ** 2, **marks)

plot.save('fitQa.pdf')

plot.clear()

plot.ylabel = r'Electron-phonon coupling (eV\textsuperscript{3/2})'
plot.ymin = 0.0
plot.ymax = 0.65
plot.ystep = 0.1
plot.axes()

for i in range(ph.size):
    plot.line(x, abs(di[:, i]) * elphmod.misc.Ry ** 1.5,
        color=orange, thick=True, label=r'$\bm Z$')

for i in range(ph.size):
    plot.line(x, abs(d[:, i]) * elphmod.misc.Ry ** 1.5,
        color=darkorange, thick=True, label=r'$\bm Z^*, Q$')

for i in range(ph.size):
    plot.line(x0, d0[:, i] * elphmod.misc.Ry ** 1.5, **marks)

plot.save('fitQb.pdf')

storylines.combine('fitQ.pdf', ['fitQ%s' % a for a in 'ab'])
