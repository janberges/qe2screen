#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize
import storylines

comm = elphmod.MPI.comm
info = elphmod.MPI.info

plot_coupling_in_mode_basis = True # otherwise in Cartesian basis

info('Load tight-binding, mass-spring, and coupling models')

el = elphmod.el.Model('TaS2')
ph = elphmod.ph.Model('dyn', lr=False)
elph = elphmod.elph.Model('work/TaS2.epmatwp', 'work/TaS2.wigner', el, ph)
elph.sample_orig()

info('Set up q-point path')

q, x, special = elphmod.bravais.path([(0.0, 0.5 / d, 0.0) for d in [6e3, 2, 1]],
    ibrav=4, N=198)

q0, x0, w0 = elphmod.el.read_bands('dynref.freq')

q0 = 2 * np.pi * np.dot(q0, ph.a.T) / np.linalg.norm(ph.a[0])
x0 += x[-1] - x0[-1]

info('Load reference data')

sqrtM = np.sqrt(np.repeat(ph.M, 3))

D0 = np.empty((len(q0), ph.size, ph.size), dtype=complex)

for iq in range(len(q0)):
    D0[iq] = elphmod.ph.read_flfrc('dynref%d' % (iq + 1))[0][1][0]

D0 /= sqrtM[np.newaxis, np.newaxis, :]
D0 /= sqrtM[np.newaxis, :, np.newaxis]

g0 = np.empty((len(q0), ph.size), dtype=complex)

iq = 0

with open('phref.out') as lines:
    for line in lines:
        if 'Printing the electron-phonon matrix elements' in line:
            next(lines)

            for i in range(ph.size):
                re, im = next(lines).split()
                g0[iq, i] = float(re) + 1j * float(im)

            iq += 1

g0 /= sqrtM[np.newaxis, :]

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

    g = elphmod.dispersion.sample(elph.g, q, elbnd=True,
        comm=elphmod.MPI.I)[:, :, 0, 0]

    return D, g

Dd, gd = sample(q)

info('Optimze quadrupole tensors')

def error():
    D, g = sample(q0)

    dD = (abs(D - D0) ** 2).sum()
    dg = (abs(abs(g) - abs(g0)) ** 2).sum()

    return dD, dg

dD0, dg0 = error()

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

    dD, dg = error()

    dD /= dD0
    dg /= dg0

    info('error(D) = %.10g%%' % (100 * dD))
    info('error(g) = %.10g%%' % (100 * dg))

    return np.sqrt(dD ** 2 + dg ** 2)

scipy.optimize.minimize(objective, np.ones(3), tol=1e-3)

info('Interpolate dynamical matrix and coupling for optimal Q')

Dq, gq = sample(q)

info('Plot results')

if comm.rank != 0:
    raise SystemExit

elphmod.ph.write_quadrupole_fmt('quadrupole.fmt', ph.Q)

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)

marks = dict(mark='*', mark_size='0.8pt', only_marks=True)

plot = storylines.Plot(
    style='APS',

    xticks=list(zip(x[special], [r'$\Gamma$', None, 'M'])),

    lpos='rt',
    lopt='below left',
    )

plot.width = plot.single / 2

plot.ylabel = 'Phonon energy (meV)'
plot.axes()

wd2, ud = np.linalg.eigh(Dd)
wq2, uq = np.linalg.eigh(Dq)
w02, u0 = np.linalg.eigh(D0)

wd = elphmod.ph.sgnsqrt(wd2)
wq = elphmod.ph.sgnsqrt(wq2)
w0 = elphmod.ph.sgnsqrt(w02)

for nu in range(ph.size):
    plot.line(x, wd[:, nu] * 1e3 * elphmod.misc.Ry, color=orange, thick=True)

for nu in range(ph.size):
    plot.line(x, wq[:, nu] * 1e3 * elphmod.misc.Ry, color=mauve, thick=True)

for nu in range(ph.size):
    plot.line(x0, w0[:, nu] * 1e3 * elphmod.misc.Ry, **marks)

plot.save('fitQa.pdf')

plot.clear()

plot.ylabel = r'Electron-phonon coupling (eV\textsuperscript{3/2})'
plot.ymin = 0.0
plot.axes()

if plot_coupling_in_mode_basis:
    gd = np.einsum('qx,qxv->qv', gd, ud)
    gq = np.einsum('qx,qxv->qv', gq, uq)
    g0 = np.einsum('qx,qxv->qv', g0, u0)

gd = np.sort(abs(gd), axis=1)
gq = np.sort(abs(gq), axis=1)
g0 = np.sort(abs(g0), axis=1)

for nu in range(ph.size):
    plot.line(x, gd[:, nu] * elphmod.misc.Ry ** 1.5, color=orange, thick=True,
        label=r'$\mathbf Z^*$ only')

for nu in range(ph.size):
    plot.line(x, gq[:, nu] * elphmod.misc.Ry ** 1.5, color=mauve, thick=True,
        label=r'$\mathbf Z^*$ and $Q$')

for nu in range(ph.size):
    plot.line(x0, g0[:, nu] * elphmod.misc.Ry ** 1.5, **marks)

plot.save('fitQb.pdf')

storylines.combine('fitQ.pdf', ['fitQ%s' % a for a in 'ab'])
