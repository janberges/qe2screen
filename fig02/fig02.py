#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize
import storylines

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)

t = 1.0
w0 = 0.05
g0 = 0.02

ldkT = np.arange(-3, 1)
kT = 10.0 ** ldkT
q = k = np.linspace(0, 2 * np.pi, 2000, endpoint=False)

e = -t * np.cos(k)
w2 = w0 ** 2 * (1 - np.cos(q))
g = 1j * g0 * (np.sin(k) - np.sin(np.add.outer(q, k)))

W = np.empty((len(kT), len(q)))
C = np.empty((len(kT), len(q)))

w2_pi = 2 * w0 ** 2
g2_pi = (2 * g0 * np.sin(k)) ** 2

def W2_pi(kT):
    return w2_pi + elphmod.diagrams.phonon_self_energy([[np.pi]],
        e.reshape((-1, 1)), g2=g2_pi, kT=kT)[0, 0]

solution = scipy.optimize.root_scalar(W2_pi, bracket=[kT[0], kT[1]])

elphmod.MPI.info(solution)
elphmod.MPI.info('T(CDW) = %.1f K' % (solution.root / elphmod.misc.kB))

# ω₀²t/g₀² < ∫[0, 2π] dk tanh(t/2T cos(k)) tan(k) sin(k)

for i in range(len(kT)):
    W2 = w2 + elphmod.diagrams.phonon_self_energy(q.reshape((-1, 1)),
        e.reshape((-1, 1)), g2=abs(g) ** 2, kT=kT[i])[:, 0]

    W[i] = elphmod.ph.sgnsqrt(W2)
    C[i] = np.fft.ifft(W2).real

if elphmod.MPI.comm.rank != 0:
    raise SystemExit

C *= 1e6

plot = storylines.Plot(
    style='APS',

    width=4.7,
    height=3.7,

    xlabel='Momentum',
    ylabel='Phonon energy (meV)',

    xclose=True,

    xmax=1.1 * np.pi,
    ymin=0.0,
    ymax=70.0,

    xticks=[0, (np.pi, r'$\pi$')],
    ystep=20.0,

    colorbar=False,

    upper=orange,
    lower=mauve,
    )

plot.axes()

for i in range(len(kT)):
    plot.line(q, 1e3 * W[i], i, cut=True, thick=True)

plot.save('fig02.in.pdf')

plot = storylines.Plot(
    style='APS',

    height=6.0,

    xlabel=r'Distance (sites)',
    ylabel=r'Force constant (meV\textsuperscript2)',
    zlabel=r'Electronic temperature (eV)',

    xmax=30.0,
    ymin=-abs(C[:, 1:]).max(),
    ymax=0.1 * abs(C[:, 1:]).max(),

    xstep=10.0,
    ystep=250.0,
    zticks=zip(ldkT, ['%g' % z for z in kT]),

    upper=orange,
    lower=mauve,
    )

plot.axes()

d = np.arange(len(q))

style = dict(thick=True, dotted=True, cut=True)

for i in range(len(kT)):
    plot.line(d, -abs(C[i]), ldkT[i], **style)
    plot.line(d, abs(C[i]), color='lightgray', zindex=0, **style)

for i in range(len(kT)):
    plot.line(d, C[i], ldkT[i], mark='*', mark_size='1pt', only_marks=True,
        cut=True)

plot.node(plot.xmax, plot.ymin, r'\includegraphics{fig02.in.pdf}',
    above_left='1mm')

plot.save('fig02.pdf')
