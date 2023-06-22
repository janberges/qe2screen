#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize
import storylines

comm = elphmod.MPI.comm
info = elphmod.MPI.info

Margin = 0.96
margin = 0.12

fast = True # reduce accuracy to speed up calculations?
correct = False # subtract errors from approximations?

white = storylines.Color(255, 255, 255)
orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

dashed = 'on 0.8mm off 0.8mm'
dotted = 'on 0mm off 0.2mm on 0.4mm off 0.2mm'

t = 1.0
wb = 0.05
g = 0.02
U0 = 1.0
V0 = 0.5

occ = elphmod.occupations.fermi_dirac
eta = 0.0004
nk = 55440

if fast:
    eta *= 10
    nk //= 10

kT = np.array([1.0, 0.1, 0.01])
k = np.linspace(0, 2 * np.pi, nk, endpoint=False)
q = k[::42]

qmin = 1 / 3 - 0.05
qmax = 1 / 2 + 0.05

e = -t * np.cos(k)
wb2 = wb ** 2 * (1 - np.cos(q))
g0 = 1j * g * (np.sin(k) - np.sin(np.add.outer(q, k)))

w, dw = np.linspace(0.0, 0.083, 726, retstep=True)

Amax = 120.0

W = np.empty((len(kT), len(q)))
C = np.empty((len(kT), len(q)))

wb2_pi = 2 * wb ** 2
g2_pi = (2 * g * np.sin(k)) ** 2

def W2_pi(kT):
    return wb2_pi + elphmod.diagrams.phonon_self_energy([[np.pi]],
        e.reshape((-1, 1)), g2=g2_pi, kT=kT)[0, 0]

solution = scipy.optimize.root_scalar(W2_pi, bracket=[0.001, 0.01])

info('T(CDW) = %.2f eV' % (1e3 * solution.root))

# ω₀²t/g₀² < ∫[0, 2π] dk tanh(t/2T cos(k)) tan(k) sin(k)

for i in range(len(kT)):
    W2 = wb2 + elphmod.diagrams.phonon_self_energy(q.reshape((-1, 1)),
        e.reshape((-1, 1)), g2=abs(g0) ** 2, kT=kT[i])[:, 0]

    W[i] = elphmod.ph.sgnsqrt(W2)
    C[i] = np.fft.ifft(W2).real

C *= 1e6

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        margin=margin,
        left=Margin,
        bottom=Margin,

        xlabel=r'Distance (sites)',
        ylabel=r'Force constant (meV\textsuperscript2)',

        xmax=33.0,
        ymin=-abs(C[:, 1:]).max(),
        ymax=0.2 * abs(C[:, 1:]).max(),

        xstep=10.0,
        ystep=500.0,

        lower=orange,
        upper=mauve,

        colorbar=False,
        )

    plot.height = 0.42 * (Margin + margin - plot.single)
    plot.width = 0.42 * (Margin + 7 * margin - plot.double) / 2 * (qmax
        / (2 * qmax - qmin))

    plot.axes()

    d = np.arange(len(q))

    style = dict(thick=True, dotted=True, cut=True)

    for i in range(len(kT)):
        plot.line(d, -abs(C[i]), i, **style)
        plot.line(d, abs(C[i]), color='lightgray', zindex=1, **style)

    for i in range(len(kT)):
        plot.line(d, C[i], i, mark='*', mark_size='1pt', only_marks=True,
            cut=True)

    plot.save('fig02a.in.pdf')

    plot = storylines.Plot(
        style='APS',

        preamble=r'\usepackage{nicefrac}',

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        ylabel='Phonon energy (meV)',

        xmax=1.1 * np.pi,
        ymin=1e3 * w[0],
        ymax=1e3 * w[-1],
        ystep=10.0,

        xticks=[0, (np.pi / 3, r'$\nicefrac16$'), (2 * np.pi / 3,
            r'$\nicefrac13$'), (np.pi, r'$\nicefrac12$')],

        colorbar=False,

        lower=orange,
        upper=mauve,

        ltop=r'$T$ (meV)',
        lali='left',
        lpos='lt',
        lopt='below right=4mm',
        )

    plot.height = Margin + margin - plot.single
    plot.width = (Margin + 7 * margin - plot.double) / 2 * (qmax
        / (2 * qmax - qmin))

    plot.axes()

    for i in range(len(kT)):
        plot.line(q, 1e3 * W[i], i, cut=True, thick=True,
            label='%g' % (1e3 * kT[i]))

    plot.node(0.0, plot.ymax, '(a)', below_right=True)

    plot.node(0.0, plot.ymin, r'\includegraphics{fig02a.in.pdf}',
        above_right='3mm', rounded_corners='1pt', draw='lightgray',
        fill='white')

    plot.save('fig02a.pdf')

kT = 0.025

kF = np.array([np.pi / 2, np.pi / 3])
qK = 2 * kF

labels = ['\pi', r'\frac23 \pi']
colors = [(darkorange, orange), (darkmauve, mauve)]

iq1 = int(round(qmin * len(q)))
iq2 = int(round(qmax * len(q)))

q = q[iq1:iq2]

sizes, bounds = elphmod.MPI.distribute(len(q), bounds=True, comm=comm)

my_wb0 = np.empty(sizes[comm.rank])
my_Abw = np.empty((sizes[comm.rank], len(w)))
my_A00 = np.empty((sizes[comm.rank], len(w)))
my_Ab0 = np.empty((sizes[comm.rank], len(w)))

wb0 = np.empty((len(kF), len(q)))
Abw = np.empty((len(kF), len(q), len(w)))
A00 = np.empty((len(kF), len(q), len(w)))
Ab0 = np.empty((len(kF), len(q), len(w)))

for ikF in range(len(kF)):
    for my_iq, iq in enumerate(range(*bounds[comm.rank:comm.rank + 2])):
        scale = 2 * np.pi / len(k)
        Q = int(round(q[iq] / scale))
        qr = Q * scale

        e = -t * np.cos(k)
        wb2 = wb ** 2 * (1 - np.cos(qr))
        g0 = 1j * g * (np.sin(k) - np.sin(k + qr))
        U = U0 - V0 * np.log(np.sin(abs(qr) / 2))

        e -= -t * np.cos(kF[ikF])

        f = occ(e / kT)
        d = occ.delta(e / kT) / (-kT)

        de = e - np.roll(e, -Q)
        df = f - np.roll(f, -Q)

        pre = 2 / len(k)

        ok = abs(de) > 1e-10

        chi0 = np.empty(len(k))
        chi0[ok] = df[ok] / de[ok]
        chi0[~ok] = d[~ok]
        chi0 *= pre

        chiw = pre * df / (de + w[:, None] + 1j * eta)

        W = U / (1 - U * chiw.sum(axis=1))

        gb = g0 - U * (chi0 * g0).sum()
        gw = gb + (W * (chiw * gb).sum(axis=1))[:, None]

        Pib00 = (gb.conj() * chi0 * g0).sum()
        Pi000 = (g0.conj() * chi0 * g0).sum()
        Pibww = (gb.conj() * chiw * gw).sum(axis=1)
        Pib0w = (gb.conj() * chiw * g0).sum(axis=1)
        Pi00w = (g0.conj() * chiw * g0).sum(axis=1)

        # To swap k and k + q:
        # print(np.allclose(g0.conj(), -np.roll(g0[::-1], 1 - Q)))
        # The first holds for g0 only, the second for both g0 and gw

        w2 = wb2 + Pib00

        Dbw = w2 - Pib00 + Pibww
        Db0 = w2 - Pib00 + Pib0w
        D00 = w2 - Pi000 + Pi00w

        if correct:
            Db0 -= (gb.conj() * chiw * (g0 - gw)).sum(axis=1)
            D00 -= np.average(g0 - gw, axis=1) ** 2 / W

        my_wb0[my_iq] = elphmod.ph.sgnsqrt(w2.real)
        my_Abw[my_iq] = -1 / np.pi * (2 * w / (w ** 2 - Dbw)).imag
        my_Ab0[my_iq] = -1 / np.pi * (2 * w / (w ** 2 - Db0)).imag
        my_A00[my_iq] = -1 / np.pi * (2 * w / (w ** 2 - D00)).imag

    comm.Allgatherv(my_wb0, (wb0[ikF], comm.allgather(my_wb0.size)))
    comm.Allgatherv(my_Abw, (Abw[ikF], comm.allgather(my_Abw.size)))
    comm.Allgatherv(my_A00, (A00[ikF], comm.allgather(my_A00.size)))
    comm.Allgatherv(my_Ab0, (Ab0[ikF], comm.allgather(my_Ab0.size)))

A = Abw.copy()
A[A > Amax] = Amax
A /= Amax

image = -255

for ikF, color in enumerate([orange, mauve]):
    image += elphmod.plot.color(A[ikF].T[::-1],
        storylines.colormap((0.0, white), (1.0, color)), 0.0, 1.0)

if comm.rank == 0:
    storylines.save('fig02b.png', image)

    plot = storylines.Plot(
        style='APS',

        preamble=r'\usepackage{nicefrac}',

        margin=margin,
        bottom=Margin / 2,

        xticks=[0, (np.pi / 3, r'$\nicefrac16$'), (2 * np.pi / 3,
            r'$\nicefrac13$'), (np.pi, r'$\nicefrac12$')],

        ymin=1e3 * w[0],
        ymax=1e3 * w[-1],
        ylabels=False,

        background='fig02b.png',

        ltop=r'$T = 25$~meV',
        lali='left',
        lpos='lt',
        lopt='below right=4mm',
        )

    plot.height = Margin + margin - plot.single
    plot.width = (Margin + 7 * margin - plot.double) / 2 * ((qmax - qmin)
        / (2 * qmax - qmin))

    for ikF in range(len(kF)):
        plot.line(q, 1e3 * wb0[ikF], label=r'$\omega = 0$')

    plot.node(q[0], plot.ymax, '(b)', below_right=True)

    plot.save('fig02b.pdf')

    plot = storylines.Plot(
        style='APS',
        preamble=r'\usepackage{mathtools}',

        margin=margin,
        bottom=Margin / 2,

        dright=2 * margin,

        ymin=1e3 * w[0],
        ymax=1e3 * w[-1],
        ylabels=False,

        xmin=0.0,
        xmax=1e-3 * Amax,
        xstep=0.03,

        lpos='rt',
        lopt='below left',
        llen='4mm',

        line_cap='butt',
        line_join='miter',
        )

    plot.height = Margin + margin - plot.single
    plot.width = (Margin + 7 * margin - plot.double) / 4

    plot.axes()

    for ikF in range(len(kF)):
        iqK = np.argmin(abs(q - qK[ikF]))

        print('q = %s' % labels[ikF])
        print('Integral(bw): %g' % (Abw[ikF, iqK].sum() * dw))
        print('Integral(00): %g' % (A00[ikF, iqK].sum() * dw))
        print('Integral(b0): %g' % (Ab0[ikF, iqK].sum() * dw))

        plot.line(1e-3 * Abw[ikF, iqK, ::-1], w[::-1] * 1e3,
            thick=True, color=colors[ikF][0])

        plot.line(1e-3 * A00[ikF, iqK, ::-1], w[::-1] * 1e3,
            thick=True, color=colors[ikF][1], dash_pattern=dashed)

        plot.line(1e-3 * Ab0[ikF, iqK, ::-1], w[::-1] * 1e3,
            thick=True, color=colors[ikF][1], dash_pattern=dotted)

        iw = np.argmax(Abw[ikF, iqK])
        iw = np.argmin(abs(Abw[ikF, iqK, iw:] - Abw[ikF, iqK, iw] / 2)) + iw

        plot.node(1e-3 * Abw[ikF, iqK, iw], w[iw] * 1e3,
            r'$q = 2 k_{\text F} = %s$' % labels[ikF], color=colors[ikF][0],
            above_right=True)

    plot.line(thick=True, color='darkgray',
        label=r'$g\mathrlap{^{\text b}}'
        r'\phantom{(0)} \chi^{\text b}(\omega) g(\omega)$')

    plot.line(thick=True, color='gray', dash_pattern=dashed,
        label=r'$g(0) \chi^{\text b}(\omega) g(0)$')

    plot.line(thick=True, color='gray', dash_pattern=dotted,
        label=r'$g\mathrlap{^{\text b}}'
        r'\phantom{(0)} \chi^{\text b}(\omega) g(0)$')

    plot.node(0.0, plot.ymax, '(c)', below_right=True, xshift='1mm')

    plot.save('fig02c.png')

Ry2meV = 1e3 * elphmod.misc.Ry

a = 10.0 / elphmod.misc.a0
meff = 2 / 2
kF = 0.5 * np.pi / a
q = 0.5 * kF
wb = 25.0 / Ry2meV
gb0 = 55.0 / (Ry2meV ** 1.5 / elphmod.misc.a0)
U0 = 5.0 / (Ry2meV / elphmod.misc.a0 ** 2)

occ = elphmod.occupations.fermi_dirac
kT = 0.003 / elphmod.misc.Ry
eta = 0.003 / elphmod.misc.Ry
nk = 240

if fast:
    eta *= 2
    nk //= 2

w /= elphmod.misc.Ry
dw /= elphmod.misc.Ry

k = np.linspace(0, 2 * np.pi, nk, endpoint=False) / a
k[k > np.pi / a] -= 2 * np.pi / a

k1, k2, k3 = np.meshgrid(k, k, k)

e = (k1 ** 2 + k2 ** 2 + k3 ** 2) / (2 * meff)
e -= kF ** 2 / (2 * meff)

f = occ(e / kT)
d = occ.delta(e / kT) / (-kT)

pre = 2 / nk ** 3
scale = 2 * np.pi / a / nk

for first in True, False:
    if first:
        Q = int(round(q / scale))
        q = k[Q]

        de = e - np.roll(e, shift=-Q, axis=0)
        df = f - np.roll(f, shift=-Q, axis=0)

        ok = abs(de) > 1e-10

        chi0 = pre * ((df[ok] / de[ok]).sum() + d[~ok].sum())

        sizes, bounds = elphmod.MPI.distribute(len(w), bounds=True, comm=comm)

        my_chiw = np.empty(sizes[comm.rank], dtype=complex)

        chiw = np.empty(len(w), dtype=complex)

        for my_iw, iw in enumerate(range(*bounds[comm.rank:comm.rank + 2])):
            my_chiw[my_iw] = pre * (df / (de + w[iw] + 1j * eta)).sum()

        comm.Allgatherv(my_chiw, (chiw, sizes))
    else:
        q = k[1:nk // 2]

        sizes, bounds = elphmod.MPI.distribute(len(q), bounds=True, comm=comm)

        my_chi0 = np.empty(sizes[comm.rank])
        my_chiw = np.empty(sizes[comm.rank], dtype=complex)

        chi0 = np.empty(len(q))
        chiw = np.empty(len(q), dtype=complex)

        for my_iq, iq in enumerate(range(*bounds[comm.rank:comm.rank + 2])):
            Q = int(round(q[iq] / scale))

            de = e - np.roll(e, shift=-Q, axis=0)
            df = f - np.roll(f, shift=-Q, axis=0)

            ok = abs(de) > 1e-10

            my_chi0[my_iq] = pre * ((df[ok] / de[ok]).sum() + d[~ok].sum())
            my_chiw[my_iq] = pre * (df / (de + wb + 1j * eta)).sum()

        comm.Allgatherv(my_chi0, (chi0, sizes))
        comm.Allgatherv(my_chiw, (chiw, sizes))

    gb = 1j * gb0 / q
    U = U0 / q ** 2

    g0 = gb / (1 - U * chi0)
    gw = gb / (1 - U * chiw)

    Pib00 = gb.conjugate() * chi0 * g0
    Pi000 = g0.conjugate() * chi0 * g0
    Pibww = gb.conjugate() * chiw * gw
    Pib0w = gb.conjugate() * chiw * g0
    Pi00w = g0.conjugate() * chiw * g0

    if correct:
        W = U / (1 - U * chiw)

        Pib0w -= gb.conjugate() * chiw * (g0 - gw)
        Pi00w -= (g0 - gw) ** 2 / W

    w2 = wb ** 2 + Pib00

    Dbw = w2 - Pib00 + Pibww
    Db0 = w2 - Pib00 + Pib0w
    D00 = w2 - Pi000 + Pi00w

    if first:
        Abw = -1 / np.pi * (2 * w / (w ** 2 - Dbw)).imag
        Ab0 = -1 / np.pi * (2 * w / (w ** 2 - Db0)).imag
        A00 = -1 / np.pi * (2 * w / (w ** 2 - D00)).imag

        info('Integral(bw): %g' % (Abw.sum() * dw))
        info('Integral(b0): %g' % (Ab0.sum() * dw))
        info('Integral(00): %g' % (A00.sum() * dw))
    else:
        gammabw = -1 / wb * Pibww.imag
        gammab0 = -1 / wb * Pib0w.imag
        gamma00 = -1 / wb * Pi00w.imag

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        width=2.5,

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        xmin=0.0,
        xmax=np.pi / a,
        ymin=0.0,
        ymax=20.0,
        ystep=6.0,

        xticks=[(0.0, r'$\Gamma$'), (np.pi / a, 'X')],

        ylabel='Linewidth (meV)',

        lpos='rt',
        lopt='below left',

        line_cap='butt',
        line_join='miter',
        )

    plot.height = 0.42 * (Margin + margin - plot.single)
    plot.width = 0.42 * (Margin + 7 * margin - plot.double) / 4

    plot.line(x=kF / 2, color='gray', dash_pattern=dashed)

    plot.axes()

    plot.line(q[::-1], gammabw[::-1] * Ry2meV,
        thick=True, color='darkgray')

    plot.line(q[::-1], gamma00[::-1] * Ry2meV,
        thick=True, color='gray', dash_pattern=dashed)

    plot.line(q[::-1], gammab0[::-1] * Ry2meV,
        thick=True, color='gray', dash_pattern=dotted)

    plot.save('fig02d.in.pdf')

    plot = storylines.Plot(
        style='APS',

        margin=margin,
        bottom=Margin / 2,

        dright=2 * margin,

        ymin=w[0] * Ry2meV,
        ymax=w[-1] * Ry2meV,
        ylabels=False,

        xmin=0.0,
        xmax=2e-3 * Amax,
        xstep=0.06,

        line_cap='butt',
        line_join='miter',
        )

    plot.height = Margin + margin - plot.single
    plot.width = (Margin + 7 * margin - plot.double) / 4

    plot.axes()

    plot.line(Abw[::-1] / Ry2meV, w[::-1] * Ry2meV,
        thick=True, color='darkgray')

    plot.line(Ab0[::-1] / Ry2meV, w[::-1] * Ry2meV,
        thick=True, color='gray', dash_pattern=dotted)

    plot.line(A00[::-1] / Ry2meV, w[::-1] * Ry2meV,
        thick=True, color='gray', dash_pattern=dashed)

    iw = np.argmax(Abw)
    iw = np.argmin(abs(Abw[iw:] - Abw[iw] / 2)) + iw

    plot.node(Abw[iw] / Ry2meV, w[iw] * Ry2meV,
        r'$q = \frac12 k_{\text F} = \frac12 \mathrm X$', above_right=True)

    plot.node(0.0, plot.ymax, '(d)', below_right=True)

    plot.node(plot.xmax, plot.ymax, r'\includegraphics{fig02d.in.pdf}',
        below_left='3mm', xshift='-%gcm' % plot.dright, rounded_corners='1pt',
        draw='lightgray', fill='white')

    plot.save('fig02d.png')

    for label, xlabel, left in [
        ('e', 'Phonon momentum ($2 \pi$)', Margin),
        ('f', 'Spectral weight (1/meV)', margin),
        ]:

        plot = storylines.Plot(
            style='APS',

            height=Margin / 2,

            margin=0.0,
            left=left,
            right=margin,

            xyaxes=False,
            )

        plot.width = (Margin + 3 * margin - plot.double) / 2

        plot.code(r'\node [below=\baselineskip] at (<x=0>, %g) {%s};'
            % (Margin - plot.tick, xlabel))

        plot.save('fig02%s.pdf' % label)

    storylines.combine('fig02.pdf', ['fig02%s' % a for a in 'abcdef'],
        columns=4)
