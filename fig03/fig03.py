#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

comm = elphmod.MPI.comm

height = 11.0 + 2.0 / 3.0
Margin = 0.96
margin = 0.12

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)

style = dict(cut=True, line_cap='butt', line_join='miter')

labels = 'abcdefghij'

nbnd = [1, 5, 13, 17, 22]

orbitals = [r'Ta-$d_{z^2, x^2 - y^2, x y}$', r'$+$\,Ta-$d_{x z, y z}$',
    r'$+$\,S-$s{,}p$', r'$+$\,Ta-$s{,}p$', r'$+$\,MLWF']

mu = -5.8726

k, x, corners = elphmod.bravais.path('GMKG', ibrav=4)

k0, x0, e0 = elphmod.el.read_bands('bands.dat')
e0 -= mu

#a = elphmod.bravais.primitives(ibrav=4)
#b = elphmod.bravais.reciprocals(*a)
#k = 2 * np.pi * elphmod.bravais.cartesian_to_crystal(k0, *b)

for i, n in enumerate(nbnd):
    #el = elphmod.el.Model('wannier%02d' % n)
    #e, order = elphmod.dispersion.dispersion(el.H, k, order=True)
    #elphmod.el.write_bands('bands%02d.dat' % n, k0, e.T)

    kw, xw, ew = elphmod.el.read_bands('bands%02d.dat' % n)
    ew -= mu

    if comm.rank != 0:
        continue

    plot = storylines.Plot(
        style='APS',

        left=Margin,
        right=margin,
        bottom=margin,
        top=Margin / 2,

        title=r'%d band' % n,

        ymax=12.0,
        ymin=-70.0,
        ystep=15.0,

        xticks=x[corners],
        xlabels=False,

        ylabel='Electron energy (eV)',
        )

    if i != 0:
        plot.left = margin
        plot.ylabel = None
        plot.ylabels = False

    if n != 1:
        plot.title += 's'

    plot.width = (Margin + (2 * len(nbnd) - 1) * margin
        - plot.double) / len(nbnd)
    plot.height = (Margin + 2 * margin - height) / 2

    plot.line(y=0.0, color='lightgray')

    for m in range(len(e0)):
        plot.line(x0, e0[m], color='lightgray', **style)

    for m in range(n):
        plot.line(xw, ew[m], color=mauve, thick=True, **style)

    if n == 1:
        for m in range(1, 3):
            plot.line(xw, ew[m], color=orange, **style)

    plot.node((x[0] + x[-1]) / 2, (e0[0, 0] + e0[1, 0]) / 2, orbitals[i])

    plot.node(x[0], plot.ymax, '(%s)' % labels[i], below_right='0.5mm',
        inner_sep='1pt', rounded_corners=True, fill='white')

    plot.save('fig03%s.pdf' % labels[i])

q0, x0, w0 = elphmod.el.read_bands('phband.freq')

for i, n in enumerate(nbnd):
    qc, xc, wc = elphmod.el.read_bands('phband%02dc.freq' % n)
    qu, xu, wu = elphmod.el.read_bands('phband%02dr.freq' % n)
    qr, xr, wr = elphmod.el.read_bands('phband%02dcr.freq' % n)

    if comm.rank != 0:
        continue

    plot = storylines.Plot(
        style='APS',

        left=Margin,
        right=margin,
        bottom=Margin / 2,
        top=margin,

        ymin=0.0,
        ymax=174.0,
        ystep=25.0,

        xticks=list(zip(x[corners], [
            r'$\Gamma$',
            r'$\mathrm M$',
            r'$\mathrm K$',
            r'$\Gamma$',
            ])),

        ylabel='Phonon energy (meV)',

        lpos='cmt',
        llen='4mm',
        )

    plot.width = (Margin + (2 * len(nbnd) - 1) * margin
        - plot.double) / len(nbnd)
    plot.height = (Margin + 2 * margin - height) / 2

    if i != 0:
        plot.left = margin
        plot.ylabel = None
        plot.ylabels = False

    for nu in range(len(w0)):
        plot.line(x0, w0[nu], thick=True, color=orange, **style)

    for nu in range(len(wr)):
        plot.line(xr, wr[nu], thick=True, color=mauve,
            dash_pattern='on 0.8mm off 0.8mm', **style)

    for nu in range(len(wc)):
        plot.line(xc, wc[nu], thick=False, color=mauve, **style)

    for nu in range(len(wc)):
        plot.line(xu, wu[nu], thick=False, color=orange,
            dash_pattern='on 0.5mm off 0.2mm', **style)

    if i == 0:
        plot.line(thick=True, color=orange, label='DFPT', **style)
        plot.line(thick=False, color=mauve, label='cDFPT', **style)
        plot.line(thick=True, color=mauve, dash_pattern='on 0.8mm off 0.8mm',
            label=r'cDFPT+$\varPi^{\text{p0}}$', **style)
        plot.line(thick=False, color=orange, dash_pattern='on 0.5mm off 0.2mm',
            label=r'DFPT\textminus$\varPi^{\text{00}}$', **style)

    plot.node(x[0], plot.ymax, '(%s)' % labels[len(nbnd) + i],
        below_right='0.5mm', inner_sep='1pt', rounded_corners=True,
        fill='white')

    plot.save('fig03%s.pdf' % labels[len(nbnd) + i])

if comm.rank == 0:
    storylines.combine('fig03.pdf', ['fig03' + a for a in labels],
        columns=len(nbnd))
