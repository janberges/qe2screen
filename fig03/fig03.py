#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

comm = elphmod.MPI.comm

labels = 'abcdefghij'

style = dict(cut=True, line_cap='butt', line_join='miter')

Margin = 0.9
margin = 0.12

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

nbnd = [1, 5, 13, 17, 22]

orbitals = [r'Ta-$d_{0, \pm 2}$', r'$+$\,Ta-$d$', r'$+$\,S-$s{,}p$',
    r'$+$\,Ta-$s{,}p$', r'$+$\,MLWF']

mu = -5.8726

k, x, GMKG = elphmod.bravais.path('GMKG', ibrav=4)

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

        height=5.0,

        margin=margin,
        left=Margin,
        top=0.4,

        title=r'%d band' % n,

        ymax=11.0,
        ystep=15.0,
        ypadding=1.0,

        xticks=x[GMKG],
        xmarks=False,

        ylabel='Electron energy (eV)',
        )

    if i != 0:
        plot.left = margin
        plot.ylabel = None
        plot.ymarks = False

    if n != 1:
        plot.title += 's'

    plot.width = (Margin + margin * (2 * len(nbnd) - 1)
        - plot.single) / len(nbnd)

    plot.line(grid=True)

    for m in range(len(e0)):
        plot.line(x0, e0[m], color='lightgray', **style)

    for m in range(1 if len(ew) == 3 else len(ew)):
        plot.line(xw, ew[m], color=mauve, thick=True, **style)

    plot.node((x[0] + x[-1]) / 2, -52.5, orbitals[i], inner_sep='2pt',
        rounded_corners='1pt', draw='lightgray', fill='white')

    plot.node((x[GMKG[0]] + x[GMKG[1]]) / 2, -20.0, '(%s)' % labels[i])

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

        height=5.0,

        margin=margin,
        left=Margin,
        bottom=0.4,

        ymin=0.0,
        ymax=174.0,
        ystep=25.0,

        xticks=zip(x[GMKG], [
            r'$\Gamma$',
            r'$\mathrm M$',
            r'$\mathrm K$',
            r'$\Gamma$',
            ]),

        ylabel='Phonon energy (meV)',
        )

    plot.width = (Margin + margin * (2 * len(nbnd) - 1)
        - plot.single) / len(nbnd)

    plot.line(grid=True)

    if i != 0:
        plot.left = margin
        plot.ylabel = None
        plot.ymarks = False

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

    plot.node((x[GMKG[0]] + x[GMKG[1]]) / 2, 165.0,
        '(%s)' % labels[len(nbnd) + i])

    plot.save('fig03%s.pdf' % labels[len(nbnd) + i])

if comm.rank != 0:
    raise SystemExit

storylines.combine('fig03.bg.pdf', ['fig03' + a for a in labels],
    columns=len(nbnd))

lx = 1.5 * plot.width - margin + Margin
ly = plot.bottom + (plot.height - plot.bottom - plot.top) * (112.5 / plot.ymax)

plot = storylines.Plot(
    style='APS',

    width=plot.single,
    height=2 * plot.height,

    margin=0.0,
    xyaxes=False,

    background='fig03.bg.pdf',

    lpos='BL',
    lcol=2,
    llen='4mm',
    lwid=4,
    lopt='xshift=%gcm, yshift=%gcm, inner ysep=3pt, rounded corners=1pt, '
        'draw=lightgray,  fill=white' % (lx, ly),
    )

plot.line(thick=False, color=mauve, label='cDFPT', **style)
plot.line(thick=True, color=orange, label='DFPT', **style)
plot.line(thick=False, color=orange, dash_pattern='on 0.5mm off 0.2mm',
    label=r'DFPT\textminus$\varPi^{\text{00}}$', **style)
plot.line(thick=True, color=mauve, dash_pattern='on 0.8mm off 0.8mm',
    label=r'cDFPT+$\varPi^{\text{\rlap{p0}}}$', **style)

plot.save('fig03.pdf')
