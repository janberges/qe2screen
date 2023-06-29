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

labels = ['abcde', 'fghij']

nbnd = [1, 5, 13, 17, 22]

q, x, corners = elphmod.bravais.path('GMKG', ibrav=4)

q0, x0, w0 = elphmod.el.read_bands('ref.freq')
qi, xi, wi = elphmod.el.read_bands('phband.freq')

w0 *= 1e3 * elphmod.misc.cmm1

for l, label in enumerate(['r', 'cr']):
    for i, n in enumerate(nbnd):
        qr, xr, wr = elphmod.el.read_bands('phband%02d%s1.freq' % (n, label))

        if comm.rank != 0:
            continue

        plot = storylines.Plot(
            style='APS',

            margin=Margin / 2,

            dbottom=margin,
            dtop=margin,

            ymax=22.5,
            ystep=15.0,

            xticks=list(zip(x[corners], [
                r'$\Gamma$',
                r'$\mathrm M$',
                r'$\mathrm K$',
                r'$\Gamma$',
                ])),

            yformat=lambda y: '$%g\,\mathrm i$' % abs(y)
                if y < 0 else '$%g$' % y,
            )

        plot.ymin = plot.ymax - plot.ystep * 3

        plot.width = (1.5 * Margin + 2 * (len(nbnd) - 1) * margin
            - plot.double) / len(nbnd)

        plot.height = 3 * (height - Margin + margin) / 16 + Margin / 2 - margin

        if l == 0:
            plot.bottom = margin
            plot.xlabels = False

            plot.title = r'%d band' % nbnd[i]

            if n != 1:
                plot.title += 's'
        else:
            plot.top = margin

        if i != 0:
            plot.left = margin
            plot.ylabel = None
            plot.ylabels = False

        plot.line(y=0, color='lightgray')

        if i == len(nbnd) - 1:
            plot.node(x[-1], (plot.ymax + plot.ymin) / 2,
                r'via $\varPi^{\text{%s}}$' % ['00', 'p0'][l],
                below=True, rotate=90)
        else:
            plot.right = margin

        for nu in range(len(wi)):
            plot.line(xi[1:-1], wi[nu, 1:-1], color='lightgray', cut=True)

        for nu in range(len(wr)):
            plot.fatband(xr[1:-1], wr[nu, 1:-1], thick=True, cut=True,
                color=[orange, mauve][l], thickness=0.8 * storylines.pt,
                protrusion=1.0)

        plot.axes()

        for nu in range(len(w0)):
            plot.line(x0, w0[nu], cut=True,
                mark='*', mark_size='0.3pt', only_marks=True)

        if l == 1 and i > 2:
            phi = 135.0 if i == 3 else 157.5

            iq = -20 + np.argmax(wr[2, -20:])
            plot.code(r'\draw [<-, thick, shorten <=1pt] '
                '(<x=%g>, <y=%g>) -- +(%g:0.4);' % (xr[iq], wr[2, iq], phi))

        plot.node(x[0], plot.ymax, '(%s)' % labels[l][i], below_right='0.5mm',
            yshift='%gcm' % -margin, inner_sep='1pt', rounded_corners=True,
            fill='white')

        plot.save('fig11%s.pdf' % labels[l][i])

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        width=Margin / 2,
        height=3.0,

        margin=0.0,
        xyaxes=False,
        )

    plot.code(r'\node [rotate=90, above=\baselineskip] at (%g, <y=0>) '
        r'{Phonon energy (meV)};' % (Margin - plot.tick))

    plot.save('fig11.l.pdf')

    storylines.combine('fig11.r.pdf',
        ['fig11%s' % a for abc in labels for a in abc], columns=len(nbnd))

    storylines.combine('fig11.pdf', ['fig11.l', 'fig11.r'], align=0.5)
