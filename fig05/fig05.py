#!/usr/bin/env python3

import elphmod
import storylines

comm = elphmod.MPI.comm

height = 11.0 + 2.0 / 3.0
Margin = 0.96
margin = 0.12

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)

labels = ['abcde', 'fghij', 'klmno', 'pqrst', 'uvwxy']

nbnd = [1, 5, 13, 17, 22]

q, x, corners = elphmod.bravais.path('GMKG', ibrav=4)

q0, x0, w0 = elphmod.el.read_bands('ref.freq')
qi, xi, wi = elphmod.el.read_bands('phband.freq')

w0 *= 1e3 * elphmod.misc.cmm1

for l, label in enumerate(['r', 'cr', 'tcr', 'bcr', 'btcr']):
    for i, n in enumerate(nbnd):
        gr, xr, wr = elphmod.el.read_bands('phband%02d%s1.freq' % (n, label))

        if comm.rank != 0:
            continue

        plot = storylines.Plot(
            style='APS',

            margin=Margin / 2,

            dbottom=0 if l == 4 else margin,
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

        plot.ymin = plot.ymax - plot.ystep * (4 if l == 3 else 3)

        plot.width = (1.5 * Margin + 2 * (len(nbnd) - 1) * margin
            - plot.double) / len(nbnd)

        plot.height = (4 if l == 3 else 3) * (height - Margin + margin) / 16

        if l == 4:
            plot.height += Margin / 2
        else:
            plot.bottom = margin
            plot.xlabels = False

        if l == 0:
            plot.height += Margin / 2 - margin

            plot.title = r'%d band' % nbnd[i]

            if n != 1:
                plot.title += 's'
        else:
            plot.top = margin

        if i != 0:
            plot.left = margin
            plot.ylabel = None
            plot.ylabels = False

        plot.line(y=0.0, color='lightgray')

        if i == len(nbnd) - 1:
            plot.node(x[-1], (plot.ymax + plot.ymin) / 2,
                r'Via $\varPi^{\text{%s}}$'
                % ['00', 'p0', 'p$T$', 'b0', 'b$T$'][l], below=True, rotate=90)
        else:
            plot.right = margin

        for nu in range(len(wi)):
            plot.line(xi[1:-1], wi[nu, 1:-1], color='lightgray', cut=True)

        plot.axes()

        for nu in range(len(wr)):
            plot.line(xr[1:-1], wr[nu, 1:-1], thick=True, cut=True,
                color=[orange, mauve, mauve, 'gray', 'gray'][l])

        for nu in range(len(w0)):
            plot.line(x0, w0[nu], mark='*', mark_size='0.3pt', only_marks=True,
                cut=True)

        plot.node(x[0], plot.ymax, '(%s)' % labels[l][i], below_right='0.5mm',
            yshift='%gcm' % -margin, inner_sep='1pt', rounded_corners=True,
            fill='white')

        plot.save('fig05%s.pdf' % labels[l][i])

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        width=Margin / 2,
        height=5.0,

        margin=0.0,
        xyaxes=False,
        )

    plot.code(r'\node [rotate=90, above=\baselineskip] at (%g, <y=0>) '
        r'{Phonon energy (meV)};' % (Margin - plot.tick))

    plot.save('fig05.l.pdf')

    storylines.combine('fig05.r.pdf',
        ['fig05%s' % a for abc in labels for a in abc], columns=len(nbnd))

    storylines.combine('fig05.pdf', ['fig05.l', 'fig05.r'], align=0.5)
