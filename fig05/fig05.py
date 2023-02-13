#!/usr/bin/env python3

import elphmod
import storylines

comm = elphmod.MPI.comm

labels = ['abcde', 'fghij', 'klmno', 'pqrst', 'uvwxy']

Margin = 0.4
margin = 0.12

nbnd = [1, 5, 13, 17, 22]

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)

q, x, GMKG = elphmod.bravais.path('GMKG', ibrav=4)

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

            margin=Margin,
            left=Margin + 0.1,
            right=Margin + 0.1,

            dbottom=0 if l == 4 else margin,
            dtop=margin,

            ymax=22.5,
            ystep=15.0,

            xticks=zip(x[GMKG], [
                r'$\Gamma$',
                r'$\mathrm M$',
                r'$\mathrm K$',
                r'$\Gamma$',
                ]),

            yformat=lambda y: '$%g\,\mathrm i$' % abs(y)
                if y < 0 else '$%g$' % y,
            )

        plot.ymin = plot.ymax - plot.ystep * (4 if l == 3 else 3)

        plot.width = (3 * Margin + 0.2 + 2 * margin * (len(nbnd) - 1)
            - plot.single) / len(nbnd)

        plot.height = ((margin if l == 4 else 2 * margin)
            + (4 if l == 3 else 3) * 0.47 * plot.width)

        plot.line(grid=True)

        if l != 4:
            plot.bottom = margin
            plot.xmarks = False

        if l == 0:
            plot.title = r'%d band' % nbnd[i]

            if n != 1:
                plot.title += 's'
        else:
            plot.top = margin

        if i != 0:
            plot.left = margin
            plot.ylabel = None
            plot.ymarks = False

        if i == len(nbnd) - 1:
            plot.node(x[-1], (plot.ymax + plot.ymin) / 2,
                r'via $\varPi^{\text{%s}}$'
                % ['00', 'p0', 'p$T$', 'b0', 'b$T$'][l], below=True, rotate=90)
        else:
            plot.right = margin

        for nu in range(len(wi)):
            plot.line(xi[1:-1], wi[nu, 1:-1], color='lightgray', cut=True)

        for nu in range(len(wr)):
            plot.fatband(xr[1:-1], wr[nu, 1:-1], thick=True,
                color=[orange, mauve, mauve, 'gray', 'gray'][l],
                thickness=0.8 * storylines.pt, protrusion=1.0, cut=True)

        plot.axes()

        for nu in range(len(w0)):
            plot.line(x0, w0[nu], mark='*', mark_size='0.3pt', only_marks=True,
                cut=True)

        plot.node(x[GMKG[2]], -30 if l == 3 else -15, '(%s)' % labels[l][i],
            above_right='-1mm')

        plot.save('fig05%s.pdf' % labels[l][i])

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        width=Margin,
        height=2 * plot.height,

        margin=0.0,
        xyaxes=False,
        )

    plot.node(0.0, 0.0, 'Phonon energy (meV)', rotate=90)
    plot.save('fig05.l.pdf')

    storylines.combine('fig05.r.pdf',
        ['fig05%s' % a for abc in labels for a in abc], columns=len(nbnd))

    storylines.combine('fig05.pdf', ['fig05.l', 'fig05.r'], align=0.5)
