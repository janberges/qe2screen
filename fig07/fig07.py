#!/usr/bin/env python3

import elphmod
import storylines
import numpy as np

comm = elphmod.MPI.comm

labels = 'ab'

cmap = storylines.colormap(
    (0.00, storylines.Color(255, 255, 255)),
    (0.05, storylines.Color(94, 60, 153)),
    (0.10, storylines.Color(230, 97, 1)),
    (1.00, storylines.Color(241, 163, 64)))

margin = 0.12

#q, x, A0 = elphmod.el.read_bands('phband0rw1.freq')
#q, x, Ad = elphmod.el.read_bands('phbanddrw1.freq')
#
#minimum = min(A0.min(), Ad.min())
#
#lgA0 = np.log10(A0 / minimum)
#lgAd = np.log10(Ad / minimum)
#
#orders = max(lgA0.max(), lgAd.max())
#
#if comm.rank == 0:
#    storylines.save('phband0rw1.png', lgA0[::-1, :, None] * 255 / orders)
#    storylines.save('phbanddrw1.png', lgAd[::-1, :, None] * 255 / orders)
#
#elphmod.MPI.info(orders)

orders = 8.678954153853722

q0, x0, w0 = elphmod.el.read_bands('phband0r1.freq')
qd, xd, wd = elphmod.el.read_bands('phbanddr1.freq')

A0 = 10 ** (np.array(storylines.load('phband0rw1.png'))[:, :, 0] * orders / 255)
Ad = 10 ** (np.array(storylines.load('phbanddrw1.png'))[:, :, 0] * orders / 255)

maximum = max(A0.max(), Ad.max())

q, x, GMKG = elphmod.bravais.path('GMKG', ibrav=4)

for n, (w, A, title) in enumerate([
        (w0, A0, 'No doping'),
        (wd, Ad, 'Van Hove filling')]):

    image = elphmod.plot.color(A, minimum=0.0, maximum=maximum, cmap=cmap)

    if comm.rank == 0:
        storylines.save('fig07%s.png' % labels[n], image)

        plot = storylines.Plot(
            style='APS',
            preamble=r'\usepackage{mathtools}',

            title='\smash{(%s)} %s' % (labels[n], title),

            left=0.9,
            right=0.8,
            bottom=0.4,
            top=0.4,

            xticks=zip(x[GMKG], [
                r'$\Gamma$',
                r'$\mathrm M$',
                r'$\mathrm K$',
                r'$\Gamma$',
                ]),

            ylabel='Phonon energy (meV)',
            yformat=lambda y: '$%g\,\mathrm i$' % abs(y)
                if y < 0 else '$%g$' % y,

            zlabel='Spectral function',
            zticks=[0, (1, r'\llap{max}')],
            zclose=True,

            ymin=-25.0,
            ymax=50.0,
            ystep=10.0,
            zmin=0.0,
            zmax=1.0,

            cmap=cmap,

            lpos='br',
            lopt='above left=1mm, inner ysep=3pt, rounded corners=1pt, '
                'draw=lightgray, fill=white',
            )

        plot.width = (plot.left + plot.right + 2 * margin - plot.single) / 2
        plot.height = plot.width * (plot.ymax - plot.ymin) / plot.ymax

        plot.image('fig07%s.png' % labels[n], 0, 0, x[-1], 50.0)

        plot.axes()

        for point in x[GMKG[1:-1]]:
            plot.line(x=point, color='gray')

        plot.line(y=0, color='gray')

        for nu in range(w.shape[0]):
            plot.line(x0, w[nu], cut=True, label=None if n else r'$\omega = 0$')

        if n == 0:
            plot.colorbar = False
            plot.right = margin
        else:
            plot.ylabel = None
            plot.ymarks = False
            plot.left = margin

        plot.save('fig07%s.pdf' % labels[n])

if comm.rank == 0:
    print('Height: %g cm (for adjacent Fig. 6)' % plot.height)

    storylines.combine('fig07.pdf', ['fig07%s' % a for a in labels])
