#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

if elphmod.MPI.comm.rank != 0:
    raise SystemExit

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

k = np.linspace(0, np.pi, 1000, endpoint=False)

plot = storylines.Plot(
    style='APS',

    left=0.5,
    right=2.5,
    bottom=0.5,

    ymin=-7.7,
    ymax=6.3,

    xticks=[(0, r'$\Gamma$')],
    yticks=[(0, 'Fermi level')],
    )

plot.height = plot.width

plot.line(grid=True)

plot.node(k[-1], plot.ymin, 'Crystal momentum', below_left=True)
plot.node(k[0], plot.ymax, 'Electron energy', above_left=True, rotate=90)

def brace(ymin, ymax, content):
    plot.node(k[-1], (ymax + ymin) / 2, r'$\underbrace{\hspace*{<dy=%g>cm}}$'
        % (ymax - ymin), yshift='0.5mm', below=True, rotate=90)
    plot.node(k[-1], (ymax + ymin) / 2, content, align='left', anchor='west',
        xshift='4mm')

brace(2.0, 6.2, r'High-energy\\empty states')
brace(-1.4, 0.6, r'Low-energy\\active states')
brace(-3.2, -1.9, r'High-energy\\occupied states')
brace(-7.7, -7.0, r'High-energy\\core states')

style = dict(draw='none', thickness=0.1)

plot.fatband(k, 4.8 - 2.5 * np.cos(k), fill=mauve, cut=True, **style)
plot.fatband(k, 4.8 + k ** 2, fill=mauve, cut=True, **style)
plot.fatband(k, 5.7 + (np.pi - k) ** 2, fill=mauve, cut=True, **style)
plot.fatband(k, 3.8 + k ** 2, fill=mauve, cut=True, **style)
plot.fatband(k, 3.8 + 1.5 * np.cos(k), fill=mauve, cut=True, **style)
plot.fatband(k, 2.8 + (np.pi - k) ** 2, fill=mauve, cut=True, **style)
plot.fatband(k, 0.6 + 0.2 * np.cos(k), fill=orange, cut=True, **style)

for cut, color in zip(
        [(None, None, 0, None), (None, None, None, 0)],
        [orange, darkorange],
        ):
    plot.fatband(k, -0.6 - 1.0 * np.cos(k), fill=color, cut=cut, **style)
    plot.fatband(k, -0.2 + 1.0 * np.cos(k), fill=color, cut=cut, **style)

plot.fatband(k, -1.4 - 0.2 * np.cos(k), fill=darkorange, cut=True, **style)
plot.fatband(k, -2.5 + 0.5 * np.cos(k), fill=darkmauve, cut=True, **style)
plot.fatband(k, -2.5 - 0.5 * np.cos(k), fill=darkmauve, cut=True, **style)
plot.fatband(k, -7.2 - 0.1 * np.cos(k), fill=darkmauve, cut=True, **style)

for color, x, ymin, ymax, angle, symbol in zip(
        [darkorange, darkmauve, darkmauve, darkmauve],
        [0.3, 1.1, 1.5, 2.6],
        [-1.3, -0.9, -2.3, -6.9],
        [0.5, 3.5, 0.4, 5.8],
        [60, 80, 70, 85],
        'ARRR',
        ):
    plot.code(r'''
\draw[<->, very thick, color=%s] (<x=%g>, <y=%g>)
    to[out=%g, in=%g] node[right=-2pt] {$\boldsymbol{{\in}\,%s}$}
    (<x=%g>, <y=%g>);''' % (color, x, ymin, angle, -angle, symbol, x, ymax))

plot.save('fig01.pdf')
