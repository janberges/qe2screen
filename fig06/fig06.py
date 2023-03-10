#!/usr/bin/env python3

import elphmod
import storylines

nbnd = [1, 5, 13, 17]

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

q, x, GMKG = elphmod.bravais.path('GMKG', ibrav=4)

q0, x0, w0 = elphmod.el.read_bands('../fig05/ref.freq')

w0 *= 1e3 * elphmod.misc.cmm1

qr, xr, wr = elphmod.el.read_bands('../fig05/phband01cr1.freq')
q1, x1, w1 = elphmod.el.read_bands('phband01corr1cr1.freq')
q2, x2, w2 = elphmod.el.read_bands('phband01corr2cr1.freq')

if elphmod.MPI.comm.rank != 0:
    raise SystemExit

plot = storylines.Plot(
    style='APS',

    width=5,
    height=2.5,

    margin=0.0,
    xyaxes=False,

    grid=True,

    xmin=0.17 * x[-1],
    xmax=0.43 * x[-1],
    ymin=-13.0,
    ymax=-3.5,
    ystep=10.0,

    xticks=list(zip(x[GMKG], [
        r'$\Gamma$',
        r'$\mathrm M$',
        r'$\mathrm K$',
        r'$\Gamma$',
        ])),

    ylabel='Phonon energy (meV)',

    lpos='bcr',
    llen='4mm',
    lput=False,
    lwid=3,
    lopt='above=2mm, inner ysep=3pt, rounded corners=1pt, '
        'draw=lightgray, fill=white',

    yformat=lambda y: '$%g\,\mathrm i$' % abs(y) if y < 0 else '$%g$' % y,

    )

plot.axes()

for nu in range(len(wr)):
    plot.line(xr, wr[nu], thick=True, color=mauve, cut=True,
        label='no correction')

for nu in range(len(w1)):
    plot.line(x1, w1[nu], thick=True, color=darkorange, cut=True,
        dash_pattern='on 0.8mm off 0.8mm', line_cap='butt', line_join='miter',
        label=r'correction I')

for nu in range(len(w2)):
    plot.line(x2[1:-1], w2[nu, 1:-1], thick=True, color=orange, cut=True,
        dash_pattern='on 0.5mm off 0.2mm', line_cap='butt', line_join='miter',
        label=r'correction II')

for nu in range(len(w0)):
    plot.line(x0, w0[nu], cut=True,
        mark='*', mark_size='0.5pt', only_marks=True)

plot.save('fig06.in.pdf')

plot.code(r'\draw [gray, thick, rounded corners=1pt] '
    '(<x=%g>, <y=%g>) rectangle (<x=%g>, <y=%g>);'
        % (plot.xmin, plot.ymin, plot.xmax, plot.ymax), zindex=1)

width = plot.single
height = 5.795

left = 0.9
right = 0.12
bottom = 0.4
top = 0.4

xscale = (plot.xmax - plot.xmin) / plot.width
yscale = (plot.ymax - plot.ymin) / plot.height

xscale_new = (x[-1] - x[0]) / (width - left - right)
yscale_new = xscale_new / xscale * yscale

plot.xmin = x[0]
plot.xmax = x[-1]
plot.ymax = plot.ymin + yscale_new * (height - bottom - top)
plot.ymin = -14.0

plot.width = width
plot.height = height

plot.left = left
plot.right = right
plot.bottom = bottom
plot.top = top

plot.xaxis = True
plot.yaxis = True
plot.frame = True

plot.lput = True

plot.node((x[0] + x[-1]) / 2, plot.ymax + yscale_new * top,
    r'\includegraphics{fig06.in.pdf}', below=True, inner_sep='0pt',
    rounded_corners='1pt', draw='gray', fill='white', thick=True)

plot.save('fig06.pdf')
