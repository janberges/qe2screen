#!/usr/bin/env python3

import copy
import elphmod
import numpy as np
import storylines

comm = elphmod.MPI.comm

Margin = 0.96
margin = 0.12

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1) / 2
darkmauve = storylines.Color(94, 60, 153) / 2

colors = [(darkmauve, mauve), (darkorange, orange)]

nel = [55, 25]

thickness = 0.05
factor = 1.5

ylower = 375.0
ybreak = 2350.0

q, x, (G, M, K, G) = elphmod.bravais.GMKG(12, mesh=True, corner_indices=True)
x = np.concatenate((x[K:] - x[K], x[-1] + x[1:M + 1] - x[K]))

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        label='a',

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        ymin=280.0,
        ymax=2820.0 - ybreak,

        xticks=list(zip(x, [
            r'$\mathrm K$',
            None, None, None,
            r'$\Gamma$',
            None, None, None, None, None,
            r'$\mathrm M\,$',
            ])),

        yticks=[300, 350, 2750, 2800],

        lpos='cm',
        llen='4mm',

        ylabel='Bare phonon energy (meV)',

        resolution=1e-4,
        )

    plot.yticks = [y if y < ylower else (y - ybreak, y) for y in plot.yticks]

    plot.width = plot.double / 3
    plot.height = plot.single

    plot.cut(y=ylower)

def fixx(q, x):
    a, b = sorted(np.argsort(np.linalg.norm(q, axis=1))[:2])
    x[b:] -= x[b] - x[a]
    x[b:] += np.linalg.norm(q[a]) + np.linalg.norm(q[b])

for i, n in enumerate(nel):
    q0, x0, w0 = elphmod.el.read_bands('phband%dr.freq' % n)
    qw, xw, w = elphmod.el.read_bands('phband%d.freq' % n)
    ql, xl, l = elphmod.el.read_bands('phband%dl.freq' % n)

    fixx(q0, x0)
    fixx(qw, xw)
    fixx(ql, xl)

    pol = np.empty(l.shape + (2,))
    pol[:, :, 0] = l
    pol[:, :, 1] = 1 - l

    w0[np.where(w0 > ylower)] -= ybreak
    w[np.where(w > ylower)] -= ybreak

    if comm.rank != 0:
        continue

    for nu in range(len(w)):
        plot.compline(xw, w[nu], pol[nu], colors=colors[i], thickness=thickness,
            cut=True)

    for nu in range(len(w0)):
        plot.line(x0, w0[nu], mark='*', mark_size='0.3pt', only_marks=True,
            cut=True)

if comm.rank == 0:
    inset = copy.deepcopy(plot)

    inset.width = 4.0

    inset.left = 0.0
    inset.right = 0.0
    inset.bottom = 0.0
    inset.top = 0.0

    inset.xaxis = False
    inset.yaxis = False
    inset.frame = False

    inset.xmin = 0.45
    inset.xmax = x[4]
    inset.ymin = 308.0
    inset.ymax = 316.0

    scale = inset.width / ((inset.xmax - inset.xmin)
        / (x[-1] - x[0]) * (plot.width - plot.left - plot.right))

    inset.height = scale * ((inset.ymax - inset.ymin)
        / (plot.ymax - plot.ymin) * (plot.height - plot.bottom - plot.top))

    inset.lpos = 'rb'
    inset.lopt = 'above left'

    inset.line(color=colors[1][0], label=r'Long.',
        line_width='%gcm' % (factor * thickness))

    inset.line(color=colors[1][1], label=r'Trans.',
        line_width='%gcm' % (factor * thickness))

    inset.label = None

    for line in inset.lines:
        if line['weights'] is not None:
            for prop in 'weights', 'shifts':
                line[prop] = np.array(line[prop]) * factor

        if 'mark_size' in line['options']:
            line['options']['mark_size'] = '%gpt' % (factor
                * float(line['options']['mark_size'][:-2]))

    inset.save('fig10a.in.pdf')

    plot.code(r'\draw [lightgray, ounded corners=1pt] '
        '(<x=%g>, <y=%g>) rectangle (<x=%g>, <y=%g>);'
            % (inset.xmin, inset.ymin, inset.xmax, inset.ymax))

    plot.node((x[-1] + x[0]) / 2, 384, r'\includegraphics{fig10a.in.pdf}',
        inner_sep='0.2pt', rounded_corners='1pt', draw='lightgray')

    plot.save('fig10a.pdf')

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        label='b',

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        ymin=0.0,
        ymax=2.25,
        ystep=0.5,

        xticks=list(zip(x, [
            r'$\mathrm K$',
            None, None, None,
            r'$\Gamma$',
            None, None, None, None, None,
            r'$\mathrm M\,$',
            ])),

        ylabel=r'Electron-phonon coupling (eV\textsuperscript{3/2})',

        resolution=1e-4,
        )

    plot.width = plot.double / 3
    plot.height = plot.single

    plot.axes()

for i, n in enumerate(nel):
    q0, x0, g0 = elphmod.el.read_bands('coupling%dr.dat' % n)
    qg, xg, g = elphmod.el.read_bands('coupling%d.dat' % n)
    ql, xl, l = elphmod.el.read_bands('coupling%dl.dat' % n)

    fixx(q0, x0)
    fixx(qw, xw)
    fixx(ql, xl)

    pol = np.empty(l.shape + (2,))
    pol[:, :, 0] = l
    pol[:, :, 1] = 1 - l

    if comm.rank != 0:
        continue

    for nu in range(len(g)):
        plot.compline(xg, g[nu], pol[nu], colors=colors[i], cut=True)

    plot.compline(xg, -g[-1], pol[-1], colors=colors[i], cut=True)

    for nu in range(len(g0)):
        plot.line(x0, g0[nu], mark='*', mark_size='0.3pt', only_marks=True,
            cut=True)

if comm.rank == 0:
    inset = copy.deepcopy(plot)

    inset.width = 2.0

    inset.left = 0
    inset.right = 0
    inset.bottom = 0
    inset.top = 0

    inset.xaxis = False
    inset.yaxis = False
    inset.frame = False

    inset.xmin = 0.45
    inset.xmax = x[4]
    inset.ymin = 0.0
    inset.ymax = 0.3

    scale = inset.width / ((inset.xmax - inset.xmin)
        / (x[-1] - x[0]) * (plot.width - plot.left - plot.right))

    inset.height = scale * ((inset.ymax - inset.ymin)
        / (plot.ymax - plot.ymin) * (plot.height - plot.bottom - plot.top))

    inset.label = None

    for line in inset.lines:
        if line['weights'] is not None:
            for prop in 'weights', 'shifts':
                line[prop] = np.array(line[prop]) * factor

        if 'mark_size' in line['options']:
            line['options']['mark_size'] = '%gpt' % (factor
                * float(line['options']['mark_size'][:-2]))

    inset.save('fig10b.in.pdf')

    plot.code(r'\draw [lightgray, rounded corners=1pt] '
        '(<x=%g>, <y=%g>) rectangle (<x=%g>, <y=%g>);'
            % (inset.xmin, inset.ymin, inset.xmax, inset.ymax))

    plot.node(x[0], plot.ymax, r'\includegraphics{fig10b.in.pdf}',
        below_right='1mm', inner_sep='0.2pt', rounded_corners='1pt',
        draw='lightgray')

    plot.node(x[0], plot.ymax, r'\phantom{\includegraphics{fig10b.in.pdf}}',
        below_right='1mm', inner_sep='0.2pt', rounded_corners='1pt',
        draw='lightgray')

    plot.save('fig10b.pdf')

q, x, corners = elphmod.bravais.path('GMKG', ibrav=4, N=198)

q0, x0, w0 = elphmod.el.read_bands('ref.freq')

w0 *= 1e3 * elphmod.misc.cmm1

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        label='c',

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        ymin=-60.0,
        ymax=+50.0,
        ystep=20.0,

        xticks=list(zip(x[corners], [
            r'$\Gamma$',
            r'$\mathrm M$',
            r'$\mathrm K$',
            r'$\Gamma$',
            ])),

        ylabel='Renormalized phonon energy (meV)',

        yformat=lambda y: '$%g\,\mathrm i$' % abs(y) if y < 0 else '$%g$' % y,

        lpos='bc',
        lopt='above=1mm',
        )

    plot.line(y=0.0, color='lightgray')
    plot.line(line_width='%gcm' % thickness, color=orange, label='Without SC')
    plot.line(line_width='%gcm' % thickness, color=mauve, label='With SC')

    plot.width = plot.double / 3
    plot.height = plot.single

    plot.axes()

for i, n in enumerate(nel):
    qw, xw, w = elphmod.el.read_bands('phband%dcr1.freq' % n)
    ql, xl, l = elphmod.el.read_bands('phband%dcr1l.freq' % n)

    pol = np.empty(l.shape + (2,))
    pol[:, :, 0] = l
    pol[:, :, 1] = 1 - l

    if comm.rank != 0:
        continue

    for nu in range(len(w)):
        plot.compline(xw[1:-1], w[nu, 1:-1], pol[nu, 1:-1], colors=colors[i],
            protrusion=1, cut=True)

    for nu in range(len(w0)):
        plot.line(x0, w0[nu], mark='*', mark_size='0.3pt', only_marks=True,
            cut=True)

if comm.rank == 0:
    plot.save('fig10c.pdf')

    storylines.combine('fig10.pdf', ['fig10%s' % a for a in 'abc'])
