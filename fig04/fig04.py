#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

comm = elphmod.MPI.comm

labels = 'abcdef'

Margin = 0.95
margin = 0.12

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)

nbnd = [0, 1, 5, 13, 17, 22]

q, x, (G, M, K, G) = elphmod.bravais.GMKG(12, mesh=True, corner_indices=True)
x = np.concatenate((x[K:] - x[K], x[-1] + x[1:M + 1] - x[K]))

def fixx(q, x):
    a, b = sorted(np.argsort(np.linalg.norm(q, axis=1))[:2])
    x[b:] -= x[b] - x[a]
    x[b:] += np.linalg.norm(q[a]) + np.linalg.norm(q[b])

for i, n in enumerate(nbnd):
    q0, x0, g0 = elphmod.el.read_bands('ref%02d.dat' % n)
    qn, xn, gn = elphmod.el.read_bands('nol%02d.dat' % n)

    fixx(q0, x0)
    fixx(qn, xn)

    if n:
        qd, xd, gd = elphmod.el.read_bands('dpl%02d.dat' % n)
        qq, xq, gq = elphmod.el.read_bands('qpl%02d.dat' % n)

        fixx(qd, xd)
        fixx(qq, xq)

    if comm.rank != 0:
        continue

    plot = storylines.Plot(
        style='APS',
        preamble=r'\usepackage{mathtools}',

        width=1.28,
        height=1.6,

        left=0.5,
        right=0.15,
        bottom=0.5,
        top=0.12,

        xmin=0,
        xmax=25,
        ymin=-5,
        ymax=0,

        grid=True,

        xticks=[0, (25, r'\llap{25}\smash\AA')],
        yticks=[(0, 1)] + [(np.log10(y), None if n < 3 else '')
            for ticks in [
                [2.5e-1, 5e-1, 7.5e-1],
                [2.5e-2, 5e-2, 7.5e-2, 1e-1],
                [2.5e-3, 5e-3, 7.5e-3, 1e-2],
                [2.5e-4, 5e-4, 7.5e-4, 1e-3],
                [2.5e-5, 5e-5, 7.5e-5, 1e-4]]
            for n, y in enumerate(ticks)] + [(-5, '$10\mathrlap{^{-5}}$')],
        )

    plot.axes()

    for label, color, mark in [
            ('nol', 'gray', '*'),
            ('dpl', orange, '+'),
            ('qpl', mauve, 'x')]:

        d, g = np.loadtxt('%s%02d.decay' % (label, n), skiprows=1).T

        unique = []

        for group in elphmod.misc.group(np.array(list(zip(d, g))), eps=1e-7):
            D = np.average(d[group], axis=0)
            G = np.average(g[group], axis=0)
            unique.append([D, G])

        unique = np.array(unique)

        plot.line(unique[:, 0], np.log10(unique[:, 1]), color=color, cut=True,
            mark=mark, mark_size='0.6pt', only_marks=True)

        if n == 0:
            break

    plot.save('fig04%s.in.pdf' % labels[i])

    plot = storylines.Plot(
        style='APS',
        preamble=r'\usepackage{mathtools}\usepackage{bm}',

        left=Margin,
        right=margin,
        bottom=0.4,
        top=0.4,

        title=r'\smash{(%s)} %d band' % (labels[i], n),

        grid=True,

        ymin=0.0,
        ymax=0.64,
        ystep=0.2,

        xticks=zip(x, [
            r'$\,\mathrm K$',
            None, None, None,
            r'$\Gamma$',
            None, None, None, None, None,
            r'$\mathrm M\,$',
            ]),

        ylabel=r'El.-ph.\@ coupling (eV\textsuperscript{3/2})',

        lpos='rt',
        llen='2mm',
        lopt='below left=0.5mm, inner sep=0.6mm, rounded corners=1pt, '
            'draw=lightgray, fill=white',
        )

    plot.width = (Margin + margin * (2 * len(nbnd) - 1)
        - plot.double) / len(nbnd)

    plot.height = plot.width

    plot.axes()

    if n != 1:
        plot.title += 's'

    if i != 0:
        plot.left = margin
        plot.ylabel = None
        plot.ymarks = False

    for nu in range(len(gn)):
        plot.line(xn, gn[nu], color='gray')

    if n:
        for nu in range(len(gd)):
            plot.line(xd, gd[nu], color=orange, thick=True)

        for nu in range(len(gq)):
            plot.line(xq, gq[nu], color=mauve, thick=True)
    else:
        plot.line(label='*next*', color='gray',
            mark='*', mark_size='0.6pt', only_marks=True)
        plot.line(label=r'no\,$\mathcal L$', color='gray')

        plot.line(label='*next*', color=orange,
            mark='+', mark_size='0.6pt', only_marks=True)
        plot.line(label=r'$\bm Z^*$', color=orange, thick=True)

        plot.line(label='*next*', color=mauve,
            mark='x', mark_size='0.6pt', only_marks=True)
        plot.line(label=r'$\bm Z\mathrlap{^*}, Q$', color=mauve, thick=True)

    for nu in range(len(g0)):
        plot.line(x0, g0[nu], mark='*', mark_size='0.3pt', only_marks=True)

    plot.node(0, plot.ymax, r'\includegraphics{fig04%s.in.pdf}' % labels[i],
        below_right='-0.6pt', inner_sep='0.3mm', rounded_corners='1pt',
        draw='gray', fill='white', thick=True)

    plot.save('fig04%s.pdf' % labels[i])

if comm.rank == 0:
    storylines.combine('fig04.pdf', ['fig04%s' % a for a in labels])
