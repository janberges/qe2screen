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

labels = 'abcdef'

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

        width=2.0,
        height=2.5,

        margin=margin,
        left=Margin / 2,
        bottom=Margin / 2,

        xmin=0.0,
        ymin=-6.0,
        ymax=0.0,

        xticks=[0, (20, '2~nm')],
        yticks=[(0, 1),
            (-3, '$10\mathrlap{^{-3}}$'), (-6, '$10\mathrlap{^{-6}}$')],
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

        margin=margin,
        left=Margin / 2,
        bottom=Margin / 2,

        ymin=0.0,
        ymax=0.16 if i < 3 else 0.64,
        ystep=0.05 if i < 3 else 0.2,

        xticks=list(zip(x, [
            r'$\,\mathrm K$',
            '', '', '',
            r'$\Gamma$',
            '', '', '', '', '',
            r'$\mathrm M\,$',
            ])),

        lpos='lctm',
        )

    plot.width = (Margin + 5 * margin - plot.double) / 3
    plot.height = (Margin / 2 + 3 * margin - height) / 2

    plot.axes()

    if i % 3:
        plot.left = margin
        plot.ylabel = None
        plot.ylabels = False

    if i < 3:
        plot.bottom = margin
        plot.xlabels = False

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
        plot.line(label=r'No\,$\mathcal L$', color='gray')

        plot.line(label='*next*', color=orange,
            mark='+', mark_size='0.6pt', only_marks=True)
        plot.line(label=r'$\mathbf Z^*$ only', color=orange, thick=True)

        plot.line(label='*next*', color=mauve,
            mark='x', mark_size='0.6pt', only_marks=True)
        plot.line(label=r'$\mathbf Z^*$ and $Q$', color=mauve, thick=True)

    for nu in range(len(g0)):
        plot.line(x0, g0[nu], mark='*', mark_size='0.3pt', only_marks=True)

    plot.node(x[0], plot.ymax, r'(%s) %d band%s'
        % (labels[i], n, '' if n == 1 else 's'), below_right='0.5mm')

    plot.node(x[-1], plot.ymax, r'\includegraphics{fig04%s.in.pdf}' % labels[i],
        below_left=True)

    plot.save('fig04%s.pdf' % labels[i])

if comm.rank == 0:
    plot = storylines.Plot(
        style='APS',

        width=Margin / 2,
        height=5.0,

        margin=0.0,
        xyaxes=False,
        )

    plot.code(r'\node [rotate=90, above=\baselineskip] at (%g, <y=0>) '
        r'{Electron-phonon coupling (eV\textsuperscript{3/2})};'
            % (Margin - plot.tick))

    plot.save('fig04.l.pdf')

    storylines.combine('fig04.r.pdf', ['fig04%s' % a for a in labels],
        columns=3)

    storylines.combine('fig04.pdf', ['fig04.l', 'fig04.r'], align=0.5)
