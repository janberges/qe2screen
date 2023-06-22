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

q, x, (G, M, K, G) = elphmod.bravais.GMKG(12, mesh=True, corner_indices=True)
x = np.concatenate((x[K:] - x[K], x[-1] + x[1:M + 1] - x[K]))

def fixx(q, x):
    a, b = sorted(np.argsort(np.linalg.norm(q, axis=1))[:2])
    x[b:] -= x[b] - x[a]
    x[b:] += np.linalg.norm(q[a]) + np.linalg.norm(q[b])

qq, xq, eigq = elphmod.el.read_bands('eig.dat')
qn, xn, eign = elphmod.el.read_bands('eignol.dat')
q0, x0, eig0 = elphmod.el.read_bands('eig00.dat')

qq, xq, cosq = elphmod.el.read_bands('cos.dat')
qn, xn, cosn = elphmod.el.read_bands('cosnol.dat')

fixx(qq, xq)
fixx(qn, xn)

if comm.rank != 0:
    raise SystemExit

def new_plot(label, **kwargs):
    plot = storylines.Plot(
        style='APS',
        preamble=r'\let\Re\relax\DeclareMathOperator\Re{Re}',

        margin=margin,

        xticks=list(zip(x, [
            r'$\,\mathrm K$',
            '', '', '',
            r'$\Gamma$',
            '', '', '', '', '',
            r'$\mathrm M\,$',
            ])),

        line_cap='butt',

        **kwargs)

    plot.width = (Margin + 3 * margin - plot.single) / 2
    plot.height = (Margin + 2 * margin - plot.single) / 2

    plot.node(0, plot.ymax, '(%s)' % label, below_right=True)

    return plot

plot = new_plot(label='a', xlabels=False, left=Margin, top=Margin / 2,
    ylabel=r'$\mathcal G^{(1, 2)}$ (eV\textsuperscript3)',
    title=r'no $\mathcal L$', ymin=-0.0017, ymax=0.0027)

plot.line(x0, eig0[0], color='lightgray')
plot.line(y=0.0, color='lightgray')

for nu in range(len(eign)):
    plot.line(xn, eign[nu], thick=True, color=mauve)

plot.save('fig12a.pdf')

plot = new_plot(label='b', xlabels=False, ylabels=False, top=Margin / 2,
    title=r'$\mathbf Z^*$ and $Q$', ymin=-0.0017, ymax=0.0027)

plot.line(x0, eig0[0], color='lightgray')
plot.line(y=0.0, color='lightgray')

for nu in range(len(eign)):
    plot.line(xq, eigq[nu], thick=True, color=mauve)

plot.save('fig12b.pdf')

plot = new_plot(label='c', left=Margin, bottom=Margin / 2,
    ylabel=r'$\Re \bar{\mathbf g}^{\text p} \cdot \mathbf g '
        r'/ (|\mathbf g^{\text p}| |\mathbf g|)$',
    yticks=[(0, r'0\rlap{ \smash{($\perp$)}}'), 0.5,
        (1, r'1\rlap{ \smash{($\parallel$)}}')],
    ymin=0.0, ymax=1.2, lpos='cbbbt')

plot.line(y=1.0, color='lightgray',
    label=r'$\mathbf g^{\text p} \rightarrow \mathbf g$')

plot.line(xn, cosn[0], thick=True, color=orange)

plot.save('fig12c.pdf')

plot = new_plot(label='d', ylabels=False, bottom=Margin / 2,
    ymin=0.0, ymax=1.2, ystep=0.5)

plot.line(y=1.0, color='lightgray')

plot.line(xq, cosq[0], thick=True, color=orange)

plot.save('fig12d.pdf')

storylines.combine('fig12.pdf', ['fig12%s' % a for a in 'abcd'], columns=2)
