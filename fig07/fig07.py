#!/usr/bin/env python3

import elphmod
import storylines
import numpy as np

comm = elphmod.MPI.comm

Margin = 0.96
margin = 0.12

white = storylines.Color(255, 255, 255)
orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

cmap_negative = storylines.colormap(
    (0.0, darkmauve),
    (0.5, mauve),
    (1.0, white))

cmap_positive = storylines.colormap(
    (0.0, white),
    (0.5, orange),
    (1.0, darkorange))

cmap = storylines.colormap(
    (0.00, darkmauve),
    (0.25, mauve),
    (0.50, white),
    (0.75, orange),
    (1.00, darkorange))

labels = 'abcd'

#for label in '0r', 'zr', '0cr', 'zcr', 'nol0cr', 'nolzcr', 'dr':
#    qA, xA, A = elphmod.el.read_bands('phband%sw1.freq' % label)
#
#    A *= 2 / (1e3 * elphmod.misc.Ry) ** 2 # correct units from old data
#    A = A[::-1]
#
#    if not 'z' in label:
#        elphmod.el.write_bands('phband%sw1.dat' % label, qA,
#            50 * A.sum(axis=0, keepdims=True) / A.shape[0])
#
#    if label in {'0cr', 'zcr'}:
#        A += 1
#        A /= 2
#
#    A = np.round(255 * np.minimum(np.maximum(A, 0), 1))
#
#    storylines.save('phband%sw1.png' % label, A[:, :, None])

A = dict()
integral = dict()
xw = dict()
w = dict()

for label in '0r', 'zr', '0cr', 'zcr', 'nol0cr', 'nolzcr', 'dr':
    qw, xw[label], w[label] = elphmod.el.read_bands('phband%s1.freq' % label)

    A[label] = np.array(storylines.load('phband%sw1.png' % label))[:, :, 0]

    if not 'z' in label:
        integral[label] = elphmod.el.read_bands('phband%sw1.dat' % label)[-1][0]

q, x, corners = elphmod.bravais.path('GMKG', ibrav=4)
qz, xz, cornerz = elphmod.bravais.path('KGM', ibrav=4)
qz /= 15
xz /= 15

for n, (label, zoom, title) in enumerate([
        ('0r', 'zr', r'Via $\varPi^{\text{00}}$'),
        ('0cr', 'zcr', r'Via $\varPi^{\text{p0}}$'),
        ('nol0cr', 'nolzcr', r'Via $\varPi^{\text{p0}}$ '
            r'\smash{(no\,$\mathcal L$)}'),
        ('dr', None, r'Via $\varPi^{\text{00}}$ \smash{(doped)}')]):

    image = elphmod.plot.color(A[label], minimum=0, maximum=255,
        cmap=cmap if label == '0cr' else cmap_positive)

    if zoom is not None:
        inset = elphmod.plot.color(A[zoom], minimum=0, maximum=255,
            cmap=cmap if label == '0cr' else cmap_positive)

    if comm.rank == 0:
        if zoom is not None:
            plot = storylines.Plot(
                style='APS',
                preamble=r'\usepackage{mathtools}',

                width=-1.3,
                height=-1.3,

                left=0.5,
                right=0.2,
                bottom=0.55,
                top=0.2,

                xticks=list(zip(xz[cornerz], [
                    r'$\smash{\mathllap{\frac 1 {15}}}\mathrm K$',
                    r'$\Gamma$',
                    r'$\smash{\mathllap{\frac 1 {15}}}\mathrm M$',
                    ])),

                ymin=25.0,
                ymax=50.0,
                ystep=10.0,

                background='fig07%s.in.png' % labels[n],
                )

            for nu in range(w[zoom].shape[0]):
                plot.line(xw[zoom], w[zoom][nu], cut=True)

            storylines.save(plot.background, inset)

            plot.save('fig07%s.in.pdf' % labels[n])

        plot = storylines.Plot(
            style='APS',

            title='\smash{(%s)} %s' % (labels[n], title),

            left=Margin,
            right=Margin,
            bottom=margin,
            top=Margin / 2,

            xticks=list(zip(x[corners], [
                r'$\Gamma$',
                r'$\mathrm M$',
                r'$\mathrm K$',
                r'$\Gamma$',
                ])),
            xlabels=False,

            ylabel='Phonon energy (meV)',
            yformat=lambda y: '$%g\,\mathrm i$' % abs(y)
                if y < 0 else '$%g$' % y,

            zlabel='Phonon spectral function (1/meV)',
            zclose=True,

            ymin=-25.0,
            ymax=50.0,
            ystep=10.0,
            zmin=0.0,
            zmax=1.0,
            zstep=1.0,

            cmap=cmap_positive,

            lpos='bl',
            lopt='above right=1mm',
            )

        plot.width = (2 * Margin + 6 * margin - plot.double) / 4
        plot.height = (Margin + 2 * margin - plot.single) * 0.9

        plot.image('fig07%s.png' % labels[n], 0.0, 0.0, x[-1], 50.0)

        plot.line(y=0, color='lightgray')

        for nu in range(w[label].shape[0]):
            plot.line(xw[label], w[label][nu], cut=True,
                label=None if n < 3 else r'$\omega = 0$')

        if n < 3:
            plot.colorbar = False
            plot.right = margin

        if n > 0:
            plot.ylabel = None
            plot.ylabels = False
            plot.left = margin

        if zoom is not None:
            plot.node(x[-1], plot.ymin, r'\includegraphics{fig07%s.in.pdf}'
                % labels[n], above_left=True)

        storylines.save('fig07%s.png' % labels[n], image)

        plot.save('fig07%s.pdf' % labels[n])

        plot.clear()

        plot.height = (Margin + 2 * margin - plot.single) * 0.1

        plot.title = None
        plot.bottom = Margin / 2
        plot.top = margin

        plot.xlabels = True

        if n == 0:
            plot.ylabel = 'Integral'

        plot.ymax = 9.5
        plot.ymin = 3.0
        plot.ystep = 3.0

        plot.cmap = cmap_negative
        plot.zmin = -1.0
        plot.zmax = 0.0
        plot.zticks = [-1]
        plot.zlabel = None

        for y in 8.0, 9.0:
            plot.line(y=y, color='lightgray')

        plot.line(xw[label], integral[label], color=darkorange, cut=True)

        plot.save('fig07%s.b.pdf' % labels[n])

if comm.rank == 0:
    storylines.combine('fig07.pdf', ['fig07%s%s' % (a, b)
        for b in ['', '.b'] for a in labels], columns=4)
