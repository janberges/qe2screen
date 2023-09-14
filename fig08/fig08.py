#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

comm = elphmod.MPI.comm

Margin = 0.96
margin = 0.12

orange = storylines.Color(241, 163, 64)
mauve = storylines.Color(153, 142, 195)
darkorange = storylines.Color(230, 97, 1)
darkmauve = storylines.Color(94, 60, 153)

K = np.array([-2.0 / 3.0, 1.0 / 3.0, 0.0])

k, x, corners = elphmod.bravais.path([n * K / 4 for n in [5, 4, 2, 1]],
    ibrav=4, N=198)

#a = elphmod.bravais.primitives(ibrav=4)
#b = elphmod.bravais.reciprocals(*a)
#k_cart = elphmod.bravais.crystal_to_cartesian(k / (2 * np.pi), *b)
#
#el = elphmod.el.Model('MoS2')
#
#e, order = elphmod.dispersion.dispersion(el.H, k, order=True)
#e -= elphmod.el.read_Fermi_level('scf.out')
#
#elphmod.el.write_bands('bands.dat', k_cart, e[:, 1:2].T)
#
#ph = elphmod.ph.Model('dyn', apply_asr_simple=True)
#
#w2, u = np.linalg.eigh(ph.D(0, np.pi, 0))
#
#for infile, outfile in ('MoS2.epmatwp', 'd0.dat'), ('MoS2.epmatwpc', 'dp.dat'):
#    elph = elphmod.elph.Model(infile, 'MoS2.wigner', el, ph)
#
#    g = [abs(np.einsum('xmn,xv->vmn', elph.g(0, np.pi, 0, k1, k2, k3,
#        elbnd=True), u)[1, 1, 1]) for k1, k2, k3 in k]
#
#    g = np.array(g) * elphmod.misc.Ry ** 1.5
#
#    elphmod.el.write_bands(outfile, k_cart, g[np.newaxis], fmt='%8.6f')

ke, xe, (e,) = elphmod.el.read_bands('bands.dat')

k0, x0, (g0,) = elphmod.el.read_bands('coupling0.dat')
kp, xp, (gp,) = elphmod.el.read_bands('couplingp.dat')

for l, (label, ylabel, ymin, ymax, ystep) in enumerate([
        ('a', r'Electron energy (eV)', -0.17, 0.37, 0.1),
        ('c', r'El.-ph.\@ coupling (eV$^{3/2}$)', 0.0, 0.033, 0.01)]):

    if comm.rank != 0:
        continue

    plot = storylines.Plot(
        style='APS',
        preamble=r'\usepackage{mathtools, nicefrac}',

        label=label,

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        ymin=ymin,
        ymax=ymax,
        ystep=ystep,

        xticks=list(zip(x[corners], [
            r'$\mathllap{\smash{\nicefrac 5 4}}\,\mathrm K$',
            r'$\mathrm K$',
            r'$\mathllap{\smash{\nicefrac 1 2}\,\mathrm K = {}} \mathrm Q$',
            r'$\mathllap{\smash{\nicefrac 1 4}}\,\mathrm K$',
            ])),

        ylabel=ylabel,

        lpos=(x[-1] / 2, 0.01),
        )

    plot.width = (2 * Margin + 2 * margin - plot.single) / 2
    plot.height = (Margin / 2 + 3 * margin - plot.single) / 2

    if l == 0:
        plot.bottom = margin
        plot.xlabels = False

        plot.line(y=0.0, color='lightgray')

        plot.line(xe[::-1], e, color='gray', densely_dashed=True, cut=True)
        plot.fatband(xe, e, thickness=0.05, fill=mauve, protrusion=1, cut=True)
        plot.fatband(xe, e, thickness=0.05, fill=darkmauve,
            cut=(None, None, None, 0))

        c = (x[corners[2]] + x[corners[1]]) / 2
        d = (x[corners[2]] - x[corners[1]]) / 2 * 0.85

        plot.code(r'\draw[<->, thick] '
            r'(<x=%g>, <y=-0.1>) '
            r'to node[below, align=center] '
            r'{$\mathbf k\!\leftrightarrow\!\mathbf k\!+\!\mathrm M$} '
            r'(<x=%g>, <y=-0.1>);' % (c - d, c + d))
    else:
        for y, label, color in (g0, 'DFPT', orange), (gp, 'cDFPT', darkorange):
            plot.fatband(x0, y, thickness=0.03, fill=color, protrusion=1,
                cut=True)

            plot.line(label=label, color=color, line_width='0.3mm')

        plot.ltop = r'$\mathbf q = \mathrm M$'

    plot.save('fig08%s.pdf' % plot.label)

q, x, corners = elphmod.bravais.path('GMKG', ibrav=4, N=198)

q0, x0, w0 = elphmod.el.read_bands('ref.freq')

w0 *= 1e3 * elphmod.misc.cmm1

qi, xi, wi = elphmod.el.read_bands('phband.freq')
qc, xc, wc = elphmod.el.read_bands('phbandc.freq')
qu, xu, wu = elphmod.el.read_bands('phbandr.freq')
qw, xw, ww = elphmod.el.read_bands('phbandcr.freq')
q1, x1, w1 = elphmod.el.read_bands('phbandr1.freq')
q2, x2, w2 = elphmod.el.read_bands('phbandcr1.freq')

if comm.rank != 0:
    raise SystemExit

for l, label in enumerate('bd'):
    plot = storylines.Plot(
        style='APS',

        label=label,

        margin=margin,
        left=Margin,
        bottom=Margin / 2,

        ymin=0,
        ymax=80,
        ystep=20,

        dtop=2 * margin,

        xticks=list(zip(x[corners], [
            r'$\Gamma$',
            r'$\mathrm M$',
            r'$\mathrm K$',
            r'$\Gamma$',
            ])),

        yformat=lambda y: '$%g\,\mathrm i$' % abs(y) if y < 0 else '$%g$' % y,

        ylabel='Phonon energy (meV)',

        lpos='rt',
        lcol=2 if l == 0 else 1,
        lwid=2.8,
        lopt='below left',
        )

    plot.width = (2 * Margin + 2 * margin - plot.single) / 2
    plot.height = (Margin / 2 + 3 * margin - plot.single) / 2

    if l == 0:
        plot.bottom = margin
        plot.xlabels = False

        style = dict(line_cap='butt', line_join='miter')

        for nu in range(len(wi)):
            plot.line(xi, wi[nu], thick=True, color=orange, **style)

        for nu in range(len(ww)):
            plot.line(xw, ww[nu], thick=True, color=mauve,
                dash_pattern='on 0.8mm off 0.8mm', **style)

        for nu in range(len(wc)):
            plot.line(xc, wc[nu], thick=False, color=mauve, **style)

        for nu in range(len(wu)):
            plot.line(xu, wu[nu], thick=False, color=orange,
                dash_pattern='on 0.5mm off 0.2mm', **style)

        plot.line(thick=True, color=orange, label='$D$', **style)
        plot.line(thick=False, color=mauve, label=r'$D^{\text p}$', **style)
        plot.line(thick=False, color=orange, dash_pattern='on 0.5mm off 0.2mm',
            label=r'$D^{\text u}$', **style)
        plot.line(thick=True, color=mauve, dash_pattern='on 0.8mm off 0.8mm',
            label=r'$D^{\text p} {+} \varPi^{\text{p0}}$', **style)
    else:
        for nu in range(len(w1)):
            plot.line(x1, w1[nu], thick=True, color=orange, cut=True,
                label=r'Via $\varPi^{\text{00}}$')

        for nu in range(len(w2)):
            plot.line(x2, w2[nu], thick=True, color=mauve, cut=True,
                label=r'Via $\varPi^{\text{p0}}$')

        for nu in range(len(w0)):
            plot.line(x0, w0[nu], cut=True,
                mark='*', mark_size='0.3pt', only_marks=True)

    plot.save('fig08%s.pdf' % plot.label)

storylines.combine('fig08.pdf', ['fig08%s' % a for a in 'abcd'], columns=2)
