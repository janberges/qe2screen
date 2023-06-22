#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

comm = elphmod.MPI.comm

Margin = 0.96
margin = 0.12

labels = ['abc', 'def']

nk = 96
kT = 300.0 * elphmod.misc.kB

q = np.array([[0.0, np.pi], [2 * np.pi / 3, 2 * np.pi / 3]])
Q = 'MK'

cmap = storylines.colormap(
    (0.00, storylines.Color(94, 60, 153)),
    (0.25, storylines.Color(153, 142, 195)),
    (0.75, storylines.Color(241, 163, 64)),
    (1.00, storylines.Color(230, 97, 1)),
    (None, storylines.Color(255, 255, 255)),
    )

BZ = dict(points=500, outside=np.nan)

#el = elphmod.el.Model('MoS2')
#mu = elphmod.el.read_Fermi_level('scf.out')
#ph = elphmod.ph.Model('dyn', apply_asr_simple=True)
#elph = elphmod.elph.Model('MoS2.epmatwp', 'MoS2.wigner', el, ph)
#
#e, U = elphmod.dispersion.dispersion_full_nosym(el.H, nk, vectors=True)
#e -= mu
#
#d02 = np.zeros((nk, nk, el.size))
#
#for a in 0, 3, 4:
#    d02 += abs(U[:, :, a, :]) ** 2
#
#for k1 in range(nk):
#    for k2 in range(nk):
#        order = np.argsort(d02[k1, k2, 1:] < 0.5) + 1
#
#        e[k1, k2, 1:] = e[k1, k2, order]
#
#        for a in range(el.size):
#            U[k1, k2, a, 1:] = U[k1, k2, a, order]
#
#e = e[..., 1:2].copy()
#U = U[..., 1:2].copy()
#
#w2, u = elphmod.dispersion.dispersion(ph.D, q, vectors=True)
#u = u[range(2), :, [1, 2], np.newaxis]
#
#g = abs(elph.sample(q, U=U, u=u))
#g *= elphmod.misc.Ry ** 1.5
#
#if comm.rank == 0:
#    with open('mesh.dat', 'w') as data:
#        for k1 in range(nk):
#            for k2 in range(nk):
#                data.write('%2d %2d %6.3f %6.4f %6.4f\n'
#                    % (k1, k2, e[k1, k2, 0],
#                    g[0, 0, k1, k2, 0, 0],
#                    g[1, 0, k1, k2, 0, 0]))

e = np.empty((nk, nk, 1))
g = np.empty((len(q), 1, nk, nk, 1, 1))

if comm.rank == 0:
    with open('mesh.dat') as data:
        for k1 in range(nk):
            for k2 in range(nk):
                (K1, K2, e[k1, k2, 0],
                g[0, 0, k1, k2, 0, 0],
                g[1, 0, k1, k2, 0, 0]) = tuple(map(float, next(data).split()))

comm.Bcast(e)
comm.Bcast(g)

g2 = g ** 2

Pi = elphmod.diagrams.phonon_self_energy(q, e[..., :1], g2=g2[:, :1], kT=kT,
    occupations=elphmod.occupations.fermi_dirac, fluctuations=True)[1]

X0 = elphmod.diagrams.phonon_self_energy(q, e[..., :1], kT=kT,
    occupations=elphmod.occupations.fermi_dirac, fluctuations=True)[1]

Pi /= Pi.min()
X0 /= X0.min()
g2 /= g2.max()

a = elphmod.bravais.translations()
b = elphmod.bravais.reciprocals(*a)

for iq in range(len(q)):
    FS_kk = [np.dot(c, b) for c in elphmod.dos.isoline(e[:, :, 0])(0)]
    FS_kq = [np.dot(c, b) for c in elphmod.dos.isoline(np.roll(np.roll(e,
        shift=-int(round(q[iq, 0] * nk / (2 * np.pi))), axis=0),
        shift=-int(round(q[iq, 1] * nk / (2 * np.pi))), axis=1)[:, :, 0])(0)]

    kxmax, kymax, kx, ky = elphmod.plot.toBZ(return_only_k=True, **BZ)

    Pi_BZ = elphmod.plot.toBZ(Pi[iq, 0, :, :, 0, 0], **BZ)
    X0_BZ = elphmod.plot.toBZ(X0[iq, 0, :, :, 0, 0], **BZ)
    g2_BZ = elphmod.plot.toBZ(g2[iq, 0, :, :, 0, 0], **BZ)

    if comm.rank == 0:
        for i, (data, title) in enumerate([
                (Pi_BZ, r'$-\varPi^{\text{00}}$'),
                (X0_BZ, r'$-\chi^{\text b}$'),
                (g2_BZ, r'$g^2$')]):

            plot = storylines.Plot(
                style='APS',
                height=0,
                margin=margin,
                xyaxes=False,
                background='fig09%s.png' % labels[iq][i],
                )

            plot.width = (Margin / 2 + 5 * margin - plot.single) / 3

            plot.line(*zip(*elphmod.bravais.BZ()), thick=True)

            image = storylines.colorize(np.round(255 * data), cmap,
                minimum=0, maximum=255)

            storylines.save(plot.background, image)

            for contour in FS_kk:
                plot.line(*list(zip(*contour)))

            for contour in FS_kq:
                plot.line(*list(zip(*contour)),
                    color='white', densely_dotted=True)

            if i == 0:
                plot.left = Margin / 2
                plot.node(-kxmax, 0, r'$\mathbf q = \mathrm %s$' % Q[iq],
                    above=True, rotate=90)

            if iq == 0:
                plot.top = Margin / 2
                plot.title = title
            else:
                plot.bottom = Margin / 2

                if i == 0:
                    plot.node(kxmax, -kymax, '0', left=True,
                        yshift='-0.35cm', inner_sep=0, outer_sep=0)

                elif i == 1:
                    colorbar = storylines.colorize([[n / (BZ['points'] - 1.0)
                        for n in range(BZ['points'])]], cmap)

                    storylines.save('fig09.bar.png', colorbar)

                    plot.node(0, -kymax, r'\includegraphics'
                        '[width=<dx=%g>cm, height=2mm]{fig09.bar.png}'
                            % (2 * kxmax), yshift='-0.35cm')

                elif i == 2:
                    plot.node(-kxmax, -kymax, 'max', right=True,
                        yshift='-0.35cm', inner_sep=0, outer_sep=0)

            plot.node(-kxmax, kymax, '(%s)' % labels[iq][i], below_right=True,
                inner_sep=0, outer_sep=0)

            plot.save(plot.background.replace('.png', '.pdf'))

if comm.rank == 0:
    storylines.combine('fig09.pdf',
        ['fig09%s' % a for abc in labels for a in abc], columns=3)
