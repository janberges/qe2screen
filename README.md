# To screen, or not to screen

This repository contains the source code and data associated with the paper:

- Jan Berges, Nina Girotto, Tim Wehling, Nicola Marzari, and Samuel Poncé,
  *Phonon self-energy corrections: To screen, or not to screen*,
  [arXiv:2212.11806](https://arxiv.org/abs/2212.11806) (2022).

## Quantum ESPRESSO patch

The data shown in our paper have been calculated with a modified version of the
PHonon and EPW codes of Quantum ESPRESSO. In particular, we have introduced the
inputs and outputs listed below.

### Installation

Attached we provide a patch with the relevant modifications that can be applied
to version 6.8 of Quantum ESPRESSO. This can either be done with Git:

    git clone https://gitlab.com/QEF/q-e.git
    cd q-e
    git checkout qe-6.8
    git apply path/to/qe2screen.patch

Or without Git:

    wget https://github.com/QEF/q-e/releases/download/qe-6.8/qe-6.8-ReleasePack.tgz
    tar -xzf qe-6.8-ReleasePack.tgz
    cd qe-6.8
    patch -p1 < path/to/qe2screen.patch

Afterward, Quantum ESPRESSO can be installed as usual.

### Changes in PHonon

We have introduced the following inputs in `ph.x`:

- `cdfpt_min`: Lower bound of cDFPT active energy window in eV
- `cdfpt_max`: Upper bound of cDFPT active energy window in eV
    - If not set, a standard DFPT calculation is done
    - Default `dis_froz_min` and `dis_froz_max` in EPW
    - The energy is *not* measured from the Fermi level
- `bare`: Suppress electronic response to atomic displacements? [`false`]
- `ibndprt`: Band index for which k = 0 vertex shall be printed

A cDFPT calculation creates an additional output file:

- `outdir/_ph0/prefix.phsave/cdfpt_subspace.xml`: Information passed on to EPW

### Changes in EPW

We have introduced the following inputs in `epw.x`:

- `cdfpt_dir`: Equivalent of `dvscf_dir` for cDFPT data
- `xdfpt_dir`: Equivalent of `dvscf_dir` for extra DFPT data (low smearing)
- `unscreen_fine`: Perform unscreening on fine mesh? [`false`]
- `geff`: Drop negative eigenvalue of outer product of vertices? [`false`]
- `bare`: Is the data in `cdfpt_dir` bare (relevant for unscreening)? [`false`]
- `types`: Smearing types corresponding to `temps` [`-99`]
    - `0`: Gaussian
    - `1`: Methfessel-Paxton
    - `-1`: Marzari-Vanderbilt
    - `-99`: Fermi-Dirac
- `temp_inf`: High smearing at which phonons are interpolated in K
- `type_inf`: Smearing type corresponding to `temp_inf` [`-99`]
- `corr`: Method to correct screened vertex [`0`]
    - `0`: No correction
    - `1`: cRPA correction I
    - `2`: Linear correction II
    - Negative band index: Correction II applied to single band
- `phlabel`: Filename stem for band plots [`phband`]
- `prtgkkc`: Equivalent of `prtgkk` for cDFPT data [`false`]
- `prtgkkx`: Equivalent of `prtgkk` for extra DFPT data [`false`]
- `lpolarc`: Equivalent of `lpolar` for cDFPT data [`true`]
- `lpolarx`: Equivalent of `lpolar` for extra DFPT data [`false`]
- `L`: Range-separation parameter [`0`]
- `perp`: Take out-of-plane polarizability into account (for `L > 0`)? [`false`]

Some inputs are used differently in the cDFPT case:

- `asr_typ`: The new option `'none'` disables the ASR correction [`'simple'`]
- `temps`: Low smearings for which the phonon frequencies are estimated in K
- `degaussw`: Self-energy smearing for phonon spectral function in eV [`0.025`]
- `degaussq`: Overall smearing for phonon spectral function in meV [`0.05`]

EPW also watches for new optional input files:

- `dipole.fmt`: Alternative DFPT Born effective charges
- `dipolec.fmt`: Alternative cDFPT Born effective charges
- `dipolex.fmt`: Alternative extra DFPT Born effective charges
- `quadrupolec.fmt`: Equivalent of `quadrupole.fmt` for cDFPT data
- `quadrupolex.fmt`: Equivalent of `quadrupole.fmt` for extra DFPT data

In a cDFPT calculation, `band_plot` produces several output files:

- `phlabel.freq`: DFPT phonon frequencies
    - Custom filename stem defined via the input `phlabel`
- `phlabelc.freq`: cDFPT phonon frequencies
- `phlabelx.freq`: Extra DFPT phonon frequencies
- `phlabelr.freq`: Renormalized DFPT phonon frequencies
    - Phonon frequencies at `temp_inf` if given
    - Otherwise unscreened DFPT phonon frequencies
- `phlabelcr.freq`: Renormalized cDFPT phonon frequencies
    - Phonon frequencies at `temp_inf` if given
    - Otherwise equal to DFPT phonon frequencies (sanity check)
- `phlabelrn.freq`: DFPT-DFPT phonon frequencies for `n`th smearing
- `phlabelcrn.freq`: cDFPT-DFPT phonon frequencies for `n`th smearing

`specfun_ph` additionally produces:

- `phlabelrwn.freq`: DFPT-DFPT phonon spectral function for `n`th smearing
- `phlabelcrwn.freq`: cDFPT-DFPT phonon spectral function for `n`th smearing

`eig_plot` additionally produces:

- `phlabelren.freq`: DFPT-DFPT phonon eigenvectors for `n`th smearing
- `phlabelcren.freq`: cDFPT-DFPT phonon eigenvectors for `n`th smearing

## Figure scripts and data

The directories `figXY` contain the Python scripts and data necessary to create
all figures shown in our paper. This can be done in a virtual environment:

    python3 -m venv venv
    source venv/bin/activate
    python3 -m pip install -r requirements.txt

A LaTeX installation, preferably TeX Live, is required to typeset the figures.

## Optimization of quadrupole tensors

The directory `fitQ` contains a minimal working example of the optimization of
quadrupole tensors as done in our paper. If both the patched version of Quantum
ESPRESSO and the above Python environment are installed, simply do `./run.sh`.

More precisely, we first calculate (bare) dynamical matrices and electron-phonon
matrix elements both on a coarse q mesh and for selected q points along a path.
Then we Fourier-interpolate the former, minimizing deviations from the latter.
Here, the free parameters are the independent elements of the quadrupole tensors
Q and the range-separation parameter L entering the formulas for the long-range
components that are subtracted and added before and after interpolation. First,
we optimize L for Q = 0, minimizing the short-range part of the force constants.
Second, we optimize Q for constant L. A simultaneous optimization would also be
possible, but it is only important that L is in the correct range.

Note that the results obtained in this example are not converged. The parameters
have been chosen such that the calculations can be done on a personal computer.

## Quantum ESPRESSO input files

The directory `input` contains the Quantum ESPRESSO input files we have used.
Rather than a separate file for each set of parameters considered, we provide
one generic example and indicate possible variants by comments.

## Licence

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation (<https://www.gnu.org/licenses>), either version 2 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

Copyright (C) 2023 J. Berges, N. Girotto, T. Wehling, N. Marzari, S. Poncé
