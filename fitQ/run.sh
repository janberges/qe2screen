#!/bin/bash

source elphmodenv

url=http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard # [1, 2]

for pp in Ta.upf S.upf
do
    test -e $pp || (wget $url/$pp.gz && gunzip $pp)
done

# [1] van Setten et al., Comput. Phys. Commun. 226, 39 (2018)
# [2] Hamann, Phys. Rev. B 88, 085117 (2013)

nk=2

mpirun pw.x -nk $nk < scf.in | tee scf.out
mpirun ph.x -nk $nk < ph.in | tee ph.out
ph2epw

mpirun pw.x -nk $nk < nscf.in | tee nscf.out
mpirun -n $nk epw.x -nk $nk < epw.in | tee epw.out

mpirun pw.x -nk $nk < scf.in | tee scf.out
mpirun ph.x -nk $nk < ph0.in | tee ph0.out

python3 fitQ.py
