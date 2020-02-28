#!/bin/csh

set cnt = 1
set cntmax = 7

set input = step7_production

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = step7_${cnt}
    set pstep = step7_${pcnt}

    if ( ${cnt} == 1 ) set pstep = step6.6_equilibration

    python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${pstep}.rst -orst ${istep}.rst -odcd ${istep}.dcd -hmr > ${istep}.out
    @ cnt += 1
end

