#!/bin/csh

set init = step5_charmm2omm
set cnt = 1

while ( ${cnt} <= 6 )
    @ pcnt = ${cnt} - 1
    set istep = step6.${cnt}_equilibration
    set pstep = step6.${pcnt}_equilibration

    if ( ${cnt} == 1 ) then
        python -u openmm_run.py -i ${istep}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str -orst ${istep}.rst -odcd ${istep}.dcd -hmr > ${istep}.out
    else
        python -u openmm_run.py -i ${istep}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${pstep}.rst -orst ${istep}.rst -odcd ${istep}.dcd -hmr > ${istep}.out
    endif
    @ cnt += 1
end

