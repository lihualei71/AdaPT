#!/bin/bash

params_filename='corr_params_simul2.txt'

while read times seed
do
    filename='../results/AdaPT_simul2_seed_"'$seed'".sh'
    echo $filename
    Rout_filename='../results/AdaPT_simul2_seed_"'$seed'".Rout'
    echo $Rout_filename
    touch $filename
    echo '#!/bin/bash' > $filename
    echo '#SBATCH --cpus-per-task 1' >> $filename
    echo 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK' >> $filename
    echo '' >> $filename
    echo 'export times="'$times'"' >> $filename
    echo 'export seed="'$seed'"' >> $filename    
    echo "R CMD BATCH --no-save simul2.R "$Rout_filename >> $filename
    chmod 755 $filename
    sbatch -p high $filename
done < $params_filename
