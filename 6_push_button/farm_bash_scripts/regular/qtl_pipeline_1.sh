#! /bin/bash -l
#SBATCH -D /home/jkta/projects/col_sha/complete_pipeline
#SBATCH -o /home/jkta/projects/col_sha/complete_pipeline/logs/output/qtl_I_%j.txt
#SBATCH -e /home/jkta/projects/col_sha/complete_pipeline/logs/error/qtl_I_%j.txt
#SBATCH -J qtl_I

module load R

R CMD BATCH --no-save --no-restore "--args ${SLURM_ARRAY_TASK_ID}" qtl_pipeline_push_button_part_I.R /home/jkta/projects/col_sha/complete_pipeline/logs/r_out/qtl_I_${SLURM_ARRAY_TASK_ID}.txt

echo "Done"
