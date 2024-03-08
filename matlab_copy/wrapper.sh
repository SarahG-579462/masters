module load matlab/2020a
echo "$SLURM_ARRAY_TASK_ID"
matlab -nosplash -nodesktop -r "run full_analysis_$SLURM_ARRAY_TASK_ID.m"
