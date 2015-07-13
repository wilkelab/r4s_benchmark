#!/bin/bash
sim_model_arr=(dN dN_dS ms_dS ms_no_dS)

for model in ${sim_model_arr[*]}
do
	rm ./src/job_script_${model}
	out=sim_${model}
	nodes=1216
	time="24:00:00"
	
	echo "#!/bin/bash" >> ./src/job_script_${model}
	echo "#SBATCH -J $out   # job name" >> ./src/job_script_${model}
	echo "#SBATCH -o ${out}.o%j       # output and error file name (%j expands to jobID)" >> ./src/job_script_${model}
	echo "#SBATCH -n $nodes              # total number of mpi tasks requested" >> ./src/job_script_${model}
	echo "#SBATCH -p normal     # queue (partition) -- normal, development, etc." >> ./src/job_script_${model}
	echo "#SBATCH -t $time       # run time (hh:mm:ss) - 1.5 hours" >> ./src/job_script_${model}
	echo "#SBATCH --mail-user=dariya.k.sydykova@gmail.com" >> ./src/job_script_${model}
	echo "#SBATCH --mail-type=begin  # email me when the job starts" >> ./src/job_script_${model}
	echo "#SBATCH --mail-type=end    # email me when the job finishes" >> ./src/job_script_${model}
	echo "#SBATCH -A A-bio7" >> ./src/job_script_${model}
	echo "module load python" >> ./src/job_script_${model}
	echo "module load R" >> ./src/job_script_${model}
	echo "ibrun ./src/run_${model}_sim_pipeline.sh" >> ./src/job_script_${model} 
done       