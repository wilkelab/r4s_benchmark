#!/bin/bash
sim_model_arr=(dN dN_dS ms_dS ms_no_dS)

for model in ${sim_model_arr[*]}
do
	rm ./src/job_script_${model}
	out=sim_${model}
	nodes=1216
	time="12:00:00"
	
	echo "#!/bin/bash" >> ./src/job_script_${model}
	echo "#$ -V	#Inherit the submission environment" >> ./src/job_script_${model}
	echo "#$ -cwd	# Start job in submission directory" >> ./src/job_script_${model}
	echo "#$ -N $out	# Job Name" >> ./src/job_script_${model}
	echo "#$ -j y	# Combine stderr and stdout" >> ./src/job_script_${model}
	echo "#$ -o $JOB_NAME.o$JOB_ID	# Name of the output file (eg. myMPI.oJobID)" >> ./src/job_script_${model}
	echo "#$ -pe 12way 96	# Requests 12 tasks/node, 24 cores total" >> ./src/job_script_${model}
	echo "#$ -q normal	# Queue name normal" >> ./src/job_script_${model}
	echo "#$ -l h_rt=15:00:00	# Run time (hh:mm:ss) - 1.5 hours" >> ./src/job_script_${model}
	echo "#$ -M dariya.k.sydykova@gmail.com	# Address for email notification" >> ./src/job_script_${model}
	echo "#$ -m be	# Email at Begin and End of job" >> ./src/job_script_${model}
 	echo "module load python" >> ./src/job_script_${model}
 	echo "module load R" >> ./src/job_script_${model}
	echo "set -x	# Echo commands, use set echo with csh" >> ./src/job_script_${model}
	echo "ibrun ./src/run_${model}_sim_pipeline.sh	# Run the MPI executable named a.out" >> ./src/job_script_${model}
	
	# echo "#!/bin/bash" >> ./src/job_script_${model}
# 	echo "#SBATCH -J $out   # job name" >> ./src/job_script_${model}
# 	echo "#SBATCH -o ${out}.o%j       # output and error file name (%j expands to jobID)" >> ./src/job_script_${model}
# 	echo "#SBATCH -n $nodes              # total number of mpi tasks requested" >> ./src/job_script_${model}
# 	echo "#SBATCH -p normal     # queue (partition) -- normal, development, etc." >> ./src/job_script_${model}
# 	echo "#SBATCH -t $time       # run time (hh:mm:ss) - 1.5 hours" >> ./src/job_script_${model}
# 	echo "#SBATCH --mail-user=dariya.k.sydykova@gmail.com" >> ./src/job_script_${model}
# 	echo "#SBATCH --mail-type=begin  # email me when the job starts" >> ./src/job_script_${model}
# 	echo "#SBATCH --mail-type=end    # email me when the job finishes" >> ./src/job_script_${model}
# 	echo "#SBATCH -A A-bio7" >> ./src/job_script_${model}
# 	echo "module load python" >> ./src/job_script_${model}
# 	echo "module load R" >> ./src/job_script_${model}
# 	echo "ibrun ./src/run_${model}_sim_pipeline.sh" >> ./src/job_script_${model} 
done       