#!/bin/bash
#clear

echo -e "\n#############################################################"
echo -e   "##### test of computational time using slurp ################"
echo -e   "#############################################################"

##########################################################################

echo -e "\n -- Setting the simulation parameters -- \n"

num_ranks=16
num_threads_per_rank=1

out_directory="output/"
out_filename_root="basic_strong_scaling"
num_qubits=28
num_gates=1

exec_args=" -nq "$num_qubits" -ng "$num_gates" -od "$out_directory" -of "$out_filename_root\
" -nt "$num_threads_per_rank

##########################################################################

exec_name="basic_scaling"
exec_file="../build/bin/basic_code_for_scaling.exe"

##########################################################################

echo -e "\n -- Create the output folder if not present -- \n"

if [ ! -d $out_directory ]; then
	mkdir $out_directory
else
	# Eliminate the summary files if they exist.
	filename=$out_directory$out_filename_root"_first_q"$num_qubits".txt"
        echo "% Time cost (in sec) for 1-qubit gate on first qubit vs num_ranks" > $filename 
	filename=$out_directory$out_filename_root"_last_q"$num_qubits".txt"
        echo "% Time cost (in sec) for 1-qubit gate on last qubit vs num_ranks" > $filename 
fi

##########################################################################

echo -e " -- Creating the slurm job file -- \n"

job_file="job_"$exec_name".slurm"
job_content=\
"#!/bin/bash"$'\n\n'\
"#SBATCH -o "$out_directory$exec_name".o%j"$'\n'\
"#SBATCH -e "$out_directory$exec_name".err%j"$'\n'\
"#SBATCH -D ./"$'\n'\
"#SBATCH --get-user-env"$'\n'\
"#SBATCH --partition=clx"$'\n'\
"#SBATCH --time=01:00:00"$'\n'\
"#SBATCH --ntasks="$num_ranks$'\n'\
"#SBATCH --cpus-per-task="$num_threads_per_rank$'\n'\
"#SBATCH -J "$exec_name$'\n'\
"#SBATCH --mail-user=gian.giacomo.guerreschi@intel.com"$'\n'\
"#SBATCH --mail-type=begin --mail-type=end"$'\n\n'\
"export INTEL_LICENSE_FILE=/swtools/intel/licenses/"$'\n'\
"export I_MPI_HYDRA_BOOTSTRAP=slurm"$'\n\n'\
"mpiexec.hydra -genv I_MPI_DEBUG 4 -genv OMP_NUM_THREADS "$num_threads_per_rank\
" -genv KMP_AFFINITY granularity=thread,1,0,verbose -ppn $num_ranks ./"$exec_file" "$exec_args

echo "$job_content"  > $job_file

##########################################################################

echo -e "\n -- Running the simulations in a single job -- \n"
sbatch $job_file

##########################################################################
