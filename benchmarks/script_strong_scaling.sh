#!/bin/bash
#clear

echo -e "\n#############################################################"
echo -e   "##### Strong scaling of Intel Quantum Simulator #############"
echo -e   "#############################################################"

# STRONG SCALING:
# The speedup is limited by the fraction of the serial part of the code that is not
# amenable to parallelization. Amdahl’s law can be formulated as follows:
#     speedup = 1 / (s + p / N)
# where s is the proportion of execution time spent on the serial part, p is the
# proportion of execution time spent on the part that can be parallelized, and
# N is the number of processors.

##########################################################################

# The script below is run for a system that uses slurm.
# Options are: "pcl-clx", "single-node".
flag_for_slurm="single-node"

declare -a set_rank_values=(1 2 4 8 16)

##########################################################################

echo -e "\n -- Setting the parameters that stay unchanged -- \n"

num_threads_per_rank=2

out_directory="output/"
out_filename_root="basic_strong_scaling"
num_qubits=26
num_gates=1

##########################################################################

# If the script is launched in a slurm machine.
# Here we considered 8 nodes with 2 sockets each.
if [ $flag_for_slurm == "pcl-clx" ]; then
	num_ranks_per_node=2
	num_threads_per_rank=28
fi

##########################################################################

exec_args=" -nq "$num_qubits" -ng "$num_gates" -od "$out_directory" -of "$out_filename_root\
" -nt "$num_threads_per_rank

exec_file="./bin/basic_code_for_scaling.exe"

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

echo -e "\n -- Run the basic executable -- \n"

for num_ranks in "${set_rank_values[@]}"
do
	if [ $flag_for_slurm == "single-node" ]; then
		num_ranks_per_node=$num_ranks
	fi
	mpiexec.hydra -n $num_ranks -ppn $num_ranks_per_node -genv I_MPI_DEBUG 4 -genv OMP_NUM_THREADS $num_threads_per_rank -genv KMP_AFFINITY granularity=fine $exec_file $exec_args
done

##########################################################################

echo -e "\n -- Graph of gate-time per qubit at highest number of ranks -- \n"

num_ranks=${set_rank_values[-1]}
data_filename=$out_directory$out_filename_root"_q"$num_qubits"_n"$num_ranks".txt"
plot_filename=${out_directory}one_qubit_gates.png

gnuplot_instructions="set xlabel \"qubits\" font \",10\";"\
"set ylabel \"time (sec)\" font \",10\";"\
"set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2;"\
"set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2;"\
"set terminal png size 800,600; set output '$plot_filename';"\
"plot '$data_filename' using 1:2 title \"$num_ranks ranks\" with linespoints ls 1;"

gnuplot --persist -e "$gnuplot_instructions"

echo -e "\n [ plot in file "$plot_filename" ] \n"
display $plot_filename &

##########################################################################

echo -e "\n -- Graph of gate-time per qubit at highest number of ranks -- \n"

data_filename_1=$out_directory$out_filename_root"_first_q"$num_qubits".txt"
data_filename_2=$out_directory$out_filename_root"_last_q"$num_qubits".txt"
plot_filename=${out_directory}strong_scaling.png

gnuplot_instructions="set xlabel \"num ranks\" font \",10\";"\
"set ylabel \"time (sec)\" font \",10\";"\
"set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2;"\
"set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2;"\
"set terminal png size 800,600; set output '$plot_filename';"\
"plot '$data_filename_1' using 1:2 title \"first qubit\" with linespoints ls 1,"\
"'$data_filename_2' using 1:2 title \"last qubit\" with linespoints ls 2;"

gnuplot --persist -e "$gnuplot_instructions"

echo -e "\n [ plot in file "$plot_filename" ] \n"
display $plot_filename &

##########################################################################
