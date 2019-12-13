#!/bin/bash
#clear

echo -e "\n#############################################################"
echo -e   "##### Weak scaling of Intel Quantum Simulator #############"
echo -e   "#############################################################"

# WEAK SCALING:
# Assume that, approximately, the parallel part scales linearly with the amount
# of resources and that the serial part does not increase with respect to the size
# of the problem. Gustafson’s law provides the formula for scaled speedup:
#     scaled speedup = s + p × N

##########################################################################

# The script below is run for a system that uses slurm.
# Options are: "pcl-clx", "single-node".
flag_for_slurm="single-node"

declare -a set_rank_values=(1 2 4 8 16)

##########################################################################

echo -e "\n -- Setting the parameters that stay unchanged -- \n"

num_threads_per_rank=2

out_directory="output/"
out_filename_root="basic_weak_scaling"
min_num_qubits=26
num_gates=1

exec_file="../build/bin/basic_code_for_scaling.exe"

##########################################################################

# If the script is launched in a slurm machine.
# Here we considered 8 nodes with 2 sockets each.
if [ $flag_for_slurm == "pcl-clx" ]; then
	num_ranks_per_node=2
	num_threads_per_rank=28
fi

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

# We need to define a function returning the log2 of an integer.
function log2 {
    local x=0
    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
        let x=$x+1
    done
    echo $x
}

file_list=()
numq_list=()
for num_ranks in "${set_rank_values[@]}"
do
	num_global_qubits=$(log2 $num_ranks)
	num_qubits=$((min_num_qubits + num_global_qubits))
	exec_args=" -nq "$num_qubits" -ng "$num_gates\
" -od "$out_directory" -of "$out_filename_root\
" -nt "$num_threads_per_rank
	if [ $flag_for_slurm == "single-node" ]; then
		num_ranks_per_node=$num_ranks
	fi
	mpiexec.hydra -n $num_ranks -ppn $num_ranks_per_node -genv I_MPI_DEBUG 4 -genv OMP_NUM_THREADS $num_threads_per_rank -genv KMP_AFFINITY granularity=fine $exec_file $exec_args
	filename=$out_directory$out_filename_root"_q"$num_qubits"_n"$num_ranks".txt"
	file_list+=($filename)
	numq_list+=($num_qubits)
done

##########################################################################

echo -e "\n -- Graph of gate-time per qubit (weak scaling) -- \n"

data_filename_1=${file_list[0]}
data_filename_2=${file_list[1]}
data_filename_3=${file_list[2]}
data_filename_4=${file_list[3]}
data_filename_5=${file_list[-1]}
plot_filename=${out_directory}weak_scaling.png

#"set style line 3 lt rgb '#0e1a8b' pt 2 ps 1 lt 1 lw 2;"\

gnuplot_instructions="set xlabel \"qubits\" font \",10\";"\
"set ylabel \"time (sec)\" font \",10\";"\
"set style line 1 lt rgb '#8b1a0e' lw 2 pt 1 ps 0.3;"\
"set style line 2 lt rgb '#5e9c36' lw 2 pt 2 ps 0.5;"\
"set style line 3 lt rgb '#0060ad' lw 2 pt 3 ps 0.6;"\
"set style line 4 lt rgb '#ff949E' lw 2 pt 4 ps 0.8;"\
"set style line 5 lt rgb '#674C87' lw 2 pt 5 ps 0.5;"\
"set style line 6 lt rgb '#777777' lw 2 pt 6 ps 0.8;"\
"set style line 7 lt rgb '#222222' lw 2 pt 7 ps 0.5;"\
"set style line 8 lt rgb '#FBDF61' lw 2 pt 8 ps 0.5;"\
"set style line 9 lt rgb '#F25900' lw 2 pt 12 ps 0.6;"\
"set style line 10 lt rgb '#014F4B' lw 2 pt 10 ps 0.5;"\
"set terminal png size 800,600; set output '$plot_filename';"\
"plot '$data_filename_1' using 1:2 title \"${set_rank_values[0]} ranks / ${numq_list[0]} qubits\" with linespoints ls 1,"\
"'$data_filename_2' using 1:2 title \"${set_rank_values[1]} ranks / ${numq_list[1]} qubits\" with linespoints ls 2,"\
"'$data_filename_3' using 1:2 title \"${set_rank_values[2]} ranks / ${numq_list[2]} qubits\" with linespoints ls 3,"\
"'$data_filename_4' using 1:2 title \"${set_rank_values[3]} ranks / ${numq_list[3]} qubits\" with linespoints ls 4,"\
"'$data_filename_5' using 1:2 title \"${set_rank_values[-1]} ranks / ${numq_list[-1]} qubits\" with linespoints ls 5;"

gnuplot --persist -e "$gnuplot_instructions"
display $plot_filename &

##########################################################################
