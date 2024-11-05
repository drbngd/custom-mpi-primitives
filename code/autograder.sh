#!/bin/bash

#SBATCH --job-name=mpi_results
#SBATCH --output=mpi_results.out
#SBATCH --error=mpi_results.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:05:30

# Maximum time allowed for the script (in seconds)
MAX_EXECUTION_TIME=300

# Run the tests within the time limit
timeout $MAX_EXECUTION_TIME bash -c '

lscpu
echo "-------------------------------------"

# Run the tests
module load openmpi

make clean
make

# Correctness tests
echo " "
echo "=================== CORRECTNESS TESTS ==================="
echo " "
echo "##Scatter"
output_16=$(srun -n 16 ./primitives -s 100000000 5)
output_24=$(srun -n 24 ./primitives -s 150000000 5)
echo "$output_16"
echo "----------------------------------"
echo "$output_24"
if [[ $output_16 != *"Implementation is correct!"* ]] || [[ $output_24 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
fi
echo "=================================="
echo " "
echo "##AllGather"
output_16=$(srun -n 16 ./primitives -g 100000000 0)
output_16=$(srun -n 16 ./primitives -g  0)
output_24=$(srun -n 24 ./primitives -g 150000000 0)
echo "$output_16"
echo "----------------------------------"
echo "$output_24"
if [[ $output_16 != *"Implementation is correct!"* ]] || [[ $output_24 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
fi
echo "=================================="
echo " "
echo "##AllReduce"
output_16=$(srun -n 16 ./primitives -r 100000000 0)
output_24=$(srun -n 24 ./primitives -r 150000000 0)
echo "$output_16"
echo "----------------------------------"
echo "$output_24"
if [[ $output_16 != *"Implementation is correct!"* ]] || [[ $output_24 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
fi
echo "=================================="
echo " "
echo "##AlltoAll_Arbitrary"
output_16=$(srun -n 16 ./primitives -a 100000000 0)
output_24=$(srun -n 24 ./primitives -a 150000000 0)
echo "$output_16"
echo "----------------------------------"
echo "$output_24"
if [[ $output_16 != *"Implementation is correct!"* ]] || [[ $output_24 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
fi
echo "=================================="
echo " "
echo "##AlltoAll_Hypercube"
output_16=$(srun -n 16 ./primitives -h 100000000 0)
echo "$output_16"
if [[ $output_16 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
fi
echo "=================================="
echo " "


echo "==================== RUNTIME TESTS ====================="
echo " "

multiplier=25
echo "##Scatter"
output_16=$(srun -n 16 ./primitives -s 100000000 0)
echo "$output_16"
if [[ $output_16 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
else
    custom_16=$(echo "$output_16" | grep "Custom_" | awk "{print \$3}")
    mpi_16=$(echo "$output_16" | grep "MPI_" | awk "{print \$3}")

    if (( $(echo "$custom_16 < $multiplier * $mpi_16" | bc -l) )); then
        echo "Runtime test passed!"
    else
    echo "Runtime test failed"
    fi
fi

echo "=================================="
echo " "
echo "##AllGather"
multiplier=3
output_16=$(srun -n 16 ./primitives -g 100000000 0)
echo "$output_16"
if [[ $output_16 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
else
	custom_16=$(echo "$output_16" | grep "Custom_" | awk "{print \$3}")
	mpi_16=$(echo "$output_16" | grep "MPI_" | awk "{print \$3}")

	if (( $(echo "$custom_16 < $multiplier * $mpi_16" | bc -l) )); then
    		echo "Runtime test passed!"
	else
    	echo "Runtime test failed"
	fi
fi

echo "=================================="
echo " "
echo "##AllReduce"
multiplier=10
output_16=$(srun -n 16 ./primitives -r 100000000 0)
echo "$output_16"
if [[ $output_16 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
else
	custom_16=$(echo "$output_16" | grep "Custom_" | awk "{print \$3}")
	mpi_16=$(echo "$output_16" | grep "MPI_" | awk "{print \$3}")

	if (( $(echo "$custom_16 < $multiplier * $mpi_16" | bc -l) )); then
    		echo "Runtime test passed!"
	else
    		echo "Runtime test failed"
	fi
fi

echo "=================================="
echo " "
echo "##AlltoAll_Arbitrary"
multiplier=2
output_16=$(srun -n 16 ./primitives -a 100000000 0)
echo "$output_16"
if [[ $output_16 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
else
	custom_16=$(echo "$output_16" | grep "Custom_" | awk "{print \$3}")
	mpi_16=$(echo "$output_16" | grep "MPI_" | awk "{print \$3}")

	if (( $(echo "$custom_16 < $multiplier * $mpi_16" | bc -l) )); then
    		echo "Runtime test passed!"
	else
    		echo "Runtime test failed"
	fi
fi

echo "=================================="
echo " "
echo "##AlltoAll_Hypercube"
multiplier=10
output_16=$(srun -n 16 ./primitives -h 100000000 0)
echo "$output_16"
if [[ $output_16 != *"Implementation is correct!"* ]]; then
    echo "ERROR"
else
	custom_16=$(echo "$output_16" | grep "Custom_" | awk "{print \$3}")
	mpi_16=$(echo "$output_16" | grep "MPI_" | awk "{print \$3}")

	if (( $(echo "$custom_16 < $multiplier * $mpi_16" | bc -l) )); then
    		echo "Runtime test passed!"
	else
    		echo "Runtime test failed"
	fi
fi

echo " "

'

if [[ $? -eq 124 ]]; then
    echo "Error: Script execution exceeded $MAX_EXECUTION_TIME seconds."
    exit 1
fi
