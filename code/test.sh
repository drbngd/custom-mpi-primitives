module load openmpi
make clean
make





#output_16=$(srun -n 16 ./primitives -g 100000000 0)
#echo "8"
#echo "$output_16"

#output_16=$(srun -n 4 ./primitives -h 100 0)
#output_16=$(srun -n 16 ./primitives -h 100000000 0)
#output_16=$(srun -n 16 ./primitives -g 32 0)
#output_16=$(srun -n 2 ./primitives -g 2 0)
#output_16=$(srun -n 16 ./primitives -g 16 0)
#output_16=$(srun -n 16 ./primitives -g 100000000 0)
#output_16=$(srun -n 16 ./primitives -g 100000000 0)
#output_16=$(srun -n 2 ./primitives -g 100000000 0)
#output_16=$(srun -n 4 ./primitives -g 100000000 0)
#output_16=$(srun -n 8 ./primitives -g 100000000 0)

output_16=$(srun -n 16 ./primitives -s 100000 0)



echo "$output_16"
#echo "$output_24"