Ниже указаны команды, использовавшиеся для сборки и запуска программ на IBM Polus:

1) Компиляция кода с последовательным выполнением:
g++ main_seq_M_N.cpp -fopenmp -o main_seq_M_N.out,
где M и N - размерности по заданию
2) Выставление в очередь на исполнение:
bsub -J NxM -q normal -W 3:00 -o main_seq_N_M_result.out -e main_seq_N_M_result.err ./main_seq_N_M.out

3) Компиляция кода с OpenMP:
g++ main_par_M_N_num.cpp -fopenmp -o main_par_M_N_num.out,
где num - количество OpenMP нитей
4) Выставление в очередь на исполнение:
bsub -J NxMxnum -q normal -W 3:00 -o main_par_M_N_num_result.out -e main_par_M_N_num_result.err ./main_par_M_N_num.out

5) Компиляция кода с MPI:
mpicxx -o main_mpi main_mpi.cpp
6) Выставление в очередь на исполнение:
bsub -n 1 -o /polusfs/home_edu/edu-cmc-skmodel24-616/edu-cmc-skmodel24-616-06/ -e /polusfs/home_edu/edu-cmc-skmodel24-616/edu-cmc-skmodel24-616-06/ mpiexec ./main_mpi

7) Компиляция кода с MPI+OpemMP:
mpicxx -qsmp=omp -o main_hybrid_num_M_N main_hybrid.cpp
8) Выставление в очередь на исполнение:
bsub -n num -R "affinity[thread(num, same=core)*1]" -o /polusfs/home_edu/edu-cmc-skmodel24-616/edu-cmc-skmodel24-616-06/ -e /polusfs/home_edu/edu-cmc-skmodel24-616/edu-cmc-skmodel24-616-06/ mpiexec ./main_hybrid_num_M_N