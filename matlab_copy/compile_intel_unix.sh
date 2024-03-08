cd ./src
module load intel/2018.3
module load openmpi/3.1.4
module load make
module load netcdf-fortran-mpi/4.5.1
setenv MPI_DSM_DISTRIBUTE yes
export NETCDF=/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/netcdf-fortran-mpi/4.5.1
export NETCDFLIB=/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/netcdf-fortran-mpi/4.5.1/lib
make
