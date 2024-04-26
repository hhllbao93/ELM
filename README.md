# Getting started with offline ELM land only (to be used with WRF)


## 1. Install ParallelIO(PIO) and ESMF, which is necessary to compile ELM with MPI


### Compile ParallelIO (PIO) library (v2.5.9 is used as an example)
Get the code
```
   git clone git@github.com:NCAR/ParallelIO.git
```
Copy this `pio-build-perlmutter.sh` script to the parent directory of the unzipped code
```
   cp /PATH/TO/ELM_CODE/tools/pio-build-perlmutter.sh ./  
```
Update the file paths of unzipped code by changing -DCMAKE_INSTALL_PREFIX=/PATH/TO/PIO/BUILD
(the -S should point to where the code is, and the -B should point to where you want it to build)
```
   vim pio-build-perlmutter.sh
```
Run pio-build-perlmutter.sh
```
   sh pio-build-perlmutter.sh   
```
After it builds successfully, change the directory to the code folder 
```
   make
   make install
```

 
### Compile ESMF: (v8.4.2 is used as an example)
Get code for ESMF
```
   git clone git@github.com:esmf-org/esmf.git
```
Copy `esmf-build-perlmutter.sh` to inside the ESMF_CODE  directory
```
   cp /PATH/TO/ELM_CODE/tools/esmf-build-perlmutter.sh` /PATH/TO/ESMF_CODE
```
Update the ESMF_DIR, PIO, ESMF_INSTALL_PREFIX, ESMF_PIO_INCLUDE, ESMF_PIO_LIBPATH in esmf-build-perlmutter.sh
``` 
   vim esmf-build-perlmutter.sh
```
With the last line saying “make”
```
   sh esmf-build-perlmutter.sh
```
With the last line saying “make install”
```
   sh esmf-build-perlmutter.sh
```
Make ESMF accessible as a module
Copy the  `esmf-perlmutter-modulefile` into a directory named `esmf`
```
   cp esmf-perlmutter-modulefile ./esmf
```
Update the name of this file to `perlmutter-8.4.2`
```
   mv esmf-perlmutter-modulefile perlmutter-8.4.2
```
Update the paths in the file to point to where you installed PIO and ESMF
```
   vim perlmutter-8.4.2 
```

## 2. Build ELM and WRF

### Clone ELM and WRF code 
Clone the WRF repository and checkout develop branch:
```
    git clone https://github.com/wrf-model/WRF.git WRF-ELM
    cd WRF-ELM
    git checkout develop
```

Clone the ELM repository:
```
    git clone https://github.com/hhllbao93/ELM.git ELM
    cd ELM
```


### Build ELM and its dependencies
In your ELM code directory, build ELM and its dependencies. Currently, we only support building WRF-ELM on Perlmutter with gnu
```
    ./lilac/build_ctsm /PATH/TO/BUILD/ELM --machine perlmutter --compiler gnu
```
Changes need to be made in /PATH/TO/ELM_CODE/ccs_config/machines/config_machines.xml, especially
```
      <modules compiler="gnu">
        <command name="load">PrgEnv-gnu/8.3.3</command>
        <command name="load">gcc/11.2.0</command>
        <command name="load">craype/2.7.20</command>
        <command name="load">cpe/23.03</command>
        <command name="load">cray-mpich/8.1.25</command>
        <command name="load">cray-hdf5-parallel/1.12.2.3</command>
        <command name="load">cray-netcdf-hdf5parallel/4.9.0.3</command>
        <command name="load">cray-parallel-netcdf/1.12.3.3</command>
        <command name="load">cmake/3.22.0</command>
        <command name="load">libfabric</command>
        <command name="load">cray-libsci/23.02.1.1</command>
        <command name="use">/global/cfs/cdirs/m3878/hhllbao/esmf/modulefiles/</command>
        <command name="load">esmf/perlmutter-8.4.2</command>
      </modules>
```


### Building WRF with ELM
Load the same modules and set the same environments as used for ELM build by sourcing elm_build_environment.sh for Bash:
```
    source elm_build_environment.sh
```
Set makefile variables from ELM needed for the WRF build by setting the following environment. For example for Bash
```
    export WRF_ELM_MKFILE=/PATH/TO/BUILD/ELM/ctsm.mk
```
configure WRF code with Compiler choice: 34; Nesting option: 1
```
    ./configure
```
Update configure.wrf based on example /PATH/TO/WRF/configure.wrf.example
Compile the code using build_WRF.sh
```
    sh build_WRF.sh
```

## RUN WRF-ELM
Copy necessary namelist for ELM to WRF run directory

lnd_in: This is the main namelist input file for CTSM

lnd_modelio.nml: This sets CTSM’s PIO (parallel I/O library) configuration settings

lilac_in: This namelist controls the operation of LILAC

The surfacedata, mesh file, and domain file needs to be interpolated to WRF grids (e.g. geo_em.d01.nc)

Sample script to run WRF-ELM on Perlmutter
```
#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=128
###SBATCH --cpus-per-task=1
#SBATCH -C cpu
#SBATCH --account=m3878

##export OMP_PROC_BIND=spread
##export OMP_PLACES=threads
##export OMP_NUM_THREADS=8
module load cpu
module load cmake
module load PrgEnv-gnu
#module load gcc/11.2.0
module load cray-hdf5-parallel/1.12.2.3
module load cray-netcdf-hdf5parallel/4.9.0.3
module load cray-parallel-netcdf
module load cray-mpich/8.1.25
module load libfabric

###export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MPICH_DIR}/lib
export MPI_C_LIBRARIES=${MPICH_DIR}/lib/libmpich.a
export MPI_Fortran_LIBRARIES=${MPICH_DIR}/lib/libmpichf90.a
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/REDIST/cuda/11.5/targets/x86_64-linux/lib:$LD_LIBRARY_PATH

srun wrf.exe
```
