# Getting started with offline ELM land only (to be used with WRF)

## 1. To compile ELM with MPI, ParallelIO(PIO) and ESMF is needed

### Compile ParallelIO (PIO) library (v2.5.9 is used as an example)
Get the code
```
   git clone git@github.com:NCAR/ParallelIO.git
```
Copy this `pio-build-perlmutter.sh` script to the parent directory of the unzipped code
```
   cp /PATH/TO/ELM_CODE/tools/pio-build-perlmutter.sh ./  
```
Update the file paths of unzipped code 
```
   vim pio-build-perlmutter.sh
   #by changing -DCMAKE_INSTALL_PREFIX=/global/cfs/cdirs/m3878/hhllbao/pio/perlmutter-2.5.9
   # the -S should point to where the code is, and the -B should point to where you want it to build
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
```
   git clone git@github.com:esmf-org/esmf.git  
   cp /PATH/TO/ELM_CODE/tools/esmf-build-perlmutter.sh` /PATH/TO/ESMF_CODE # Copy this `esmf-build-perlmutter.sh` to inside the unzipped directory
   vim esmf-build-perlmutter.sh  # Update the file paths in esmf-build-perlmutter.sh
   sh esmf-build-perlmutter.sh   # with the last line saying “make”
   sh esmf-build-perlmutter.sh   # with the last line saying “make install”
   # make ESMF accessible as a module:
   cp esmf-perlmutter-modulefile ./esmf # Copy this `esmf-perlmutter-modulefile` into a directory named `esmf`
   mv esmf-perlmutter-modulefile perlmutter-8.4.2 # Update the name of this file to `perlmutter-8.4.2`
   vim perlmutter-8.4.2 # Update the paths in the file to point to where you installed PIO and ESMF
```

## 2. Clone ELM and WRF code 
### Clone the WRF repository and checkout develop branch:
```
    git clone https://github.com/wrf-model/WRF.git WRF-ELM
    cd WRF-ELM
    git checkout develop
```

### Clone the ELM repository:
```
    git clone https://github.com/hhllbao93/E3SM.git ELM
    cd ELM
    ./manage_externals/checkout_externals 
```

## 3. Build ELM and its dependencies
### In your ELM directory, build ELM and its dependencies. Currently we only support build WRF-ELM on Perlmutter with gnu
```
    ./lilac/build_ctsm /PATH/TO/ELM/BUILD --machine perlmutter --compiler gnu
```

## 4. Building WRF with ELM
### Load the same modules and set the same environments as used for ELM build by sourcing elm_build_environment.sh for Bash:
```
    source elm_build_dir/elm_build_environment.sh
```
### Set makefile variables from ELM needed for the WRF build by setting the following environment. For example for Bash
```
    export WRF_ELM_MKFILE=/glade/scratch/$USER/WRF-ELM/ELM/elm_build_dir/bld/elm.mk
```
### Compile the code using build_WRF.sh
```
    sh build_WRF.sh
```
