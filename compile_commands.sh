#/bin/bash

#compile and run init program
gfortran -O3 Labinit_mphy.f  -o Labinit_mphy.exe
./Labinit_mphy.exe

#compile 
gfortran -O3 -g -fbacktrace array.f90 const.f90 rand1.f90 rkqs.f90 fcnkb.f90 fcnkb0.f90 rkck.f90 ew.f90 odeint.f90 size_dis.f90 microphysics.f90 LabExp_mphy.f90 -o LabExp_mphy.exe

#spinup to steady state
nohup ./LabExp_mphy.exe default nml_mpy steady_state 001 &> steady_state.log &
mkdir -p input/steady_state

#copy steady state profile and input params
cp  output/steady_state/001/*   input/steady_state/
cp  nml_mphy input/steady_state/

#edit LabExppar.dat to enable microphysics and change simulation time (for 1 hour its ~ 0.52D-1 in non-dim time)

#run odt model with microphysics
nohup ./LabExp_mphy.exe steady_state nml_mphy test 001 &> 001.log &


#run multiple realization of a ODT experiment
./odt_multi_realz.sh steady_state nml_mphy test &> multi_realz.log &

