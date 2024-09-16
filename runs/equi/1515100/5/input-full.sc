## General parameters
tasks 16 
cluster sparky
lammps_exe /home/saturas/development/lammps-eam-20Aug13/lammps-20Aug13/src/lmp_openmpi

## Build parameters
seed 123456789
a 8.99180 
b 5.18985 
c 2.47976
nx 15 
ny 15 
nz 100 
crystallinity 0.6 
initial_chain_cut_ratio 0.24
removal_ratio 0.07  #(0.8045-0.7483)/0.8045 # crystalline and amorphous phase density at 460K
stem 2,0,1
bond_stiffness 113.746321 
bond_length 2.573622
block E

## MC parameters
# when resuming a run, make sure the following files
# {amorphous_atoms.txt, root_atoms.txt, movable_atoms, 
#  accepted_steps, accepted.lammps
#  step-######-I.log, step-######-II.log}
# are in the corresponding folders.
# total number of melt, anneal, and sample steps
steps 20000 30000 10000
# temperature parameters.
# melt_T is determined by eq (1) below,
# need to perform some test runs to get a reasonable scope of max_energy_drop
# for the system under investigation.
melt_T 350000.0
target_T 460.0
# factor for the melt stage,
# max_energy_drop == factor*kB*melt_T
# => melt_T = max_energy_drop/(factor*kB)   (1)
# due to the accept ratio is mainly determined by ratio == exp(-energy_drop/(kB*melt_T))
# the optimal vaule of factor is about 2.0,
# so the smallest ratio is exp(-2.0)==0.135 and not too many cycles will be wasted.
max_drop_factor 2.0
# defines start and target search radii
search_radius 9.0 6.0
# threshold to check for stretched loops.
length_factor 1.2

## Parameters for saving accepted states
save_frequency 50

## MD parameters
md_dt 5.0 
md_steps 1000   # 1000 steps is sufficient for a 15x15x80
                # system with crystallinity 0.5
bead_mass 28.0538 
crystal_buffer_depth 2
cutoff_distance 16.0

## Post processing parameters
#e2e
e2e_start_index 0
#pp_steps 0 30 2   # do every other step
pp_steps all
slice_dz_factor 1.970
n_neighbor 18
