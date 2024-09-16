#!/usr/bin/env python
import numpy
import os

def equilibration_input(md_dt, velocity_seed, run_steps, current_stage, MD_id):
    """ Writes lammps input script for equilibration step """
    # Generate keys for script_I or script_II.
    if current_stage == 'melt':
        nve_setting  = 'nve/limit 4.0'
        temp_setting = 'fix 3 movable_atoms temp/rescale 1 460 460 10 1.0'
        temp_unfix   = 'unfix 3'
    else:
        nve_setting  = 'nve'
        temp_setting = ''
        temp_unfix   = ''
    # Generate input script for LAMMPS run.
    if MD_id == 'I':
        script = script_I.format(**locals())
        script = script.format(nve_setting, temp_setting)
        data_file = 'accepted.lammps'
    elif MD_id == 'II':
        script = script_II.format(**locals())
        script = script.format(nve_setting, temp_setting, temp_unfix)
        data_file = 'trial_EB_start.lammps'
    else:
        print 'Wrong MD_id provided. I or II.'
    in_script = main_script.format(data_file, script)
    open('in.lammps', 'w').write(in_script)

def first_step_energy(step, MD_id):
    if MD_id == 'I':
        log = os.path.join('data', 'log', 'step-{:06d}-I.log'.format(step))
    else:
        print 'Wrong paramwter (I).'
    try:
        # Reads LAMMPS log file and returns the total energy of the first time step.
        fid = open(log,'r')
        while True:
            line = fid.readline()
            if (line == ''): break
            if line.startswith('Step'):
                variables = line.split()
                ij = variables.index('TotEng')
                data = fid.readline().split()
                H = float(data[ij])
                break
    except:
        return None
    return H

def last_step_energy(step, MD_id):
    if MD_id == 'I':
        log = os.path.join('data', 'log', 'step-{:06d}-I.log'.format(step))
    elif MD_id == 'II':
        log = os.path.join('data', 'log', 'step-{:06d}-II.log'.format(step))
    else:
        print 'Wrong paramwter (I or II).'
    # Reads LAMMPS log file and returns the total energy of the last time step.
    steps_completed = 0
    try:
        with open(log, 'r') as f:
            variables = None
            for line in f:
                if line.startswith('Warning') or line.startswith('WARNING'):
                    continue
                elif line.startswith('Step'):
                    variables = line.split()
                    ij = variables.index('TotEng')
                    #ij = [variables.index(s) for s in ['PotEng','KinEng']]
                elif line.startswith('Loop'):
                    variables = None
                    steps_completed += 1
                elif variables:
                    data = line.split()
                    if len(data)==len(variables):
                        H = float(data[ij])
                        #H = sum(float(data[i]) for i in ij)
    except:
        return None
    # If less than two runs completed, the run has failed.
    if steps_completed < 2:
        return None
    # Return average of last half of H.
    return H

main_script = """
units           real
atom_style      angle
boundary        p p p 
bond_style      harmonic
pair_style      table linear 1001 
angle_style     table linear 1001
read_data       {0}

bond_coeff      1 113.746321 2.573622
pair_coeff      1 1 ../pair.table.EE.44 EE
angle_coeff     1 ../angle.table.EEE.44 EEE
neighbor        3.0 bin
special_bonds   lj/coul 0.0 0.0 1.0 

# Freeze top crystalline and bottom crystalline regions.
include         movable_atoms
group           frozen_atoms subtract all movable_atoms

# Accelerate computation by turning off the pairwise interactions among frozen atoms.
neigh_modify    delay 0 every 1 check yes exclude group frozen_atoms frozen_atoms
fix             1 frozen_atoms setforce 0.0 0.0 0.0

{1}
"""

script_I = """
# Equilibration step - NVT
velocity        movable_atoms create 460 {velocity_seed}
fix             2 movable_atoms nvt temp 460 460 100
timestep        {md_dt}
thermo_style    custom step temp press etotal ke pe epair ebond eangle enthalpy vol
thermo          50
run             {run_steps}
unfix 2

# Equilibration step - NVE
fix             2 movable_atoms {nve_setting}
{temp_setting}
thermo_style    custom step temp press etotal ke pe epair ebond eangle vol
thermo          20
timestep        1
run             500
write_data      trial.lammps
"""

script_II = """
# Equilibration step - NVE
fix             2 movable_atoms {nve_setting}
{temp_setting}
thermo_style    custom step temp press etotal ke pe epair ebond eangle vol
thermo          20
timestep        1
reset_timestep  0
run             500
unfix 2
{temp_unfix}

# Equilibration step - NVT
fix             2 movable_atoms nvt temp 460 460 100
timestep        {md_dt}
thermo_style    custom step temp press etotal ke pe epair ebond eangle enthalpy vol
thermo          50
run             {run_steps}
write_data      trial.lammps
"""

if __name__ == '__main__':
    import sys
    print mean_energy(sys.argv[1], 1)
