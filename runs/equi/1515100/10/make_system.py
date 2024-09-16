#!/usr/bin/env python
""" Runs a MC algorithm to generate an layered semicrystalline PE domain.
    NOTE: you must create a symlink to LAMMPS in this folder prior to running.

    Any single cycle is indicated as:
        the first MD part (I, NVT+NVE)
        end-bridging (EB)
        the second MD part (II, NVE+NVT).
    
    Three files are iteratively updated in the run folder:
        trial_EB_start.lammps
        trial.lammps
        accepted.lammps

    Main procedure:

    1) EB -> initialize model to generate trial_EB_start.lammps
       then copy it to step_0.lammps in data folder.
    2) Generate in.lammps for (II), then un LAMMPS to to read
       trial_EB_start.lammps and generate trial.lammps,
       then copy trial.lammps to accepted.lammps.
    3) For step in steps:
        3.1) Determine which stage (melt,anneal,sample) the current step is in.
        3.2) Generate in.lammps for (I), run LAMMPS to read accepted.lammps and
             update trial.lammps (overwrite).
        3.3) EB -> read trial.lammps to update trial_EB_start.lammps 
             (overwrite).
        3.4) Generate in.lammps for (II), run LAMMPS to read
             trial_EB_start.lammps and update trial.lammps (overwrite).
        3.5) Perform MC acceptance test for trial.lammps:
             If accepted:
                Copy trial.lammps to accepted.lammps (overwrite).
                Copy trial.lammps to step_###.lammps in data folder
                if specified.
             If rejected:
                Do nothing except copy trial_EB_start.lammps to the
                low_energy_state folder when condition is satisfied.
"""
import input_file
import subprocess
import math
import numpy
import os
import shutil
import glob
import lammps
import sys
import re
import metropolis
import time

molbuild_exe = './molbuilder'
# Boltzmann constant in kcal/mol/K
kB = 1.9872156e-3    

def main():
    """ Parses input file and runs either the builder or post processing. """
    try:
        mode  = sys.argv[1]
        ifile = sys.argv[-1]
    except IndexError:
        print './make_system (build) [input file]'
        print './make_system (resume) [input file]'
        print './make_system (post) [input file]'
        print './make_system (slip) [input file]'
        exit(1)
    opt = input_file.InputFile(ifile)

    # For now, don't let user change stem option.
    if not opt.stem == '2,0,1':
        print 'Stem vector cannot be changed yet.'
        exit(1)

    if mode == 'build':
        run_mc(opt, 0)
    elif mode == 'resume':
        step_begin = read_begin_step()
        run_mc(opt, step_begin)
    elif mode == 'post': 
        run_mc_post(opt)
    elif mode == 'slip': 
        run_mc_slip(opt)
    else: 
        print 'Unknown mode', mode

def run_mc(opt, step_begin):
    """ Runs the Monte Carlo algorithm with alternating calls to lammps and
    the MC c++ code """
    mc = metropolis.MetropolisAlgorithm(opt)
    if step_begin > 0:
        mc.lookup_last_accepted_step(step_begin)
        print 'Resuming from a previous run.'

    step_end = sum(int(s) for s in opt.steps.split())

    # Perform the build stage to generate data files for later use.
    print 'Running MC step {} to {}'.format(step_begin, step_end-1)
    if step_begin == 0:
        init_run_directory(opt)
        print 'Initilializing system'
        execute(molbuild_exe + ' build 0 ' + opt.input_file + ' melt',
                quit_on_error=True)
        # The first state is always accepted.
        subprocess.call(['cp', 'trial_EB_start.lammps',
                         os.path.join('data', 'step_0.lammps')])

    # Initialize parameters for generating lammps script.
    if opt.cluster == 'sparky':
        lmp_cmd = 'mpiexec {lammps_exe} -in in.lammps -sc none'
    else:
        lmp_cmd = 'mpiexec -n {tasks} {lammps_exe} -in in.lammps -sc none'
    lmp_cmd = lmp_cmd.format(**opt.options)

    # Generate lammps script for the initial (II) run.
    velocity_seed = int(opt.seed)
    md_dt = 1.0 # Use small dt for the intial system in NVT dynamics.
    lammps.equilibration_input(md_dt, velocity_seed, opt.md_steps, 'melt', 'II')

    # Perform hybrid Monte Carlo.
    counter = 0
    max_energy_drop = float(opt.max_drop_factor)*kB*mc.melt_temperature
    for step in range(step_begin, step_end):
        # Determine current stage.
        if step <= int(opt.steps.split()[0]):
            current_stage = 'melt'
        else:
            current_stage = 'non-melt'

        # Perform a single HMC cycle.
        if step > 0:
            # Generate lammps script for (I) and run LAMMPS to update trial.lammps.
            print '\n' + 40*'-'
            print 'Running MC step: {}'.format(step)
            velocity_seed = int(opt.seed) + step
            lammps.equilibration_input(opt.md_dt, velocity_seed, opt.md_steps,
                                       current_stage, 'I')
            print 'Calling LAMMPS for (I) ......'
            log = os.path.join('data', 'log', 'step-{:06d}-I.log'.format(step))
            execute(lmp_cmd + ' -l ' + log)
            # Test if the LAMMPS run for (I) is complete.
            if lammps.last_step_energy(step,'I') == None:
                print 5*'*', 'Rejected MC move during (I)', 5*'*'
                continue
            # Perform end-bridging to update trial_EB_start.lammps.
            Ctime_start = time.time()
            execute(molbuild_exe + ' build {} {} {}'.format(step,
                                                opt.input_file, current_stage),
                    quit_on_error=True)
            Ctime_end = time.time()
            # Generate LAMMPS script for (II).
            lammps.equilibration_input(opt.md_dt, velocity_seed, opt.md_steps,
                                       current_stage, 'II')

        # Run LAMMPS for (II) to update trial.lammps.
        print 'Calling LAMMPS for (II) ......'
        log = os.path.join('data', 'log', 'step-{:06d}-II.log'.format(step))
        LMPtime_start = time.time()
        execute(lmp_cmd + ' -l ' + log)
        LMPtime_end = time.time()
        if step > 0:
            print 'End-bridging takes', Ctime_end-Ctime_start, 's'
            print 'MD simulation for (II) takes', LMPtime_end-LMPtime_start, 's'

        # Perform MC acceptance test and postprocessing of LAMMPS files.
        E_start = lammps.first_step_energy(step,'I')
        E_trial = lammps.last_step_energy(step,'II')
        if mc.is_move_accepted(step, E_start, E_trial, max_energy_drop):
            print 5*'*', 'Accepted MC move', 5*'*'
            save_state(opt, step, counter)
            counter += 1
            mc.last_accepted_step = step
            open('accepted_steps', 'a').write('{}\n'.format(step))
        else:
            print 5*'*', 'Rejected MC move', 5*'*'
            if (not E_trial==None) and E_trial-E_start<-max_energy_drop:
                subprocess.call(['cp', 'trial_EB_start.lammps',
                                 'low_energy_state/step_{}_EB_start.lammps'.format(step)])
    # Clean files after cycle loop.
    subprocess.call('rm -f num_candidates trial.lammps trial_EB_start.lammps'.split())

def run_mc_post(opt):
    ''' Postprocess LAMMPS files. '''
    if not os.path.exists('postprocessing'): 
        os.mkdir('postprocessing')
    all_steps = glob.glob(os.path.join('data', '*.lammps'))
    all_steps = sorted([int(re.search(r'_(\d+)',s).group(1)) for s in all_steps])
    if opt.is_option_set('e2e'):
        # end-to-end calculation is fast, so do all steps >= e2e_start_index.
        start_index = int(opt.e2e_start_index)
        steps = [i for i in all_steps if i >= start_index]
        if start_index == 0:
            subprocess.call(['rm','-f', os.path.join('postprocessing', 'e2e.txt')])
    elif opt.pp_steps == 'all':
        steps = all_steps
    else:
        steps = range(*[int(s) for s in opt.pp_steps.split()])
        steps = list(set(steps) & set(all_steps))
    for step in steps:
        print '\n' + 40*'-'
        print 'Postprocessing MC step:', step
        # The current_stage parameter is not used by postprocessing,
        # just put 'melt' here so the C++ code can run in the postprocessing case.
        execute('{} post {} {} melt'.format(molbuild_exe, step, opt.input_file),
                quit_on_error=True)

def run_mc_slip(opt):
    ''' Postprocess LAMMPS data file and deformation trajectory file to compute slip. 
        Connectivity information is obtained from LAMMPS data file,
        Coordinate information is obtained from LAMMPS trajectory file.
        The trajectory file is initially put in the 'slip' folder.
    '''
    steps = glob.glob(os.path.join('slip', '*.lammps'))
    steps = sorted([int(re.search(r'_(\d+)',s).group(1)) for s in steps])
    for step in steps:
        print '\n' + 40*'-'
        print 'Postprocessing slip for MC step:', step
        # The current_stage parameter is not used by postprocessing,
        # just put 'melt' here so the C++ code can run in the postprocessing case.
        execute('{} slip {} {} melt'.format(molbuild_exe, step, opt.input_file),
                quit_on_error=True)

def execute(cmd, quit_on_error=False):
    """ Executes a system command.  If an error is returned from the command
    print an error message and quit.  Command output is written to stdout.
    """
    # Make sure that output is flushed.
    sys.stdout.flush()
    error = subprocess.call(cmd.split())
    if quit_on_error and error:
        print 'Error calling:\n{}'.format(cmd)
        exit(1)
    return error

def init_run_directory(opt):
    ''' Deletes all files generated by builder '''
    subprocess.call(['rm', '-rf', 'data', 'low_energy_state'])
    subprocess.call('rm -f *.lammps *.log accepted_steps'.split())
    os.mkdir('data')
    os.mkdir(os.path.join('data','log'))
    os.mkdir('low_energy_state')

def save_state(opt, step, counter):
    ''' Copies trial state to accepted state and archives every N states.
    '''
    shutil.copy('trial.lammps', 'accepted.lammps')
    sf = int(opt.save_frequency)
    [melt_steps,anneal_steps,sample_steps] = [int(s) for s in opt.steps.split()]
    if step > melt_steps+int(0.95*anneal_steps) and step < melt_steps+anneal_steps:
        sf = 10
    if step >= melt_steps+anneal_steps:
        sf = 2
    if counter%sf == 0:
        shutil.copy('accepted.lammps',
                    os.path.join('data', 'step_{}.lammps'.format(step)))

def read_begin_step():
    ''' Read the last #.out log file, return the number of
        the last performed cycle.
    '''
    # Wait for 5 mins before reading from #.out file.
    time.sleep(300)
    # Determine the last log file.
    path = os.getcwd()
    logs = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    logs = [f for f in logs if 'out' in f]
    submit_id = max([int(f.split('.')[0]) for f in logs])
    last_log = '{0}.out'.format(submit_id)
    # Read from the last log file to get the start cycle number
    # for the new submit.
    fid = open(last_log)
    while True:
        line = fid.readline()
        if line == '': break
        if 'Running MC step:' in line:
            step_begin = int(line.split()[-1])
        elif 'mpiexec: killing job' in line:
            break
    return step_begin

if __name__ == '__main__':
    main()
