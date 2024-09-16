import math
import random
import bisect
import lammps

# Boltzmann constant in kcal/mol/K
kB   = 1.9872156e-3    

class MetropolisAlgorithm:
    ''' 
    '''
    def __init__(self, opt):
        ''' Initializes parameters of the Metropolis algorithm from input '''
        self.melt_temperature   = float(opt.melt_T)
        self.target_temperature = float(opt.target_T)

        steps = (int(s) for s in opt.steps.split())
        self.melt_steps,self.anneal_steps,self.sample_steps = steps

        self.last_accepted_step = None

    def temperature(self, step):
        ''' Computes the temperature based on the current step '''
        if step < self.melt_steps:
            return self.melt_temperature
        elif step - self.melt_steps < self.anneal_steps:
            x = float(step-self.melt_steps) / self.anneal_steps
            dT = self.melt_temperature - self.target_temperature
            return self.melt_temperature - x*dT
        return self.target_temperature

    def is_move_accepted(self, step, E_start, E_trial, max_energy_drop):
        ''' Returns true or false, whether or not to accept this step. '''
        # Always accept the first step.
        if self.last_accepted_step is None:
            self.last_accepted_step = step
            return True

        # Handle case where LAMMPS failed to complete step.
        if E_trial is None:
            print 'LAMMPS did not complete successfully'
            return False

        print 'Current H (End_cycle H):  {:.3f} kcal/mol'.format(E_trial)
        print 'Previous H (Start_cycle H): {:.3f} kcal/mol'.format(E_start)
        delta_E = E_trial - E_start
      
        nA,nS = (float(s) for s in open('num_candidates').read().split())
        T = self.temperature(step)
        print 'Temperature is', T, 'K'
        q = min(-delta_E/(kB*T), 10.0)
        p = min(math.exp(q)*(nA/nS), 1.0)
        
        print 'MC move acceptance probability: {:.6f}'.format(p)
        # Always accept the very first few steps.
        if step <= 10:
            self.last_accepted_step = step
            return True
        # Perform metropolis criterion.
        if p > random.random():
            # Unfavor lower energy in melt stage only.
            if step < self.melt_steps:
                # Skip the very first few cycles.
                if step > 20:
                    if delta_E < -max_energy_drop:
                        print 'max energy drop test failed.'
                        return False
            self.last_accepted_step = step
            return True
        else:
            return False

    def lookup_last_accepted_step(self, step):
        ''' Reads accepted step log and last prior accepted step. '''
        A = [int(s) for s in open('accepted_steps').read().split()] 
        self.last_accepted_step = max(0, A[bisect.bisect_left(A, step)-1])
        print 'Last accepted step was {0}'.format(self.last_accepted_step)


if __name__ == '__main__':
    ''' Debugging code '''
    class Options: pass
    opt = Options()
    opt.melt_T   = '3000'
    opt.target_T = '300'
    opt.steps    = '5 40 5'
    mc = MetropolisAlgorithm(opt)

    for i in range(50):
        print i, mc.temperature(i)
