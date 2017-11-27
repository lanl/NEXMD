#!/usr/bin/python

'''

The header is defined as a class and the attributes
of header are all input parameters.

'''

import numpy as np

class header(object):
    
    def __init__(self, path):
        self.path = path
        file = open(path,'r')
        #self.file = file = file.readlines()
        self.file = file.readlines()
        
        line_num = 0
        for line in self.file:
            ## Determine lines defining the qmmm block ##
            if '&qmmm' in line:
                qmmm_block_start = line_num
            if '&endqmmm' in line:
                qmmm_block_end  = line_num
            ## Determine lines defining the moldyn block ##
            if '&moldyn' in line:
                moldyn_block_start = line_num
            if '&endmoldyn' in line
                moldyn_block_end = line_num
            line_num += 1
            
        line_num = 0
        for line in self.file:
            '''
            ## Start geometry optimization ##
            if 'maxcyc=' in line and 'dav_maxcyc=' not in line:
                self.maxcyc = np.int(line.split()[0][len('maxcyc='):-1])
            if 'ntpr=' in line:
                self.ntpr = np.int(line.split()[0][len('ntpr='):-1])
            if 'grms_tol=' in line:
                self.grms_tol = np.float(line.split()[0][len('grms_tol='):-1])
            ## End geometry optimization ##

            ## Start ground-state and output parameters ##
            if 'qm_theory=' in line:
                self.qm_theory = np.str(line.split()[0][len('qm_theory='):-1])
            if 'scfconv=' in line:
                self.scfconv = np.float(line.split()[0][len('scfconv='):-1])
            if 'verbosity=' in line:
                self.verbosity = np.int(line.split()[0][len('verbosity='):-1])
            if 'printdipole=' in line:
                self.printdipole = np.int(line.split()[0][len('printdipole='):-1])
            if 'printbondorders=' in line:
                self.printbondorders = np.int(line.split()[0][len('printbondorders='):-1])
            if 'density_predict=' in line:
                self.density_predict = np.int(line.split()[0][len('density_predict='):-1])
            if 'itrmax=' in line:
                self.itrmax = np.int(line.split()[0][len('itrmax='):-1])
            ## End ground-state and output parameters ##

            ## Start excited-state parameters ##
            if 'exst_method=' in line:
                self.exst_method = np.int(line.split()[0][len('exst_method='):-1])
            if 'dav_guess=' in line:
                self.dav_guess = np.int(line.split()[0][len('dav_guess='):-1])
            if 'ftol0=' in  line:
                self.ftol0 = np.float(line.split()[0][len('ftol0='):-1])
            if 'ftol1=' in line:
                self.ftol1 = np.float(line.split()[0][len('ftol1='):-1])
            if 'dav_maxcyc=' in line:
                self.dav_maxcyc = np.int(line.split()[0][len('dav_maxcyc='):-1])
            if 'printcharges=' in line:
                self.printcharges = np.int(line.split()[0][len('printcharges='):-1])
            if 'calcxdens=' in line:
                self.calcxdens = np.str(line.split()[0][len('calcxdens='):-1])
            ## End excited-state parameters ##

            ## Start solvent models and external electric fields ##
            if 'solvent_model=' in line:
                self.solvent_model = np.int(line.split()[0][len('solvent_model='):-1])
            if 'potential_type=' in line:
                self.potential_type = np.int(line.split()[0][len('potential_type='):-1])
            if 'onsager_radius=' in line:
                self.onsager_radius = np.float(line.split()[0][len('onsager_radius='):-1])
            if 'ceps=' in line:
                self.ceps = np.float(line.split()[0][len('ceps='):-1])
            if 'linmixparam=' in line:
                self.linmixparam = np.float(line.split()[0][len('linmixparam='):-1])
            if 'cosmo_scf_ftol=' in line:
                self.cosmo_scf_ftol = np.float(line.split()[0][len('cosmo_scf_ftol='):-1])
            if 'doZ=' in line:
                self.doZ = np.str(line.split()[0][len('doZ='):-1])
            if 'index_of_refraction=' in line:
                self.index_of_refraction = np.float(line.split()[0][len('index_of_refraction='):-1])
            if 'EF=' in line:
                self.EF = np.int(line.split()[0][len('EF='):-1])
            if 'Ex=' in line:
                self.Ex = np.float(line.split()[0][len('Ex='):-1])
            if 'Ey=' in line:
                self.Ey = np.float(line.split()[0][len('Ey='):-1])
            if 'Ez=' in line:
                self.Ez = np.float(line.split()[0][len('Ez='):-1])
            ## End solvent models and external electric fields ##
            '''
            ## Start general parameters ##
            if 'natoms=' in line:
                self.natoms = np.int(line.split()[0][len('natoms='):-1])
            '''
            if 'rnd_seed=' in line:
                self.rnd_seed = np.int(line.split()[0][len('rnd_seed='):-1])
            if 'bo_dynamics_flag=' in line:
                self.bo_dynamics = np.int(line.split()[0][len('bo_dynamics_flag='):-1])
            if 'exc_state_init=' in line:
                self.exc_state_init = np.int(line.split()[0][len('exc_state_init='):-1])
            '''
            if 'n_exc_states_propagate=' in line:
                self.n_exc_states_propagate = np.int(line.split()[0][len('n_exc_states_propagate='):-1])
            ## End general parameters ###
            
            ## Start dynamics parameters ##
            if 'time_init=' in line:
                self.time_init = np.float(line.split()[0][len('time_init='):-1])
            if 'time_step=' in line:
                self.time_step = np.float(line.split()[0][len('time_step='):-1])
            if 'n_class_steps=' in line:
                self.n_class_steps = np.int(line.split()[0][len('n_class_steps='):-1])
            if 'n_quant_steps=' in line:
                self.n_quant_steps = np.int(line.split()[0][len('n_quant_steps='):-1])
            '''
            if 'moldyn_deriv_flag=' in line:
                self.moldyn_deriv_flag = np.int(line.split()[0][len('moldyn_deriv_flag='):-1])
            if 'num_deriv_step=' in line:
                self.num_deriv_step = np.float(line.split()[0][len('num_deriv_step='):-1])
            if 'rk_tolerance=' in line:
                self.rk_tolerance = np.float(line.split()[0][len('rk_tolerance='):-1])
            ## End dynamics parameters ##

            ## Start Nonadiabatic parameters ##
            if 'decoher_type=' in line:
                self.decoher_type = np.int(line.split()[0][len('decoher_type='):-1])
            if 'decoher_e0=' in line:
                self.decoher_e0 = np.float(line.split()[0][len('decoher_e0='):-1])
            if 'decoher_c=' in line:
                self.decoher_c = np.float(line.split()[0][len('decoher_c='):-1])
            if 'dotrivial=' in line:
                self.dotrivial = np.int(line.split()[0][len('dotrivial='):-1])
            if 'quant_step_reduction_factor=' in line:
                self.quant_step_reduction_factor = np.float(line.split()[0][len('quant_step_reduction_factor='):-1])
            ## End Nonadiabatic parameters ##

            ## Start thermostat parameters ##
            if 'therm_type=' in line:
                self.therm_type = np.int(line.split()[0][len('therm_type='):-1])
            if 'therm_temperature=' in line:
                self.therm_temperature = np.float(line.split()[0][len('therm_temperature='):-1])
            if 'therm_friction=' in line:
                self.therm_friction = np.float(line.split()[0][len('therm_friction='):-1])
            if 'berendsen_relax_const=' in line:
                self.berendsen_relax_const = np.float(line.split()[0][len('berendsen_relax_const='):-1])
            if 'heating=' in line:
                self.heating = np.int(line.split()[0][len('heating='):-1])
            if 'heating_steps_per_degree=' in line:
                self.heating_steps_per_degree = np.int(line.split()[0][len('heating_steps_per_degree='):-1])
            ## End thermostat parameters ##
            '''
            ## Start output and log parameters ##
            if 'verbosity=' in line and moldyn_block_start < line_num < moldyn_block_end:
                self.moldyn_verbosity = np.int(line.split()[0][len('verbosity='):-1])
            if 'out_data_steps=' in line:
                self.out_data_steps = np.int(line.split()[0][len('out_data_steps='):-1])
            if 'out_coords_steps=' in line:
                self.out_coords_steps = np.int(line.split()[0][len('out_coords_steps='):-1])
            '''
            if 'out_data_cube=' in line:
                self.out_data_cube = np.int(line.split()[0][len('out_data_cube='):-1])
            if 'out_count_init=' in line:
                self.out_count_init = np.int(line.split()[0][len('out_count_init='):-1])
            ## End output and log parameters ##
            '''
            ## Check is coefficients are set ##
            if 'quant_amp_phase' in line:
                header.quant_amp_phase = 'The quant_amp_phase flag is in the header.'

            line_num += 1
