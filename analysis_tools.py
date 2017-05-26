import os, itertools
import numpy as np
import cobra
from read_excel import read_excel


class ReactionVariances(object):
    """Holds results of testing the model by varying the reaction bounds of a number of reactions.

    self.data: numpy.Matrix. The main data structure is an (n+1)-dimensional matrix, where n is the number of modified reactions. Each position in the matrix holds a 1-dimensional array, whose 0-th entry is the objective function flux at that condition, the 1-th entry is the total flux in the entire system at that condition, and the following values are the flux through each reaction in self.measured, sorted by the name of the reaction ID.
    """
    def __init__(self, model, to_modify, to_measure, solver=None):
        self.name = str(model) # Name of the model.
        if solver != None:
            self.solver = solver
        else:
            self.solver = cobra.solvers.cglpk
        self.modified = [] # List of reaction ids that were each modified and tested.
        self.measured = [] # List of reaction ids for which data are collected.
        self.measured_descrs = [] # List of human-readable descriptions.
        self._rev_suffix = '_reverse'
        self._epsilon = 1E-5
        self._parse_inputs(to_modify, to_measure)
        self._setup_models(model) # Sets up self.irr_lp, self.irr_model, and self.rev_model
        self.data = np.zeros([m[3] for m in to_modify] + [len(to_measure)+2])
        self._run() # Add print indicating how many conditions will be run.

    # # # # #  Public methods. # # # # #
    def get_flux(self, rxn_id):
        ind = self._measured_ind[rxn_id]
        return self._get_ith_measurement(ind)
    def get_objective_flux(self):
        return self._get_ith_measurement(0)
    def get_total_flux(self):
        return self._get_ith_measurement(1)

    def save(self, file_path):
        pass

    # # # # #  Private methods: input and setup. # # # # #
    def _parse_inputs(self, to_modify, to_measure):
        # Sets self.modified, self.measured, self.measured_descrs, self._measured_ind, self._steps
        if not to_modify or not to_measure:
            print('At least one reaction must be given to modify, and at least one to measure')
            exit()
        self.modified = [m[0] for m in to_modify] # Retains the same order.
        self.measured = list(sorted(to_measure.keys())) # Sorted by reaction id.
        self.measured_descrs = [to_measure[m] for m in self.measured] # Same order as self.measured.
        self._measured_ind = {m:i+2 for i,m in enumerate(self.measured)} # Gives index in self.data.
        self._steps = [self._expand_steps(*m[1:]) for m in to_modify] # A list of lists, where each sublist contains the flux steps for one modify reaction series.
    def _setup_models(self, model):
        # Sets self.irr_lp, self.rev_lp, self.irr_model, self.rev_model, self._obj_ind, self._obj_ub, self._irr_m_rxn_ind, self._rev_m_rxn_ind
        if len(model.objective) != 1:
            print('The model can only have 1 objective function.')
            exit()
        self.rev_model = model.copy()
        self.rev_lp = self.solver.create_problem(self.rev_model, objective_sense='maximize')
        self.rev_lp.solve_problem(objective_sense='maximize')
        if self.solver.get_status(self.rev_lp) != 'optimal':
            print('The model could not be initially optimized')
            exit()
        self.irr_model = model.copy()
        cobra.manipulation.modify.convert_to_irreversible(self.irr_model)
        self.irr_lp = self.solver.create_problem(self.irr_model, objective_sense='minimize')
        self.irr_lp.solve_problem(objective_sense='maximize')
        if self.solver.get_status(self.irr_lp) != 'optimal':
            print('The irreversible model could not be initially optimized')
            exit()
        irr_m_rxn_ind, rev_m_rxn_ind, obj_ind = {}, {}, None
        for r_id in self.modified + self.measured:
            rev_id = r_id + self._rev_suffix
            if r_id not in self.irr_model.reactions:
                print('Could not find reaction "%s" in the model.' % r_id)
                exit()
            elif r_id in irr_m_rxn_ind:
                print('Reaction "%s" was not unique in the model.' % r_id)
                exit()
            irr_m_rxn_ind[r_id] = None
            rev_m_rxn_ind[r_id] = None
            if rev_id in self.irr_model.reactions:
                irr_m_rxn_ind[rev_id] = None
        for i, rxn in enumerate(self.rev_model.reactions):
            if rxn.id in rev_m_rxn_ind:
                rev_m_rxn_ind[rxn.id] = i
        for i, rxn in enumerate(self.irr_model.reactions):
            if rxn.objective_coefficient != 0:
                if rxn.upper_bound <= 0:
                    print('The objective function must have an upper bound > 0.')
                    exit()
                self._obj_ub = rxn.upper_bound
                obj_ind = i
                self.irr_lp.change_variable_objective(i, 0)
            else:
                self.irr_lp.change_variable_objective(i, 1)
            if rxn.id in irr_m_rxn_ind:
                irr_m_rxn_ind[rxn.id] = i
        for r_id, ind in irr_m_rxn_ind.items():
            if ind == None:
                print('"%s" in modify or measured did not have an index in the model' % r_id)
                exit()
        if obj_ind == None:
            print('Could not find an index for the objective function')
            exit()
        self._obj_ind = obj_ind
        self._irr_m_rxn_ind = irr_m_rxn_ind
        self._rev_m_rxn_ind = rev_m_rxn_ind

    # # # # #  Private methods: running the analysis. # # # # #
    def _run(self):
        print('\nMeasuring parsimonious flux in %i reactions from %s, while modifying %i reactions:' % (len(self.measured), self.rev_model, len(self.modified)))
        total_conditions = 1
        for mod, steps in zip(self.modified, self._steps):
            print('  %s: %i conditions from %.1f to %.1f.' % (mod, len(steps), steps[0], steps[-1]))
            total_conditions *= len(steps)
        print('Analysis requires %i optimizations...\n' % (total_conditions))
        for rxn_steps, ind in self._run_steps_iter():
            self._set_reaction_bounds(rxn_steps)
            self._optimize_and_measure(ind, rxn_steps)
    def _run_steps_iter(self):
        # Generator that yields (list, index) at each iteration. The list: [(rxn1_id, rxn1_flux), (rxn2_id, rxn2_flux), ...] for each rxn in to_modify. Every reaction should be constrained to these flux values before one optimization and measurement. The index indicates where the measurements should be stored in self.data.
        for r_fs, ind in itertools.izip(itertools.product(*self._steps), itertools.product( *[range(len(s)) for s in self._steps] )):
            yield zip(self.modified, r_fs), ind
    def _optimize_and_measure(self, data_ind, rxn_steps):
        infeasible = 0.1
        obj_f, total_f = 0.0, 0.0

        self.rev_lp.solve_problem(objective_sense='maximize')

        self.rev_model.solution = self.solver.format_solution(self.rev_lp, self.rev_model)
        if self.rev_model.solution.status == "optimal" and self.rev_model.solution.f > infeasible:
            obj_f = self.rev_model.solution.f
            self.irr_lp.change_variable_bounds(self._obj_ind, obj_f-self._epsilon, self._obj_ub)
            self.irr_lp.solve_problem(objective_sense='minimize') # Minimized for parsimony
            self.irr_model.solution = self.solver.format_solution(self.irr_lp, self.irr_model)
            if self.irr_model.solution.status == "optimal":
                total_f = self.irr_model.solution.f
                self.irr_model.solution.f = obj_f
        self._measure_reactions(data_ind, obj_f, total_f)

    # # # # #  Private methods: getting and setting values. # # # # #
    def _set_reaction_bounds(self, rxn_steps):
        for rxn_id, rxn_f in rxn_steps:
            fwd_f, rev_f = 0.0, 0.0
            if rxn_f > 0:
                fwd_f = rxn_f
            elif rxn_f < 0:
                rev_f = abs(rxn_f)
                rxn_f -= self._epsilon
            self.rev_lp.change_variable_bounds(self._rev_m_rxn_ind[rxn_id], rxn_f, rxn_f+self._epsilon)
            self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rxn_id], fwd_f, fwd_f+self._epsilon)
            rev_id = rxn_id + self._rev_suffix
            if rev_id in self._irr_m_rxn_ind:
                self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rev_id], rev_f, rev_f+self._epsilon)
    def _measure_reactions(self, data_ind, obj_f, total_f):
        self.data[data_ind][0] = obj_f
        self.data[data_ind][1] = total_f
        for i, rxn_id in enumerate(self.measured):
            if obj_f == 0:
                val = 0
            else:
                val = self.irr_model.reactions.get_by_id(rxn_id).x
                rev_id = rxn_id+self._rev_suffix
                if val == 0 and rev_id in self.irr_model.reactions:
                    val = -self.irr_model.reactions.get_by_id(rev_id).x
            self.data[data_ind][i+2] = val
    def _get_ith_measurement(self, ind):
        slc = [slice(None)] * len(self.data.shape)
        slc[-1] = ind
        return self.data[slc]

    # # # # #  Private methods: misc. # # # # #
    def _expand_steps(self, _min, _max, _steps):
        if _steps == 1:
            return [(_min+_max)/2.0]
        delta = (_max - _min) / float(_steps - 1)
        return [_min + delta*x for x in range(_steps)]


if __name__ == '__main__':
    files_dir = '/mnt/hgfs/win_projects/brugia_project'
    model_names = ['model_o_vol_3.5.xlsx', 'model_b_mal_3.5.xlsx']

    to_modify = [('CARBON_SOURCE', 25, 200, 3), ('DIFFUSION_2', -80, -500, 3)]
    #to_modify = [('CARBON_SOURCE', 50, 200, 20), ('FA_SOURCE', 0, 50, 20), ('DIFFUSION_2', 0, -500, 20)]
    to_measure = {'M_TRANS_5':'Mito Pyruvate', 'R01082_M':'Fum -> mal'}

    model_files = [os.path.join(files_dir, m_file) for m_file in model_names]
    models = [read_excel(m_file, verbose=False) for m_file in model_files]
    rv = ReactionVariances(models[0], to_modify, to_measure)

    if False:
        print '\nobjective'
        print rv.get_objective_flux()
        print '\nmito pyr'
        print rv.get_flux('M_TRANS_5')
        print '\nfum -> mal'
        print rv.get_flux('R01082_M')
