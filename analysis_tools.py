import os, itertools
import numpy as np
import cobra
from read_excel import read_excel


class ReactionVariances(object):
    """Holds results of testing the model by varying the reaction bounds of a number of reactions.

    self.data: numpy.Matrix. The main data structure is an (n+1)-dimensional matrix, where n is the number of modified reactions. Each position in the matrix holds a 1-dimensional array, whose 0-th entry is the objective function flux at that condition, the 1-th entry is the total flux in the entire system at that condition, and the following values are the flux through each reaction in self.measured, sorted by the name of the reaction ID.
    """
    def __init__(self, model, to_modify, to_measure, solver):
        self.name = str(model) # Name of the model.
        self.solver = solver
        self.modified = [] # List of reaction ids that were each modified and tested.
        self.measured = [] # List of reaction ids for which data are collected.
        self.measured_descrs = [] # List of human-readable descriptions.
        self._rev_suffix = '_reverse'
        self._parse_inputs(to_modify, to_measure)
        self._setup_models(model) # Sets up self.lp, self.model, and self.obj_model
        self.data = np.zeros([m[3] for m in to_modify] + [len(to_measure)+2])
        self._run()

    # # # # #  Public methods. # # # # #
    def get(self, rxn_id):
        ind = self._measured_ind[rxn_id]
        # return slice of self.data
    def get_objective_flux(self):
        pass
    def get_total_flux(self):
        pass
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
        # Sets self.lp, self.obj_lp, self.model, self.obj_model, self._obj_ind, self._obj_ub, self._rxn_ind, self._obj_rxn_ind
        if len(model.objective) != 1:
            print('The model can only have 1 objective function.')
            exit()
        self.obj_model = model.copy()
        self.obj_lp = self.solver.create_problem(self.obj_model, objective_sense='maximize')
        self.solver.solve_problem(self.obj_lp, objective_sense='maximize')
        if self.solver.get_status(self.obj_lp) != 'optimal':
            print('The model could not be initially optimized')
            exit()
        m = model.copy()
        cobra.manipulation.modify.convert_to_irreversible(m)
        self.lp = self.solver.create_problem(m, objective_sense='maximize')
        self.solver.solve_problem(self.lp, objective_sense='maximize')
        if self.solver.get_status(self.lp) != 'optimal':
            print('The irreversible model could not be initially optimized')
            exit()
        rxn_ind, obj_rxn_ind, obj_ind = {}, {}, None
        for r_id in self.modified + self.measured:
            rev_id = r_id + self._rev_suffix
            if r_id not in m.reactions:
                print('Could not find reaction "%s" in the model.' % r_id)
                exit()
            elif r_id in rxn_ind:
                print('Reaction "%s" was not unique in the model.' % r_id)
                exit()
            rxn_ind[r_id] = None
            obj_rxn_ind[r_id] = None
            if rev_id in m.reactions:
                rxn_ind[rev_id] = None
        for i, rxn in enumerate(self.obj_model.reactions):
            if rxn.id in obj_rxn_ind:
                obj_rxn_ind[rxn.id] = i
        for i, rxn in enumerate(m.reactions):
            if rxn.objective_coefficient != 0:
                if rxn.upper_bound <= 0:
                    print('The objective function must have an upper bound > 0.')
                    exit()
                self._obj_ub = rxn.upper_bound
                obj_ind = i
                self.solver.change_variable_objective(self.lp, i, 0)
            else:
                self.solver.change_variable_objective(self.lp, i, 1)
            if rxn.id in rxn_ind:
                rxn_ind[rxn.id] = i
        for r_id, ind in rxn_ind.items():
            if ind == None:
                print('"%s" in modify or measured did not have an index in the model' % r_id)
                exit()
        if obj_ind == None:
            print('Could not find an index for the objective function')
            exit()
        self.model = m
        self._obj_ind = obj_ind
        self._rxn_ind = rxn_ind
        self._obj_rxn_ind = obj_rxn_ind

    # # # # #  Private methods: running the analysis. # # # # #
    def _run(self):
        print self.modified, self._steps
        for rxn_steps, ind in self._run_steps_iter():
            self._set_reaction_bounds(rxn_steps)
            self._optimize_and_measure(ind)
        print self.data
    def _run_steps_iter(self):
        # Generator that yields (list, index) at each iteration. The list: [(rxn1_id, rxn1_flux), (rxn2_id, rxn2_flux), ...] for each rxn in to_modify. Every reaction should be constrained to these flux values before one optimization and measurement. The index indicates where the measurements should be stored in self.data.
        for r_fs, ind in itertools.izip(itertools.product(*self._steps), itertools.product( *[range(len(s)) for s in self._steps] )):
            yield zip(self.modified, r_fs), ind
    def _optimize_and_measure(self, data_ind):
        infeasible = 0.01
        obj_f, total_f = 0.0, 0.0
        self.solver.solve_problem(self.obj_lp, objective_sense='maximize')
        self.obj_model.solution = self.solver.format_solution(self.obj_lp, self.obj_model)
        if self.obj_model.solution.status == "optimal" and abs(self.obj_model.solution.f) > infeasible:
            obj_f = self.obj_model.solution.f
            self.solver.change_variable_bounds(self.lp, self._obj_ind, obj_f, self._obj_ub)
            self.solver.solve_problem(self.lp, objective_sense='minimize') # Minimized for parsimony
            self.model.solution = self.solver.format_solution(self.lp, self.model)
            if self.model.solution.status == "optimal":
                total_f = self.model.solution.f
                self.model.solution.f = obj_f
        self._measure_reactions(data_ind, obj_f, total_f)

    # # # # #  Private methods: getting and setting reaction values. # # # # #
    def _set_reaction_bounds(self, rxn_steps):
        for rxn_id, rxn_f in rxn_steps:
            self.solver.change_variable_bounds(self.obj_lp, self._obj_rxn_ind[rxn_id], rxn_f, rxn_f)
            fwd_f, rev_f = 0.0, 0.0
            if rxn_f > 0:
                fwd_f = rxn_f
            elif rxn_f < 0:
                rev_f = abs(rxn_f)
            self.solver.change_variable_bounds(self.lp, self._rxn_ind[rxn_id], fwd_f, fwd_f)
            rev_id = rxn_id + self._rev_suffix
            if rev_id in self._rxn_ind:
                self.solver.change_variable_bounds(self.lp, self._rxn_ind[rev_id], rev_f, rev_f)
    def _measure_reactions(self, data_ind, obj_f, total_f):
        self.data[data_ind][0] = obj_f
        self.data[data_ind][1] = total_f
        for i, rxn_id in enumerate(self.measured):
            if obj_f == 0:
                val = 0
            else:
                val = self.model.reactions.get_by_id(rxn_id).x
                rev_id = rxn_id+self._rev_suffix
                if val == 0 and rev_id in self.model.reactions:
                    val = self.model.reactions.get_by_id(rev_id).x
            self.data[data_ind][i+2] = val

    # # # # #  Private methods: misc. # # # # #
    def _expand_steps(self, _min, _max, _steps):
        delta = (_max - _min) / float(_steps - 1)
        return [_min + delta*x for x in range(_steps)]


def analyze_reaction_bounds(model, to_modify, to_measure, solver):
    # Adapted from version 0.5.11 of cobra.flux_analysis.parsimonious
    m = model.copy()
    cobra.manipulation.modify.convert_to_irreversible(m) # So all reactions have positive flux(?)
    lp = solver.create_problem(m, objective_sense="maximize")
    solver.solve_problem(lp, objective_sense='maximize')
    if solver.get_status(lp) != 'optimal':
        print('The model could not be initially optimized')
        exit()
    orig_obj_f = solver.get_objective_value(lp)
    #_setup_solver_for_parsimony(m, orig_obj_f, lp, solver)

    for i, rxn in enumerate(m.reactions):
        if rxn.id == to_modify[0][0]:
            r_ind = i
            break
    for r_id, min_f, max_f, steps in to_modify:
        for rxn_f in _expand_steps(min_f, max_f, steps):
            solver.change_variable_bounds(lp, r_ind, rxn_f, rxn_f)
            #solver.solve_problem(lp, objective_sense='minimize')
            solver.solve_problem(lp, objective_sense='maximize')
            solution = solver.format_solution(lp, m)
            m.solution = solution
            print rxn_f, m.reactions.get_by_id(to_measure[0]).x

    return
    solver.solve_problem(lp, objective_sense='minimize') # Must be minimize for the parsimony
    solution = solver.format_solution(lp, m)
    m.solution = solution
    if solution.status == "optimal":
        m.solution.f = sum([coeff * rxn.x for rxn, coeff in m.objective.items()])
    print m.solution.f, m.reactions.R00431.x, m.reactions.R00431_reverse.x

def _setup_solver_for_parsimony(m, orig_obj_f, lp, solver):
    if len(m.objective) != 1:
        print('The model can only have 1 objective function.')
        exit()
    for i, rxn in enumerate(m.reactions):
        if rxn.objective_coefficient != 0:
            # Sets objective function lower bound to original objective flux.
            target_flux = orig_obj_f / rxn.objective_coefficient
            solver.change_variable_bounds(lp, i, target_flux, rxn.upper_bound)
            solver.change_variable_objective(lp, i, 0)
        else:
            solver.change_variable_objective(lp, i, 1) # Adds all reactions as part of the objective function.
def _expand_steps(_min, _max, _steps):
    delta = (_max - _min) / float(_steps - 1)
    return [_min + delta*x for x in range(_steps)]

if __name__ == '__main__':
    files_dir = '/mnt/hgfs/win_projects/brugia_project'
    model_names = ['model_o_vol_3.5.xlsx', 'model_b_mal_3.5.xlsx']
    solver = cobra.solvers.cglpk

    model_files = [os.path.join(files_dir, m_file) for m_file in model_names]
    models = [read_excel(m_file, verbose=False) for m_file in model_files]

    to_modify = [('CARBON_SOURCE', 50, 200, 3)]
    to_modify = [('CARBON_SOURCE', 0, 200, 20), ('FA_SOURCE', 0, 50, 20)]
    to_measure = {'M_TRANS_5':'Mito Pyruvate'}
    #analyze_reaction_bounds(models[0], to_modify, to_measure, solver)

    rv = ReactionVariances(models[0], to_modify, to_measure, solver)
