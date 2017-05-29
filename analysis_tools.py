import os, itertools
import numpy as np
import cobra
from read_excel import read_excel
import matplotlib.pyplot as plt


class ReactionVariances(object):
    """Holds results of testing the model by varying the reaction bounds of a number of reactions.

    self.data: numpy.Matrix. The main data structure is an (n+1)-dimensional matrix, where n is the number of modified reactions. Each position in the matrix holds a 1-dimensional array, whose 0-th entry is the objective function flux at that condition, the 1-th entry is the total flux in the entire system at that condition, and the following values are the flux through each reaction in self.measured, sorted by the name of the reaction ID.
    """
    def __init__(self, model, to_modify, to_measure, solver=None):
        self.name = str(model) # Name of the model.
        self.modified = [] # List of reaction ids that were each modified and tested.
        self.modified_attrs = {} # Dict of properties used for graphing.
        self.measured = [] # List of reaction ids for which data are collected.
        self.measured_attrs = {} # List of human-readable descriptions.
        self.infeasible = 0.1 # An objective value below this is considered infeasible.
        self.epsilon = 1E-6 # Required to deal with floating point errors.
        self._rev_suffix = '_reverse'
        if solver != None:
            self.solver = solver
        else:
            self.solver = cobra.solvers.cglpk
        self._parse_inputs(to_modify, to_measure)
        self._setup_models(model) # Sets up self.irr_lp, self.irr_model, and self.rev_model
        self.data = np.zeros([m[4] for m in to_modify] + [len(to_measure)+2])
        self._run()

    # # # # #  Public methods. # # # # #
    def get_flux(self, rxn_id):
        ind = self._measured_ind[rxn_id]
        return self._get_ith_measurement(ind)
    def get_objective_flux(self):
        return self._get_ith_measurement(0)
    def get_total_flux(self):
        return self._get_ith_measurement(1)

    def negative_modified(self, modified_rxn):
        """Report negative values for this reaction."""
        self.modified_attrs[modified_rxn]['coefficient'] = -1.0
    def negative_measured(self, measured_rxn):
        """Report negative values for this reaction."""
        self.measured_attrs[measured_rxn]['coefficient'] = -1.0
    def update_modified_attrs(self, modified_rxn, attr_dict):
        self.modified_attrs[modified_rxn].update(attr_dict)

    def heatmaps(self, measured_rxns, colormap='RdBu', interpolation='quadric'):
        # interpolation = 'spline36', 'hanning'
        figsize = (len(measured_rxns)*3.5, 3)
        if len(self.modified) != 2:
            print('Error: heatmaps() can only be called with 2 dimensional data.')
            exit()
        for m_rxn in measured_rxns:
            if m_rxn not in self.measured:
                print('Error: the given reaction "%s" was not measured in the current data.' % m_rxn)
                exit()
        fig, axs = plt.subplots(1, len(measured_rxns), figsize=figsize, sharey=True)
        if len(measured_rxns) == 1:
            axs = [axs]
        for ind, (m_rxn, ax) in enumerate(zip(measured_rxns, axs)):
            self.draw_heatmap(m_rxn, fig, ax, colormap, interpolation)
            if ind == 0: # Specific to left-most plot.
                pass
            else: # All others.
                ax.tick_params('y', which='both', left='off', labelleft='off')
        plt.tight_layout()
        plt.show()
    def draw_heatmap(self, measured_rxn, fig, axis, colormap, interpolation):
        min_flux_range = 5.0 # Affects the colouring range.
        fluxes = self.get_flux(measured_rxn) * self.measured_attrs[measured_rxn]['coefficient']
        max_flux = max(abs(fluxes.max()), abs(fluxes.min()), min_flux_range)
        x_attrs = self.modified_attrs[self.modified[1]]
        y_attrs = self.modified_attrs[self.modified[0]]
        img = axis.imshow(fluxes, origin='lower', aspect='auto', vmin=-max_flux, vmax=max_flux, cmap=colormap, interpolation=interpolation)
        fig.colorbar(img, ax=axis, shrink=0.75, aspect=10, ticks=[-int(max_flux), 0, int(max_flux)])
        axis.set_title(self.measured_attrs[measured_rxn]['label'])
        axis.set_xlabel(x_attrs['label'])
        axis.set_ylabel(y_attrs['label'])
        x_vals = [round(s*x_attrs['coefficient'], 1) for s in self._steps[1]]
        x_mid, x_mid_ind = ((len(x_vals)+1) // 2) - 1, (len(x_vals) - 1) / 2.0
        x_mid_val = round((x_vals[x_mid] + x_vals[-x_mid-1]) / 2.0, 1)
        y_vals = [round(s*y_attrs['coefficient'], 1) for s in self._steps[0]]
        y_mid, y_mid_ind = ((len(y_vals)+1) // 2) - 1, (len(y_vals) - 1) / 2.0
        y_mid_val = round((y_vals[y_mid] + y_vals[-y_mid-1]) / 2.0, 1)
        plt.sca(axis)
        plt.xticks([0, x_mid_ind, len(x_vals)-1], [x_vals[0], x_mid_val, x_vals[-1]])
        plt.yticks([0, y_mid_ind, len(y_vals)-1], [y_vals[0], y_mid_val, y_vals[-1]])

    def save(self, file_path):
        pass

    # # # # #  Private methods: running the analysis. # # # # #
    def _run(self):
        print('\nMeasuring parsimonious flux in %i reactions from %s, while modifying %i reactions:' % (len(self.measured), self.name, len(self.modified)))
        total_conditions = 1
        for mod, steps in zip(self.modified, self._steps):
            print('- %s: %i conditions from %.1f to %.1f.' % (mod, len(steps), steps[0], steps[-1]))
            total_conditions *= len(steps)
        print('Analysis requires %i optimizations...' % (total_conditions))
        for rxn_steps, ind in self._run_steps_iter():
            self._set_reaction_bounds(rxn_steps)
            self._optimize_and_measure(ind, rxn_steps)
        print('Analysis complete.\n')
    def _optimize_and_measure(self, data_ind, rxn_steps):
        obj_f, total_f = 0.0, 0.0
        self.rev_lp.solve_problem(objective_sense='maximize')
        self.rev_model.solution = self.solver.format_solution(self.rev_lp, self.rev_model)
        if self.rev_model.solution.status == "optimal" and self.rev_model.solution.f > self.infeasible:
            obj_f = self.rev_model.solution.f
            self.irr_lp.change_variable_bounds(self._obj_ind, obj_f-self.epsilon, self._obj_ub)
            self.irr_lp.solve_problem(objective_sense='minimize') # Minimized for parsimony
            self.irr_model.solution = self.solver.format_solution(self.irr_lp, self.irr_model)
            if self.irr_model.solution.status == "optimal":
                total_f = self.irr_model.solution.f
                self.irr_model.solution.f = obj_f
            else:
                print('Error: The internal linear program solvers disagree about feasibility, which is normally due to floating point errors. This may be solved by increasing the self.epsilon value.')
                exit()
        self._measure_reactions(data_ind, obj_f, total_f)

    # # # # #  Private methods: getting and setting values. # # # # #
    def _set_reaction_bounds(self, rxn_steps):
        for rxn_id, rxn_f in rxn_steps:
            fwd_f, rev_f = 0.0, 0.0
            if rxn_f > 0:
                fwd_f = rxn_f
            elif rxn_f < 0:
                rev_f = abs(rxn_f)
                rxn_f -= self.epsilon
            self.rev_lp.change_variable_bounds(self._rev_m_rxn_ind[rxn_id], rxn_f, rxn_f+self.epsilon)
            self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rxn_id], fwd_f, fwd_f+self.epsilon)
            rev_id = rxn_id + self._rev_suffix
            if rev_id in self._irr_m_rxn_ind:
                self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rev_id], rev_f, rev_f+self.epsilon)
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

    # # # # #  Private methods: input and setup. # # # # #
    def _parse_inputs(self, to_modify, to_measure):
        # Sets self.modified, self.measured, self.measured_attrs, self._measured_ind, self._steps
        if not to_modify or not to_measure:
            print('At least one reaction must be given to modify, and at least one to measure')
            exit()
        self.modified = [m[0] for m in to_modify] # Retains the same order.
        self.modified_attrs = {m[0]:{'label':m[1], 'coefficient':1.0} for m in to_modify}
        self.measured = list(sorted(to_measure.keys())) # Sorted by reaction id.
        self.measured_attrs = {m:{'label':n, 'coefficient':1.0} for m,n in to_measure.items()}
        self._measured_ind = {m:i+2 for i,m in enumerate(self.measured)} # Gives index in self.data.
        self._steps = [self._expand_steps(*m[2:]) for m in to_modify] # A list of lists, where each sublist contains the flux steps for one modify reaction series.
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

    # # # # #  Private methods: misc. # # # # #
    def _expand_steps(self, _min, _max, _steps):
        if _steps == 1:
            return [(_min+_max)/2.0]
        delta = (_max - _min) / float(_steps - 1)
        return [_min + delta*x for x in range(_steps)]
    def _run_steps_iter(self):
        # Generator that yields (list, index) at each iteration. The list: [(rxn1_id, rxn1_flux), (rxn2_id, rxn2_flux), ...] for each rxn in to_modify. Every reaction should be constrained to these flux values before one optimization and measurement. The index indicates where the measurements should be stored in self.data.
        modify_conditions = itertools.product(*self._steps)
        data_inds = itertools.product( *[range(len(s)) for s in self._steps] )
        for r_fs, ind in itertools.izip(modify_conditions, data_inds):
            labeled_conditions = zip(self.modified, r_fs)
            yield labeled_conditions, ind


if __name__ == '__main__':
    files_dir = '/mnt/hgfs/win_projects/brugia_project'
    model_names = ['model_o_vol_3.5.xlsx', 'model_b_mal_3.5.xlsx']

    to_modify = [('CARBON_SOURCE', 'Glucose', 0, 200, 20), ('DIFFUSION_2', 'Oxygen', 0, -500, 20)]
    #to_modify = [('CARBON_SOURCE', 'Glucose', 50, 200, 10), ('FA_SOURCE', 'Fatty acids', 0, 50, 10), ('DIFFUSION_2', 'Oxygen', 0, -500, 10)]
    to_measure = {'M_TRANS_5':'Mitochondrial pyruvate', 'R01082_M':'Fum -> mal', 'SINK_3':'Acetate waste', 'SINK_4':'Propanoate waste'}
    to_display = ['M_TRANS_5', 'R01082_M', 'SINK_3', 'SINK_4']
    #to_display = ['R01082_M']

    model_files = [os.path.join(files_dir, m_file) for m_file in model_names]
    models = [read_excel(m_file, verbose=False) for m_file in model_files]
    rv = ReactionVariances(models[0], to_modify, to_measure)
    rv.negative_modified('DIFFUSION_2')
    rv.negative_measured('R01082_M')

    rv.heatmaps(to_display)
