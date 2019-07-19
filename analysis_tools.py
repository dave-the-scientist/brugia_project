import os, itertools
from math import ceil, floor, sqrt
import numpy as np
import cobra
from model_tools import load_model, set_model_wolbachia
from read_excel import read_excel
import matplotlib
import matplotlib.pyplot as plt


class DiscreteModelResults(object):
    """Holds results of testing the model by varying the reaction bounds of a number of reactions.

    to_modify: A list of tuples.
    self.data: numpy.Matrix. The main data structure is an (n+1)-dimensional matrix, where n is the number of modified reactions. Each position in the matrix holds a 1-dimensional array, whose 0-th entry is the objective function flux at that condition, the 1-th entry is the total flux in the entire system at that condition, and the following values are the flux through each reaction in self.measured, sorted by the name of the reaction ID.
    """
    def __init__(self, model, to_modify, to_measure, mtb_update=None, tight_bounds=False, solver=None):
        self.name = str(model) # Name of the model.
        self.modified = [] # List of reaction groups, each a list of ids that will be modified and tested.
        self.modified_attrs = [] # List of dicts, with properties used for graphing.
        self.measured = [] # List of reaction ids for which data are collected.
        self.measured_attrs = {} # Dict of properties used for graphing.
        self.tight_bounds = tight_bounds # If True lb=ub on modified reactions; else lb=0 (ub=0 for reverse reactions).
        self.infeasible = 0.1 # An objective value below this is considered infeasible.
        self.epsilon = 1E-6 # Required to deal with floating point errors.
        self._rxn_rev_suffix = '_reverse'
        self._obj_fxn_id = '_objective_function'
        self._obj_fxn_label = 'Objective function'
        self._total_flux_rxn_id = '_total_system_flux'
        self._total_flux_rxn_label = 'Total system flux'
        if solver != None:
            self.solver = solver
        else:
            self.solver = cobra.solvers.cglpk
        self._parse_inputs(to_modify, to_measure)
        self._setup_models(model, mtb_update) # Sets up self.irr_lp, self.irr_model, and self.rev_model
        self.data = np.zeros([m[2][2] for m in to_modify] + [len(to_measure)+2])
        self._run()

    # # # # #  Public methods. # # # # #
    def get_fluxes(self, measured_rxn, subslice=None):
        if measured_rxn == self._obj_fxn_id:
            fluxes = self.get_objective_fluxes() * self.measured_attrs[measured_rxn]['coefficient']
        elif measured_rxn == self._total_flux_rxn_id:
            fluxes = self.get_total_fluxes() * self.measured_attrs[measured_rxn]['coefficient']
        else:
            rxn_ind = self._measured_ind[measured_rxn]
            fluxes = self._get_ith_measurement(rxn_ind) * self.measured_attrs[measured_rxn]['coefficient']
        if subslice:
            fluxes = fluxes[subslice]
        return fluxes
    def get_objective_fluxes(self):
        return self._get_ith_measurement(0)
    def get_total_fluxes(self):
        return self._get_ith_measurement(1)

    def negative_modified(self, modified_label):
        """Takes either a reaction ID as a string, or a sequence of reaction ID strings. Causes all reaction IDs to be reported as negatives."""
        if isinstance(modified_label, basestring):
            self._set_modified_attr(modified_label, 'coefficient', -1.0)
        else:
            for m_lab in modified_label:
                self._set_modified_attr(m_lab, 'coefficient', -1.0)
    def negative_measured(self, measured_rxn):
        """Takes either a reaction ID as a string, or a sequence of reaction ID strings. Causes all reaction IDs to be reported as negatives."""
        if isinstance(measured_rxn, basestring):
            if measured_rxn in self.measured_attrs:
                self.measured_attrs[measured_rxn]['coefficient'] = -1.0
        else:
            for m_rxn in measured_rxn:
                if m_rxn in self.measured_attrs:
                    self.measured_attrs[m_rxn]['coefficient'] = -1.0

    def save(self, file_path):
        pass

    # # # # #  Private methods: running the analysis. # # # # #
    def _run(self):
        print('\nMeasuring parsimonious flux in %i reactions from %s, while modifying %i reaction groups:' % (len(self.measured), self.name, len(self.modified)))
        total_conditions = 1
        for mod_attr, steps in zip(self.modified_attrs, self._steps):
            print('- %s: %i conditions from %.1f to %.1f.' % (mod_attr['label'], len(steps), steps[0], steps[-1]))
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
        if self.rev_model.solution.status == 'infeasible':
            # Sometimes the solution will return infeasible for no apparent reason: it would using [('CARBON_SOURCE', 'Glucose', (1, 251, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25)), ('FA_SOURCE', 'Fatty acids', (40, 80, 6))], but not for ('FA_SOURCE', 'Fatty acids', (40, 80, 5)). This is a slower way to solve, but doesn't return infeasible for those same parameters. I no longer think this is needed, and it throws errors sometimes if the configuration is infeasible.
            """for rxn_id, rxn_f in rxn_steps:
                rxn_lb, rxn_ub, _, _, _, _ = self._calc_new_rxn_bounds(rxn_id, rxn_f)
                self.rev_model.reactions.get_by_id(rxn_id).bounds = (rxn_lb, rxn_ub)
            cobra.flux_analysis.optimize_minimal_flux(self.rev_model).status"""
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
    def _calc_new_rxn_bounds(self, rxn_id, rxn_f):
        if rxn_id in self._p_orig_bounds:
            rxn_lb, rxn_ub = self._p_orig_bounds[rxn_id][0]*rxn_f/100.0, self._p_orig_bounds[rxn_id][1]*rxn_f/100.0
            fwd_f = max(rxn_ub, 0.0)
            fwd_f_lb = max(rxn_lb, 0.0)
            rev_f = abs(min(rxn_lb, 0.0))
            rev_f_lb = abs(min(rxn_ub, 0.0))
        else:
            fwd_f, rev_f = 0.0, 0.0
            if rxn_f >= 0:
                fwd_f = rxn_f + self.epsilon
            elif rxn_f < 0:
                rev_f = abs(rxn_f) + self.epsilon
                rxn_f -= self.epsilon
            if self.tight_bounds:
                fwd_f_lb, rev_f_lb = fwd_f-self.epsilon, rev_f-self.epsilon
                rxn_lb, rxn_ub = rxn_f, rxn_f+self.epsilon
            else:
                fwd_f_lb, rev_f_lb = 0.0, 0.0
                rxn_lb, rxn_ub = min(rxn_f, 0.0), max(rxn_f+self.epsilon, 0.0)
        return rxn_lb, rxn_ub, fwd_f_lb, fwd_f, rev_f_lb, rev_f
    def _set_reaction_bounds(self, rxn_steps):
        for rxn_id, rxn_f in rxn_steps:
            rxn_lb, rxn_ub, fwd_f_lb, fwd_f, rev_f_lb, rev_f = self._calc_new_rxn_bounds(rxn_id, rxn_f)
            self.rev_lp.change_variable_bounds(self._rev_m_rxn_ind[rxn_id], rxn_lb, rxn_ub)
            self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rxn_id], fwd_f_lb, fwd_f)
            rev_id = rxn_id + self._rxn_rev_suffix
            if rev_id in self._irr_m_rxn_ind:
                self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rev_id], rev_f_lb, rev_f)
            if rxn_id in self._mtb_update:
                new_coef = self._mtb_update[rxn_id]['mtb_orig_coef']*rxn_f/100.0
                self.rev_lp.change_coefficient(self._mtb_update[rxn_id]['mtb_rev_ind'], self._mtb_update[rxn_id]['rev_ind'], new_coef)
                self.irr_lp.change_coefficient(self._mtb_update[rxn_id]['mtb_irr_ind'], self._mtb_update[rxn_id]['irr_ind'], new_coef)

    def _measure_reactions(self, data_ind, obj_f, total_f):
        self.data[data_ind][0] = obj_f
        self.data[data_ind][1] = total_f
        for i, rxn_id in enumerate(self.measured):
            if obj_f == 0:
                val = 0.0
            else:
                val = self.irr_model.reactions.get_by_id(rxn_id).x
                rev_id = rxn_id + self._rxn_rev_suffix
                if abs(val) < self.epsilon and rev_id in self.irr_model.reactions:
                    val = -self.irr_model.reactions.get_by_id(rev_id).x
            if abs(val) < self.epsilon:
                val = 0.0
            self.data[data_ind][i+2] = val
    def _get_ith_measurement(self, ind):
        slc = [slice(None)] * len(self.data.shape)
        slc[-1] = ind
        return self.data[slc]
    def _set_modified_attr(self, _label, _property, _value):
        for group_prop in self.modified_attrs:
            if group_prop['label'] == _label:
                group_prop[_property] = _value
                break
        else:
            print('Error: reaction group "%s" was not found so could not be updated.' % (_label))
            exit()

    # # # # #  Private methods: input and setup. # # # # #
    def _parse_inputs(self, to_modify, to_measure):
        # Sets self.modified, self.measured, self.measured_attrs, self._p_orig_bounds, self._mtb_update, self._measured_ind, self._steps
        default_modified_attrs = {'coefficient':1.0}
        default_measured_attrs = {'coefficient':1.0}
        if not to_modify or not to_measure:
            print('Error: at least one reaction must be given to modify, and at least one to measure')
            exit()
        self.modified = []
        self.modified_attrs, self.measured_attrs, self._p_orig_bounds = [], {}, {}
        for m in to_modify: # Retains the given order.
            if isinstance(m[0], basestring): # If it is a rxnID string.
                self.modified.append( [m[0]] )
            else: # If it is a list of rxnID strings.
                self.modified.append( m[0] )
            self.modified_attrs.append(dict(default_modified_attrs, label=m[1]))
            if len(m[2]) == 4 and m[2][3] == 'percent':
                for r_id in self.modified[-1]:
                    self._p_orig_bounds[r_id] = None
        self.measured = list(sorted(to_measure.keys())) # Sorted by reaction id.
        while self._obj_fxn_id in self.measured:
            self._obj_fxn_id = '_' + self._obj_fxn_id
        while self._total_flux_rxn_id in self.measured:
            self._total_flux_rxn_id = '_' + self._total_flux_rxn_id
        for m, n in to_measure.items():
            self.measured_attrs.setdefault(m, {'label':n}).update(default_measured_attrs)
        self.measured_attrs.setdefault(self._obj_fxn_id, {'label':self._obj_fxn_label}).update(default_measured_attrs)
        self.measured_attrs.setdefault(self._total_flux_rxn_id, {'label':self._total_flux_rxn_label}).update(default_measured_attrs)
        self._measured_ind = {m:i+2 for i,m in enumerate(self.measured)} # Gives index in self.data.
        self._steps = [self._expand_steps(*m[2][:3]) for m in to_modify] # A list of lists, where each sublist contains the flux steps for one modify reaction series.

    def _setup_models(self, model, mtb_update):
        # Sets self.irr_lp, self.rev_lp, self.irr_model, self.rev_model, self._obj_ind, self._obj_ub, self._irr_m_rxn_ind, self._rev_m_rxn_ind
        if len(model.objective) != 1:
            print('The model can only have 1 objective function.')
            exit()
        self.rev_model = model.copy()
        self.irr_model = model.copy()

        self.rev_lp = self.solver.create_problem(self.rev_model, objective_sense='maximize')
        self.rev_lp.solve_problem(objective_sense='maximize') # Make sure it's feasible.
        if self.solver.get_status(self.rev_lp) != 'optimal':
            print('The model could not be initially optimized')
            exit()
        cobra.manipulation.modify.convert_to_irreversible(self.irr_model)
        self.irr_lp = self.solver.create_problem(self.irr_model, objective_sense='minimize')
        self.irr_lp.solve_problem(objective_sense='maximize') # Make sure it's feasible.
        if self.solver.get_status(self.irr_lp) != 'optimal':
            print('The irreversible model could not be initially optimized')
            exit()

        # While I don't know why, the initial bounds of the to_modify reactions can impact the pars flux values of certain reactions (only those that are pixely, and behaving stochasitcally). The loop below should standardize the behaviour between runs.
        for rxn_group, fluxes in zip(self.modified, self._steps):
            init_flux = fluxes[0]
            if self.tight_bounds:
                bounds = (init_flux-self.epsilon, init_flux+self.epsilon)
            else:
                bounds = tuple(sorted([0.0, init_flux]))
            for r_id in rxn_group:
                if r_id in self._p_orig_bounds: # Means it gets modified by relative amount, not set to an absolute.
                    bnds = model.reactions.get_by_id(r_id).bounds
                    self._p_orig_bounds[r_id] = bnds
                    bounds = (bnds[0]*init_flux/100.0, bnds[1]*init_flux/100.0)
                self.rev_model.reactions.get_by_id(r_id).bounds = bounds
                self.irr_model.reactions.get_by_id(r_id).bounds = bounds


        irr_m_rxn_ind, rev_m_rxn_ind, obj_ind = {}, {}, None
        for rxn_group in self.modified + [self.measured]: # Get LP indices of the reactions to be measured and modified.
            for r_id in rxn_group:
                rev_id = r_id + self._rxn_rev_suffix
                if r_id not in self.irr_model.reactions:
                    print('Could not find reaction "%s" in the model.' % r_id)
                    exit()
                """elif r_id in irr_m_rxn_ind:
                    print('Reaction "%s" was not unique in the model.' % r_id)
                    exit()"""
                irr_m_rxn_ind[r_id] = None
                rev_m_rxn_ind[r_id] = None
                if rev_id in self.irr_model.reactions:
                    irr_m_rxn_ind[rev_id] = None

        self._mtb_update = {}
        if mtb_update == None:
            mtb_update = (None, None, None)
        else:
            rxn_dict = self.rev_model.reactions.get_by_id(mtb_update[1]).metabolites
            orig_coef = rxn_dict[self.rev_model.metabolites.get_by_id(mtb_update[0])]
            self._mtb_update[mtb_update[2]] = {'rev_ind':None, 'irr_ind':None, 'mtb_rev_ind':None, 'mtb_irr_ind':None, 'mtb_orig_coef':orig_coef}
        for i, rxn in enumerate(self.rev_model.reactions):
            if rxn.id in rev_m_rxn_ind:
                rev_m_rxn_ind[rxn.id] = i
            if rxn.id == mtb_update[1]:
                self._mtb_update[mtb_update[2]]['rev_ind'] = i
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
            if rxn.id == mtb_update[1]:
                self._mtb_update[mtb_update[2]]['irr_ind'] = i
        for r_id, ind in irr_m_rxn_ind.items():
            if ind == None:
                print('"%s" in modify or measured did not have an index in the model' % r_id)
                exit()
        if obj_ind == None:
            print('Could not find an index for the objective function')
            exit()
        if self._mtb_update:
            for i, mtb in enumerate(self.rev_model.metabolites):
                if mtb.id == mtb_update[0]:
                    self._mtb_update[mtb_update[2]]['mtb_rev_ind'] = i
            for i, mtb in enumerate(self.irr_model.metabolites):
                if mtb.id == mtb_update[0]:
                    self._mtb_update[mtb_update[2]]['mtb_irr_ind'] = i
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
        modify_conditions = itertools.product(*self._steps) # All combinations of steps.
        data_inds = itertools.product( *[range(len(s)) for s in self._steps] ) # Where to store the result.
        for r_fs, ind in itertools.izip(modify_conditions, data_inds):
            labeled_conditions = [(r_id, rxn_flux) for rxn_group, rxn_flux in zip(self.modified, r_fs) for r_id in rxn_group]
            yield labeled_conditions, ind


class DmrVisualization(object):
    def __init__(self, colormap='RdBu', interpolation=False):
        """colormap should be one of the 'Diverging colormaps' from https://matplotlib.org/examples/color/colormaps_reference.html. interpolation is 'none' by default or 'spline36' if interpolation==True; a custom value can also be passed as the value of smoothed."""
        self.colormap = getattr(matplotlib.cm, colormap) if isinstance(colormap, basestring) else colormap
        self.interpolation = 'spline36' if interpolation==True else 'none' if interpolation==False else interpolation
    # # #  Graphing functions
    def heatmaps_2var(self, dmr, measured_rxns, include_objective=False, include_total_flux=False, graph_width=3.0, graph_height=2.5, per_row=6, shift_colorbar=True, show_all_x_axes=False, show_all_y_axes=False):
        if len(dmr.modified) != 2:
            print('Error: heatmaps_2var() can only be called with 2 dimensional data.')
            exit()
        measured_rxns = self._expand_measured_reactions(dmr, measured_rxns, include_objective, include_total_flux)
        num_graphs = len(measured_rxns)
        nrows, ncols = (num_graphs-1)//per_row + 1, min(num_graphs, per_row)
        figsize = (ncols*graph_width, nrows*graph_height)
        fig = plt.figure(figsize=figsize)
        fig.canvas.set_window_title(dmr.name)
        first_ax, min_flx, max_flx, cmap = None, None, None, self.colormap
        x_tick_inds, x_tick_vals = self._generate_axis_ticks(dmr.modified_attrs[1]['coefficient'], dmr._steps[1])
        y_tick_inds, y_tick_vals = self._generate_axis_ticks(dmr.modified_attrs[0]['coefficient'], dmr._steps[0])
        for ind, m_rxn in enumerate(measured_rxns):
            flx = dmr.get_fluxes(m_rxn)
            show_x_label, show_x_ticks, show_y_label, show_y_ticks = True, True, True, True
            if not show_all_x_axes and ind < num_graphs - per_row: # Applies to all graphs except bottom in each column.
                show_x_label, show_x_ticks = False, False
            if not show_all_y_axes and ind % per_row: # Applies to all graphs except left-most in each row.
                show_y_label, show_y_ticks = False, False
            if shift_colorbar:
                min_flx, max_flx = flx.min(), flx.max()
                cmap = self._shiftedColorMap(self.colormap, min_flx, max_flx, name='%s_colormap'%m_rxn)
            ax = fig.add_subplot(nrows, ncols, ind+1, sharex=first_ax, sharey=first_ax)
            #self._draw_heatmap(dmr, m_rxn, fig, ax, cmap, fluxes=flx, min_flux=min_flx, max_flux=max_flx)
            self._draw_heatmap(dmr, m_rxn, fig, ax, cmap, fluxes=flx, min_flux=min_flx, max_flux=max_flx, show_title=True, show_x_label=show_x_label, show_x_ticks=show_x_ticks, show_y_label=show_y_label, show_y_ticks=show_y_ticks)
            if ind == 0:
                first_ax = ax
                plt.sca(ax)
                plt.xticks(x_tick_inds, x_tick_vals) # These only need to be called...
                plt.yticks(y_tick_inds, y_tick_vals) # ... once for all shared graphs.
        plt.tight_layout()
        plt.show()

    def heatmaps_3var(self, dmr, measured_rxns, include_objective=False, include_total_flux=False, colorbar_ranges=None, graph_width=1.65, graph_height=1.65, shift_colorbar=True, shared_col_colorbar=True, show_all_titles=False, show_all_x_axes=False, show_all_y_axes=False, graph_colorbar_width_frac=None, row_label_width=None, row_squish=None, row_label_rpad=7, fig_top_padding=None):
        """The final 'modified_rxn' in self.modified is the one that defines the rows. It should have a relatively small number of steps. graph_colorbar_width_frac is a proportion of graph_width; so 0.25 means it is set to 25% of graph_width. row_label_width and row_squish do X. If not given a crude estimation will be used. row_label_rpad is the space between the label and the y-axis of the graph."""
        if len(dmr.modified) != 3:
            print('Error: heatmaps_3var() can only be called with 3 dimensional data.')
            exit()
        measured_rxns = self._expand_measured_reactions(dmr, measured_rxns, include_objective, include_total_flux)
        dim3_label = dmr.modified_attrs[2]['label']
        dim3_coef = dmr.modified_attrs[2]['coefficient']
        dim3_steps = [s*dim3_coef for s in dmr._steps[2]][::-1] # reversed so bottom row is first dim3 step, isntead of top row.
        nrows, ncols = len(dim3_steps), len(measured_rxns)
        num_graphs = nrows * ncols
        figsize, row_squish, fig_top_padding, wspace, hspace = self._calculate_3var_size_parameters(nrows, ncols, graph_width, graph_height, dim3_label, row_label_width, graph_colorbar_width_frac, row_squish, fig_top_padding)
        # Document what the size parameters do, inc wspace & hspace (should those be modifiable?)
        if shared_col_colorbar:
            colorbars = []
            if colorbar_ranges:
                for m_rxn, _range in zip(measured_rxns, colorbar_ranges):
                    min_flx, max_flx = _range
                    cmap = self._shiftedColorMap(self.colormap, min_flx, max_flx, name='%s_colormap'%m_rxn)
                    colorbars.append((cmap, min_flx, max_flx))
            else:
                for m_rxn in measured_rxns:
                    flx = dmr.get_fluxes(m_rxn)
                    min_flx, max_flx = flx.min(), flx.max()
                    cmap = self._shiftedColorMap(self.colormap, min_flx, max_flx, name='%s_colormap'%m_rxn)
                    colorbars.append((cmap, min_flx, max_flx))
        else:
            colorbars = [(self.colormap, None, None)] * ncols
        fig = plt.figure(figsize=figsize)
        fig.canvas.set_window_title(dmr.name)
        first_ax, min_flx, max_flx = None, None, None
        x_tick_inds, x_tick_vals = self._generate_axis_ticks(dmr.modified_attrs[1]['coefficient'], dmr._steps[1])
        y_tick_inds, y_tick_vals = self._generate_axis_ticks(dmr.modified_attrs[0]['coefficient'], dmr._steps[0])
        for ind, m_rxn in enumerate(measured_rxns * len(dim3_steps)):
            m_rxn_ind = ind % ncols
            dim3_ind = ind//len(measured_rxns)
            ax = fig.add_subplot(nrows, ncols, ind+1, sharex=first_ax, sharey=first_ax)
            flx = dmr.get_fluxes(m_rxn, [slice(None), slice(None), len(dim3_steps)-1-dim3_ind]) # the len-1-ind is needed so bottom row is first dim3 step, instead of top row.
            show_title, show_x_label, show_x_ticks, show_y_label, show_y_ticks = True, True, True, True, True
            if m_rxn_ind == 0:
                row_label = '%s\n[%.1f]' % (dim3_label, dim3_steps[dim3_ind])
                ax.annotate(row_label, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - row_label_rpad, 0), xycoords=ax.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
            if not show_all_titles and ind >= ncols: # Applies to all graphs except top row.
                show_title = False
            if not show_all_x_axes and ind < num_graphs - ncols: # Applies to all graphs except bottom in each column.
                show_x_label, show_x_ticks = False, False
            if not show_all_y_axes and ind % ncols: # Applies to all graphs except left-most in each row.
                show_y_label, show_y_ticks = False, False
            cmap, min_flx, max_flx = colorbars[m_rxn_ind]
            if not shared_col_colorbar and shift_colorbar:
                min_flx, max_flx = flx.min(), flx.max()
                cmap = self._shiftedColorMap(self.colormap, min_flx, max_flx, name='shiftedcmap_%i'%ind)
            self._draw_heatmap(dmr, m_rxn, fig, ax, cmap, fluxes=flx, min_flux=min_flx, max_flux=max_flx, show_title=show_title, show_x_label=show_x_label, show_x_ticks=show_x_ticks, show_y_label=show_y_label, show_y_ticks=show_y_ticks)
            if ind == 0:
                first_ax = ax
                plt.sca(ax)
                plt.xticks(x_tick_inds, x_tick_vals) # These only need to be called...
                plt.yticks(y_tick_inds, y_tick_vals) # ... once for all shared graphs.
        plt.tight_layout()
        fig.subplots_adjust(left=row_squish, top=1-fig_top_padding, wspace=wspace, hspace=hspace)
        plt.show()

    # # #  Main private graphing function
    def _draw_heatmap(self, dmr, measured_rxn, fig, axis, colormap, fluxes=None, min_flux=None, max_flux=None, show_title=True, show_x_label=True, show_x_ticks=True, show_y_label=True, show_y_ticks=True):
        min_flux_range = 5.0 # Affects the colouring range.
        colorbar_label_size = 8 # Default was 10.
        if fluxes is None:
            fluxes = dmr.get_fluxes(measured_rxn)
        if min_flux == max_flux == None:
            max_flux = max(abs(fluxes.max()), abs(fluxes.min()), min_flux_range)
            min_flux = -max_flux
        if show_title == True:
            title = dmr.measured_attrs[measured_rxn]['label']
        else:
            title = None
        colorbar_ticks = [int(ceil(fluxes.min())), int(floor(fluxes.max()))] # Annotates each subplot with its own min and max value ticks, even if drawn on the same scaled colorbar.
        # Zero is added to the colorbar ticks, unless one of the other values is close enough that the text would overlap.
        cbar_min_zero_dist = (colorbar_ticks[1] - colorbar_ticks[0]) / 10.0
        if abs(colorbar_ticks[0]) > cbar_min_zero_dist and abs(colorbar_ticks[1]) > cbar_min_zero_dist:
            colorbar_ticks.append(0)

        img = axis.imshow(fluxes, origin='lower', aspect='auto', vmin=min_flux, vmax=max_flux, cmap=colormap, interpolation=self.interpolation)
        cbar = fig.colorbar(img, ax=axis, shrink=0.94, aspect=15, pad=0.03, ticks=colorbar_ticks)
        cbar.ax.tick_params(labelsize=colorbar_label_size)
        self._format_graph(dmr, axis, title, show_x_label, show_x_ticks, show_y_label, show_y_ticks)

    def _format_graph(self, dmr, axis, title, show_x_label, show_x_ticks, show_y_label, show_y_ticks):
        x_attrs = dmr.modified_attrs[1]
        y_attrs = dmr.modified_attrs[0]
        if title != None:
            axis.set_title(title)
        if show_x_label:
            axis.set_xlabel(x_attrs['label'])
        if show_y_label:
            axis.set_ylabel(y_attrs['label'])
        if not show_x_ticks:
            axis.tick_params('x', which='both', bottom='off', labelbottom='off')
        if not show_y_ticks:
            axis.tick_params('y', which='both', left='off', labelleft='off')
        x_max_ind, y_max_ind = len(dmr._steps[1]) - 1, len(dmr._steps[0]) - 1
        def format_imshow_coords(x, y):
            # Defines the mouseover text.
            x_ind, y_ind = int(x + 0.5), int(y + 0.5)
            if x_ind > x_max_ind or y_ind > y_max_ind:
                return 'Index out of range'
            x_val, y_val = dmr._steps[1][x_ind], dmr._steps[0][y_ind]
            x_label, y_label = x_attrs['label'], y_attrs['label']
            return '%s=%.2f, %s=%.2f' % (x_label, x_val, y_label, y_val)
        axis.format_coord = format_imshow_coords # Sets the mouseover text.
    def _generate_axis_ticks(self, coef, steps):
        ax_vals = [round(s*coef, 1) for s in steps]
        ax_mid, ax_mid_ind = ((len(ax_vals)+1) // 2) - 1, (len(ax_vals) - 1) / 2.0
        ax_mid_val = round((ax_vals[ax_mid] + ax_vals[-ax_mid-1]) / 2.0, 1)
        return [0, ax_mid_ind, len(ax_vals)-1], [ax_vals[0], ax_mid_val, ax_vals[-1]]

    # # #  Misc private graphing functions
    def _expand_measured_reactions(self, dmr, measured_rxns, include_objective, include_total_flux):
        for m_rxn in measured_rxns:
            if m_rxn not in dmr.measured:
                print('Error: reaction "%s" was not measured in the given DiscreteModelResults.' % m_rxn)
                exit()
        if include_total_flux:
            measured_rxns = [dmr._total_flux_rxn_id] + measured_rxns
        if include_objective:
            measured_rxns = [dmr._obj_fxn_id] + measured_rxns
        return measured_rxns
    def _calculate_3var_size_parameters(self, nrows, ncols, graph_width, graph_height, dim3_label, row_label_width, graph_colorbar_width_frac, row_squish, fig_top_padding):
        if row_label_width == None:
            row_label_width = 0.50 + len(dim3_label) / 16.0
        if graph_colorbar_width_frac == None:
            graph_colorbar_width_frac = 0.27 + (ncols+1)**-1.85
        figsize = (ncols*graph_width*(1+graph_colorbar_width_frac)+row_label_width, nrows*graph_height)
        if row_squish == None:
            row_label_frac = row_label_width / figsize[0]
            row_squish = row_label_frac * 1.7
        if fig_top_padding == None:
            fig_top_padding = 0.4 / figsize[1]
        wspace, hspace = 0.20, 0.06
        return figsize, row_squish, fig_top_padding, wspace, hspace
    def _shiftedColorMap(self, cmap, min_val, max_val, name, epsilon=0.001):
        '''Function to offset the "center" of a colormap. Useful for data with a negative min and positive max and you want the middle of the colormap's dynamic range to be at zero. Adapted from https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

        Input
        -----
          cmap : The matplotlib colormap to be altered.
          start : Offset from lowest point in the colormap's range.
              Defaults to 0.0 (no lower ofset). Should be between
              0.0 and `midpoint`.
          midpoint : The new center of the colormap. Defaults to
              0.5 (no shift). Should be between 0.0 and 1.0. In
              general, this should be  1 - vmax/(vmax + abs(vmin))
              For example if your data range from -15.0 to +5.0 and
              you want the center of the colormap at 0.0, `midpoint`
              should be set to  1 - 5/(5 + 15)) or 0.75
          stop : Offset from highets point in the colormap's range.
              Defaults to 1.0 (no upper ofset). Should be between
              `midpoint` and 1.0.
        Note: If areas with value of zero are appearing as an extreme colour instead of white
          (or the middle colour), try increasing the value of epsilon.'''
        start, stop = 0.0, 1.0
        min_val, max_val = min(0.0, min_val), max(0.0, max_val) # Ensures 0 is included on the map.
        if max_val == min_val == 0:
            midpoint = 0.5
        else:
            midpoint = 1.0 - max_val/(max_val + abs(min_val))
            # don't think works if min and max are > 0. EDIT: What? Seems like it does.
        cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
        # regular index to compute the colors
        reg_index = np.linspace(start, stop, 257)
        # shifted index to match the data
        shift_index = np.hstack([np.linspace(0.0, midpoint, 128, endpoint=False), np.linspace(midpoint, 1.0, 129, endpoint=True)])
        for ri, si in zip(reg_index, shift_index):
            if abs(si - midpoint) < epsilon:
                r, g, b, a = cmap(0.5) # 0.5 = original midpoint.
            else:
                r, g, b, a = cmap(ri)
            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))
        newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
        plt.register_cmap(cmap=newcmap)
        return newcmap


def _flux_fva_to_rgb_frac(val, fva_range, min_val, max_val, max_fva_range, low_col, mid_col, high_col, var_col):
    """imshow function expects a set of float values, representing the 0-255 r/g/b integers.
    Uses transform on the fva_p which accelerates the colour towards var_col: fva_t=0.5 is a sqrt, 0.8 is less drastic, 0.3 is very drastic."""
    fva_t = 0.5
    if val == 0.0:
        val_rgb = mid_col
    else:
        if val > 0:
            val_p = val / float(max_val)
            delta_col = [float(c-m) for c,m in zip(high_col, mid_col)]
        else:
            val_p = abs(val) / abs(float(min_val))
            delta_col = [float(c-m) for c,m in zip(low_col, mid_col)]
        val_rgb = [d*val_p + c for d,c in zip(delta_col, mid_col)]
    fva_p = fva_range / float(max_fva_range)
    fva_p = fva_p**(fva_t)
    delta_var = [float(v-c) for v,c in zip(var_col, val_rgb)]
    rgb = [(d*fva_p + c)/255.0 for d,c in zip(delta_var, val_rgb)]
    return rgb
def _expand_steps(_min, _max, _steps):
    if _steps == 1:
        return [(_min+_max)/2.0]
    delta = (_max - _min) / float(_steps - 1)
    return [_min + delta*x for x in range(_steps)]

# TODO:
# - If the bounds on to_modify are illegal based on reaction bounds, get a ZeroDivisionError. Should be checked and reported more sanely.

if __name__ == '__main__':
    init_wol_percent = 10
    files_dir = '/mnt/hgfs/win_projects/brugia_project'
    #model_names = ['model_b_mal_4.5-wip.xlsx', 'model_b_mal_5_L3.xlsx', 'model_b_mal_5_L3D6.xlsx', 'model_b_mal_5_L3D9.xlsx', 'model_b_mal_5_L4.xlsx', 'model_b_mal_5_F30.xlsx', 'model_b_mal_5_M30.xlsx', 'model_b_mal_5_M30-overconstrained.xlsx']
    model_names = ['model_bm_5.xlsx']
    model_files = [os.path.join(files_dir, m_file) for m_file in model_names]
    #models = [read_excel(m_file, verbose=False) for m_file in model_files]
    models = [load_model(m_file, wol_ratio=init_wol_percent/100.0) for m_file in model_files]
    model = models[0]

    to_measure = {'M_TRANS_5':'Mitochondrial pyruvate', 'R00351_M':'Oaa -> Citrate', 'R01900_M':'Citrate -> Isocitrate', 'R00709_M':'Icit -> Akg (NADH)', 'R00267_M':'Icit -> Akg (NADPH)', 'R02570_M':'Akg -> Succinyl-CoA', 'R00405_M':'Suc-CoA -> Succinate', 'R02164_M':'Succinate -> Fumarate', 'R01082_M':'Fumarase', 'R00342_M':'Malate -> Oaa', 'R00086_M':'ATP synthase', 'RMC0184_M':'Rhod complex I', 'RMC0183_M':'Reverse complex II', 'R00479_M':'Glyoxylate pathway', 'SINK_2':'Succinate waste', 'SINK_3':'Acetate waste', 'SINK_4':'Propanoate waste', 'DIFFUSION_9':'H+ import', 'R05688_W':'Fosmidomycin target', 'R00127':'Tenofovir target', 'R00127_W':'Tenofovir Wolbachia', 'R00762':'MDL-29951 target'}
    negative_modified = [] # Must use the labels here, not the reaction IDs.
    negative_measured = ['R01082_M', 'R00086_M', 'R00351_M', 'R01900_M', 'R02570_M', 'R00405_M', 'R05688_W'] # Must use reaction IDs.
    # # #  Common reaction sets
    aa_ids = ['NUTRIENTS_%i'%i for i in range(1,21)]
    wol_scale_ids = []
    for rxn in model.reactions:
        if rxn.id.startswith('W_TRANS_') or rxn.id in ('BIO_NGAM_W',) or (rxn.id.startswith(('DIFFUSION_','CARBON_SOURCE_','FA_SOURCE_','SINK_','R')) and rxn.id.endswith('_W')):
            wol_scale_ids.append(rxn.id)
    tca_ids = ['R00351_M','R01900_M','R02570_M','R00405_M','R02164_M','R01082_M','R00342_M']
    drug_targets = ['R05688_W', 'R00127', 'R00127_W', 'R00762']
    # # #  Run options
    tight_bounds = False # Doesn't work with percent in to_modify
    modify_mtb_ratio_in_rxn_along = ('Bio_unscaled_w', 'BIOMASS_SCALED', 'DIFFUSION_1_W') # (mtb_id, in_rxn_id, along_rxn_id). Every time along_rxn_id has its bounds modified, the coefficient of mtb_id in in_rxn_id is also modified. along_rxn_id MUST be modified by 'percent', not absolute value.

    show_wol_ratio_graph = False
    show_2var_heatmap = False
    show_3var_heatmap = False
    show_knockout_heatmap = True
    show_3var_fva_heatmap = False


    if show_2var_heatmap:
        inc_total_flux = False
        inc_obj_fxn = True
        #to_modify = [('CARBON_SOURCE', 'Glucose', 1, 251, 25), ('DIFFUSION_2', 'Oxygen', 1.0, 701, 25)]
        to_modify = [('CARBON_SOURCE', 'Glucose', (1, 251, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25))]
        to_display = ['R01082_M', 'R00479_M', 'SINK_3', 'SINK_2', 'SINK_4']
        #to_modify = [('R05688_W', 'Fosmidomycin target', (-0.0, -1.5, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25))]
        #to_modify = [('R00127', 'Tenofovir target', (0, 60, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25))]
        #to_modify = [('R00762', 'MDL-29951 target', (0, 0.6, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25))]

        results = DiscreteModelResults(model, to_modify, to_measure, tight_bounds=tight_bounds)
        results.negative_modified(negative_modified)
        results.negative_measured(negative_measured)
        vis = DmrVisualization()
        vis.heatmaps_2var(results, to_display, include_objective=inc_obj_fxn, include_total_flux=inc_total_flux)
    if show_3var_heatmap:
        inc_total_flux = False
        inc_obj_fxn = True
        #to_modify = [('CARBON_SOURCE', 'Glucose', (1, 251, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25)), ('FA_SOURCE', 'Fatty acids', (25, 50, 6))]
        to_modify = [('CARBON_SOURCE', 'Glucose', (1, 251, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25)), ('FA_SOURCE', 'Fatty acids', (40, 50, 2))]
        #to_modify = [('CARBON_SOURCE', 'Glucose', (1, 251, 25)), ('DIFFUSION_2', 'Oxygen', (1.0, 701, 25)), (wol_scale_ids, 'wBm percent', (5, 20, 4, 'percent'))] # Set init_wol_percent to 100

        to_display = ['R01082_M', 'R00479_M', 'SINK_3', 'SINK_2', 'SINK_4']
        #to_display = ['R01082_M', 'SINK_3', 'SINK_2', 'SINK_4'] # Waste
        #to_display = drug_targets
        #to_display = tca_ids

        colorbar_ranges = None # If given, must be a list of tuples of floats, and match the order of to_display. If inc_obj_fxn==True, that must be the first range given.
        colorbar_ranges = [(0.0,58.0), (-332.0,213.0), (0.0,220.0), (0.0,296.0), (0.0,361.0), (0.0,412.0)]

        results = DiscreteModelResults(model, to_modify, to_measure, mtb_update=modify_mtb_ratio_in_rxn_along, tight_bounds=tight_bounds)
        results.negative_modified(negative_modified)
        results.negative_measured(negative_measured)
        vis = DmrVisualization()
        vis.heatmaps_3var(results, to_display, include_objective=inc_obj_fxn, include_total_flux=inc_total_flux, colorbar_ranges=colorbar_ranges)
    if show_3var_fva_heatmap:
        to_modify = [('CARBON_SOURCE', 'Glucose', 1, 251, 10), ('DIFFUSION_2', 'Oxygen', 1.0, 701, 10), ('FA_SOURCE', 'Fatty acids', 10, 35, 4)]
        to_display = ['M_TRANS_5', 'R01082_M', 'RMC0183_M', 'R00479_M', 'SINK_2', 'SINK_3', 'SINK_4']
        low_col, mid_col, high_col, var_col = [255,0,0], [255,255,255], [0,0,255], [0,0,0]
        min_val, max_val, max_fva_range = -1000, 1000, 2000
        num_val, num_fva = 21, 21
        m = np.zeros([num_val,num_fva,3]) # the 3 is for RGB
        for x in range(num_fva):
            fva_range = x / float(num_fva - 1) * max_fva_range
            for y, val in enumerate(_expand_steps(min_val, max_val, num_val)):
                rgba = _flux_fva_to_rgb_frac(val, fva_range, min_val, max_val, max_fva_range, low_col, mid_col, high_col, var_col)
                m[y,x] = rgba
        plt.imshow(m, origin='lower', extent=(0,max_fva_range,min_val,max_val))
        plt.show()
    if show_wol_ratio_graph:
        model_file = model_files[0]
        orig_model = load_model(model_file)
        wol_ratios = _expand_steps(0, 1.0, 101)
        to_modify = [('CARBON_SOURCE', 'Glucose', (45, 250, 2)), ('DIFFUSION_2', 'Oxygen', (90, 580, 2))]
        condition_names = ['LOLG', 'HOLG', 'LOHG', 'HOHG']
        ratio_vals = []
        for weight in wol_ratios:
            model = orig_model.copy()
            set_model_wolbachia(model, wol_ratio=weight)
            results = DiscreteModelResults(model, to_modify, to_measure, mtb_update=modify_mtb_ratio_in_rxn_along, tight_bounds=tight_bounds)
            ratio_vals.append(results.get_objective_fluxes().flatten())
        fig, ax = plt.subplots()
        if not condition_names:
            condition_names = [', '.join('{} ({:.1f})'.format(cond[1], val) for cond, val in zip(to_modify, conditions)) for conditions in itertools.product(*results._steps)]
        for name, obj_vals in zip(condition_names, zip(*ratio_vals)):
            ax.plot(wol_ratios, obj_vals, '-', label=name)
        ax.set(xlabel='$\it{Wolbachia}$ weight', ylabel='Objective function')
        ax.set_yticks([0, 10, 20, 30, 40])
        plt.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', fontsize=8, framealpha=0)
        plt.subplots_adjust(right=0.66)
        plt.show()
    if show_knockout_heatmap:
        from cobra.flux_analysis import single_reaction_deletion, double_reaction_deletion
        from scipy.cluster.hierarchy import linkage, leaves_list
        max_obj = float(model.optimize().f)
        threshold = max_obj / 2.0
        # Should I filter by only reactions with gene evidence? 849 vs 1011 to begin with.
        #to_knockout = [rxn for rxn in model.reactions if rxn.id.startswith(('ACYL', 'N000', 'R')) and rxn.gene_names]
        to_knockout = [rxn for rxn in model.reactions if rxn.id.startswith(('ACYL', 'N000', 'R'))]
        s_rates, _ = single_reaction_deletion(model, to_knockout)
        s_rxns = set(name for name, obj in s_rates.items() if obj < threshold)
        print('Considering {} reactions from the model...'.format(len(to_knockout)))
        print('{} single knockouts reduced biomass production below the threshold'.format(len(s_rxns)))
        d_to_knockout = [rxn for rxn in to_knockout if rxn.id not in s_rxns]
        double_ko_data = double_reaction_deletion(model, d_to_knockout, number_of_processes=3)
        d_r1, d_r2, d_rates = double_ko_data['y'], double_ko_data['x'], double_ko_data['data']
        hits = set()
        for i in range(len(d_r1)):
            for j in range(len(d_r2)):
                if d_rates[i,j] < threshold:
                    hits.add(i)
                    hits.add(j)
        print('{} double knockouts reduced biomass production below the threshold'.format(len(hits)))
        hits = sorted(hits)
        ko_data = d_rates[ np.ix_(hits, hits) ] # Selects only rows and columns in hits
        hit_order = [hits[ind] for ind in leaves_list(linkage(ko_data, method='ward'))]
        reaction_labels = [d_r1[ind] for ind in hit_order]
        graph_data = d_rates[ np.ix_(hit_order, hit_order) ] / max_obj

        # # #  Heatmap visual settings
        label_reactions = False
        title = None
        colourmap = 'inferno'
        colorbar_label_size = 8
        # # #  Setup the heatmap
        fig, ax = plt.subplots()
        im = ax.imshow(graph_data, cmap=colourmap, vmin=0, vmax=1, origin='lower')
        def format_imshow_coords(x, y):
            # Defines the mouseover text.
            x_val = int(floor(x + 0.5))
            y_val = int(floor(y + 0.5))
            if x_val < len(reaction_labels) and y_val < len(reaction_labels):
                return 'x={}, y={}'.format(reaction_labels[x_val], reaction_labels[y_val])
            else:
                return ''
        ax.format_coord = format_imshow_coords # Sets the mouseover text.
        cbar = fig.colorbar(im, ax=ax, ticks=[0, 0.5, 1.0])
        cbar.ax.tick_params(labelsize=colorbar_label_size)
        if title != None:
            ax.set_title(title)
        if label_reactions:
            ax.set_xticks(range(len(reaction_labels)))
            ax.set_yticks(np.arange(len(reaction_labels)))
            ax.set_xticklabels(reaction_labels)
            ax.set_yticklabels(reaction_labels)
            plt.xticks(fontsize=5, rotation=90)
            plt.yticks(fontsize=5)
        else:
            ax.tick_params('x', which='both', bottom='off', labelbottom='off')
            ax.tick_params('y', which='both', left='off', labelleft='off')
        fig.tight_layout()
        plt.show()
