import os, itertools
from math import ceil, floor
import numpy as np
import cobra
from read_excel import read_excel
import matplotlib
import matplotlib.pyplot as plt


class ReactionVariances(object):
    """Holds results of testing the model by varying the reaction bounds of a number of reactions.

    self.data: numpy.Matrix. The main data structure is an (n+1)-dimensional matrix, where n is the number of modified reactions. Each position in the matrix holds a 1-dimensional array, whose 0-th entry is the objective function flux at that condition, the 1-th entry is the total flux in the entire system at that condition, and the following values are the flux through each reaction in self.measured, sorted by the name of the reaction ID.
    """
    def __init__(self, model, to_modify, to_measure, tight_bounds=False, solver=None):
        self.name = str(model) # Name of the model.
        self.modified = [] # List of reaction ids that were each modified and tested.
        self.modified_attrs = {} # Dict of properties used for graphing.
        self.measured = [] # List of reaction ids for which data are collected.
        self.measured_attrs = {} # List of human-readable descriptions.
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
        if isinstance(modified_rxn, basestring):
            if modified_rxn in self.modified_attrs:
                self.modified_attrs[modified_rxn]['coefficient'] = -1.0
        else:
            for m_rxn in modified_rxn:
                if m_rxn in self.modified_attrs:
                    self.modified_attrs[m_rxn]['coefficient'] = -1.0
    def negative_measured(self, measured_rxn):
        """Report negative values for this reaction."""
        if isinstance(measured_rxn, basestring):
            if measured_rxn in self.measured_attrs:
                self.measured_attrs[measured_rxn]['coefficient'] = -1.0
        else:
            for m_rxn in measured_rxn:
                if m_rxn in self.measured_attrs:
                    self.measured_attrs[m_rxn]['coefficient'] = -1.0
    def update_modified_attrs(self, modified_rxn, attr_dict):
        self.modified_attrs[modified_rxn].update(attr_dict)

    def heatmaps_2D(self, measured_rxns, include_objective=False, include_total_flux=False, graph_width=3.0, graph_height=2.5, per_row=6, colormap='RdBu', smoothed=False, show_all_x_axes=False, show_all_y_axes=False):
        if len(self.modified) != 2:
            print('Error: heatmaps_2D() can only be called with 2 dimensional data.')
            exit()
        measured_rxns, interpolation = self._initialize_heatmap_options(measured_rxns, include_objective, include_total_flux, smoothed)
        num_graphs = len(measured_rxns)
        nrows, ncols = (num_graphs-1)//per_row + 1, min(num_graphs, per_row)
        figsize = (ncols*graph_width, nrows*graph_height)
        fig = plt.figure(figsize=figsize)
        first_ax = None
        for ind, m_rxn in enumerate(measured_rxns):
            ax = fig.add_subplot(nrows, ncols, ind+1, sharex=first_ax, sharey=first_ax)
            self._draw_heatmap(m_rxn, fig, ax, colormap, interpolation)
            if ind == 0:
                first_ax = ax
            if not show_all_x_axes and ind < num_graphs - per_row: # Applies to all graphs except bottom in each column.
                ax.tick_params('x', which='both', bottom='off', labelbottom='off')
            if not show_all_y_axes and ind % per_row: # Applies to all graphs except left-most in each row.
                ax.tick_params('y', which='both', left='off', labelleft='off')
        plt.tight_layout()
        plt.show()
    def heatmaps_3D(self, measured_rxns, include_objective=False, include_total_flux=False, graph_width=1.65, graph_height=1.65, colormap='RdBu', smoothed=False, shift_colorbar=True, shared_col_colorbar=True, show_all_titles=False, show_all_x_axes=False, show_all_y_axes=False, graph_colorbar_width_frac=None, row_label_width=None, row_squish=None, row_label_rpad=10, fig_top_padding=None):
        """The final 'modified_rxn' in self.modified is the one that defines the rows. It should have a relatively small number of steps. colormap should be one of the 'Diverging colormaps' from https://matplotlib.org/examples/color/colormaps_reference.html. graph_colorbar_width_frac is a proportion of graph_width; so 0.25 means it is set to 25% of graph_width. row_label_width and row_squish do X. If not given a crude estimation will be used. row_label_rpad is the space between the label and the y-axis of the graph."""
        if len(self.modified) != 3:
            print('Error: heatmaps_3D() can only be called with 3 dimensional data.')
            exit()
        if isinstance(colormap, basestring):
            colormap = getattr(matplotlib.cm, colormap)
        measured_rxns, interpolation = self._initialize_heatmap_options(measured_rxns, include_objective, include_total_flux, smoothed)
        dim3_label = self.modified_attrs[self.modified[2]]['label']
        dim3_coef = self.modified_attrs[self.modified[2]]['coefficient']
        dim3_steps = [s*dim3_coef for s in self._steps[2]]
        nrows, ncols = len(dim3_steps), len(measured_rxns)
        if row_label_width == None:
            row_label_width = 0.48 + len(dim3_label) / 17.0
        if graph_colorbar_width_frac == None:
            graph_colorbar_width_frac = 0.27 + (ncols+1)**-1.85
        num_graphs = nrows * ncols
        figsize = (ncols*graph_width*(1+graph_colorbar_width_frac)+row_label_width, nrows*graph_height)

        row_label_frac = row_label_width / figsize[0]
        if row_squish == None:
            row_squish = row_label_frac * 1.7
        if fig_top_padding == None:
            fig_top_padding = 0.4 / figsize[1]
        wspace, hspace = 0.20, 0.06

        if shared_col_colorbar:
            colorbars = []
            for m_rxn in measured_rxns:
                flx = self._get_measured_fluxes(m_rxn)
                min_flx, max_flx = flx.min(), flx.max()
                #if m_rxn == '_objective_function': # TESTING
                #    min_flx = 1.0
                cmap = self._shiftedColorMap(colormap, min_flx, max_flx, name='%s_colormap'%m_rxn)
                colorbars.append((cmap, min_flx, max_flx))
        else:
            colorbars = [(colormap, None, None)] * ncols
        fig = plt.figure(figsize=figsize)
        first_ax, min_flx, max_flx = None, None, None
        for ind, m_rxn in enumerate(measured_rxns * len(dim3_steps)):
            dim3_ind = ind//len(measured_rxns)
            flx = self._get_measured_fluxes(m_rxn, [slice(None), slice(None), dim3_ind])
            m_rxn_ind = ind % ncols
            cmap, min_flx, max_flx = colorbars[m_rxn_ind]
            if not shared_col_colorbar and shift_colorbar:
                min_flx, max_flx = flx.min(), flx.max()
                cmap = self._shiftedColorMap(colormap, min_flx, max_flx, name='shiftedcmap_%i'%ind)
            ax = fig.add_subplot(nrows, ncols, ind+1, sharex=first_ax, sharey=first_ax)
            self._draw_heatmap(m_rxn, fig, ax, cmap, interpolation, fluxes=flx, min_flux=min_flx, max_flux=max_flx)
            if ind == 0:
                first_ax = ax
            if m_rxn_ind == 0:
                row_label = '%s\n[%.1f]' % (dim3_label, dim3_steps[dim3_ind])
                ax.annotate(row_label, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - row_label_rpad, 0), xycoords=ax.yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
            if not show_all_titles and ind >= ncols: # Applies to all graphs except top row.
                ax.axes.set_title('')
            if not show_all_x_axes and ind < num_graphs - ncols: # Applies to all graphs except bottom in each column.
                ax.axes.get_xaxis().set_visible(False)
                #ax.tick_params('x', which='both', bottom='off', labelbottom='off') # Hides the axis values, but keeps the axis title.
            if not show_all_y_axes and ind % ncols: # Applies to all graphs except left-most in each row.
                ax.axes.get_yaxis().set_visible(False)
                #ax.tick_params('y', which='both', left='off', labelleft='off') # Hides the axis values, but keeps the axis title.
        plt.tight_layout()
        fig.subplots_adjust(left=row_squish, top=1-fig_top_padding, wspace=wspace, hspace=hspace)
        plt.show()
    def _initialize_heatmap_options(self, measured_rxns, include_objective, include_total_flux, smoothed):
        for m_rxn in measured_rxns:
            if m_rxn not in self.measured:
                print('Error: the given reaction "%s" was not measured in the current data.' % m_rxn)
                exit()
        if include_total_flux:
            measured_rxns = [self._total_flux_rxn_id] + measured_rxns
        if include_objective:
            measured_rxns = [self._obj_fxn_id] + measured_rxns
        interpolation = 'spline36' if smoothed else 'none'
        return measured_rxns, interpolation
    def _draw_heatmap(self, measured_rxn, fig, axis, colormap, interpolation, fluxes=None, min_flux=None, max_flux=None):
        min_flux_range = 5.0 # Affects the colouring range.
        colorbar_label_size = 8 # Default was 10.
        if fluxes is None:
            fluxes = self._get_measured_fluxes(measured_rxn, subslice)
        if min_flux == max_flux == None:
            max_flux = max(abs(fluxes.max()), abs(fluxes.min()), min_flux_range)
            min_flux = -max_flux
        colorbar_ticks = [int(ceil(fluxes.min())), 0, int(floor(fluxes.max()))] # Annotates each subplot with its own min and max value ticks, even if drawn on the same scaled colorbar.
        #colorbar_ticks = [int(min_flux), 0, int(max_flux)] # This if you want colorbars to all have the same ticks with shared_col_colorbar=True.
        img = axis.imshow(fluxes, origin='lower', aspect='auto', vmin=min_flux, vmax=max_flux, cmap=colormap, interpolation=interpolation)
        cbar = fig.colorbar(img, ax=axis, shrink=0.94, aspect=15, pad=0.03, ticks=colorbar_ticks)
        cbar.ax.tick_params(labelsize=colorbar_label_size)
        x_attrs = self.modified_attrs[self.modified[1]]
        y_attrs = self.modified_attrs[self.modified[0]]
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
    def _get_measured_fluxes(self, measured_rxn, subslice=None):
        if measured_rxn == self._obj_fxn_id:
            fluxes = self.get_objective_flux() * self.measured_attrs[measured_rxn]['coefficient']
        elif measured_rxn == self._total_flux_rxn_id:
            fluxes = self.get_total_flux() * self.measured_attrs[measured_rxn]['coefficient']
        else:
            fluxes = self.get_flux(measured_rxn) * self.measured_attrs[measured_rxn]['coefficient']
        if subslice:
            fluxes = fluxes[subslice]
        return fluxes


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
        if not self.tight_bounds:
            fwd_f_lb, rev_f_lb = 0.0, 0.0
        for rxn_id, rxn_f in rxn_steps:
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
                rxn_lb, rxn_ub = min(rxn_f, 0.0), max(rxn_f+self.epsilon, 0.0)
            self.rev_lp.change_variable_bounds(self._rev_m_rxn_ind[rxn_id], rxn_lb, rxn_ub)
            self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rxn_id], fwd_f_lb, fwd_f)
            rev_id = rxn_id + self._rxn_rev_suffix
            if rev_id in self._irr_m_rxn_ind:
                self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rev_id], rev_f_lb, rev_f)
            """
            self.rev_lp.change_variable_bounds(self._rev_m_rxn_ind[rxn_id], rxn_f, rxn_f+self.epsilon)
            self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rxn_id], fwd_f, fwd_f+self.epsilon)
            rev_id = rxn_id + self._rxn_rev_suffix
            if rev_id in self._irr_m_rxn_ind:
                self.irr_lp.change_variable_bounds(self._irr_m_rxn_ind[rev_id], rev_f, rev_f+self.epsilon)
            """

    def _measure_reactions(self, data_ind, obj_f, total_f):
        self.data[data_ind][0] = obj_f
        self.data[data_ind][1] = total_f
        for i, rxn_id in enumerate(self.measured):
            if obj_f == 0:
                val = 0
            else:
                val = self.irr_model.reactions.get_by_id(rxn_id).x
                rev_id = rxn_id + self._rxn_rev_suffix
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
        default_modified_attrs = {'coefficient':1.0}
        default_measured_attrs = {'coefficient':1.0}
        if not to_modify or not to_measure:
            print('At least one reaction must be given to modify, and at least one to measure')
            exit()
        self.modified = [m[0] for m in to_modify] # Retains the same order.
        self.measured = list(sorted(to_measure.keys())) # Sorted by reaction id.
        while self._obj_fxn_id in self.measured:
            self._obj_fxn_id = '_' + self._obj_fxn_id
        while self._total_flux_rxn_id in self.measured:
            self._total_flux_rxn_id = '_' + self._total_flux_rxn_id
        self.modified_attrs, self.measured_attrs = {}, {}
        for m in to_modify:
            self.modified_attrs.setdefault(m[0], {'label':m[1]}).update(default_modified_attrs)
        for m, n in to_measure.items():
            self.measured_attrs.setdefault(m, {'label':n}).update(default_measured_attrs)
        self.measured_attrs.setdefault(self._obj_fxn_id, {'label':self._obj_fxn_label}).update(default_measured_attrs)
        self.measured_attrs.setdefault(self._total_flux_rxn_id, {'label':self._total_flux_rxn_label}).update(default_measured_attrs)
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
            rev_id = r_id + self._rxn_rev_suffix
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

    # # # # #  Private methods: graphing. # # # # #
    def _shiftedColorMap(self, cmap, min_val, max_val, name):
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
              `midpoint` and 1.0.'''
        epsilon = 0.0001
        start, stop = 0.0, 1.0
        min_val, max_val = min(0.0, min_val), max(0.0, max_val)
        midpoint = 1.0 - max_val/(max_val + abs(min_val))


        # don't think works if min and max are > 0.

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

    to_measure = {'M_TRANS_5':'Mitochondrial pyruvate', 'R01082_M':'Fumarate -> malate', 'R00086_M':'ATP synthase', 'RMC0184_M':'Rhod complex I', 'RMC0183_M':'Reverse complex II', 'R00479_M':'Glyoxylate pathway', 'SINK_2':'Succinate waste', 'SINK_3':'Acetate waste', 'SINK_4':'Propanoate waste'}
    negative_modified = 'DIFFUSION_2'
    negative_measured = ['R01082_M', 'R00086_M']
    tight_bounds = False

    show_2D_heatmap = False
    show_3D_heatmap = True

    model_files = [os.path.join(files_dir, m_file) for m_file in model_names]
    models = [read_excel(m_file, verbose=False) for m_file in model_files]

    if show_2D_heatmap:
        to_modify = [('CARBON_SOURCE', 'Glucose', 0, 200, 20), ('DIFFUSION_2', 'Oxygen', -0.0, -500, 20)]
        to_display = ['M_TRANS_5', 'R01082_M', 'R00086_M', 'RMC0183_M', 'R00479_M']
        rv = ReactionVariances(models[0], to_modify, to_measure, tight_bounds=tight_bounds)
        rv.negative_modified(negative_modified)
        rv.negative_measured(negative_measured)
        rv.heatmaps_2D(to_display, include_objective=True, include_total_flux=True)
    if show_3D_heatmap:
        to_modify = [('CARBON_SOURCE', 'Glucose', 1, 251, 25), ('DIFFUSION_2', 'Oxygen', -1.0, -501, 25), ('FA_SOURCE', 'Fatty acids', 5, 25, 6)]
        to_display = ['M_TRANS_5', 'R01082_M', 'RMC0183_M', 'R00479_M', 'SINK_2', 'SINK_3', 'SINK_4']
        rv = ReactionVariances(models[0], to_modify, to_measure, tight_bounds=tight_bounds)
        rv.negative_modified(negative_modified)
        rv.negative_measured(negative_measured)
        rv.heatmaps_3D(to_display, include_objective=True, include_total_flux=False)
