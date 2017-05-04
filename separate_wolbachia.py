from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.flux_analysis import flux_variability_analysis
from model_tools import test_changed_constraints
from read_excel import read_excel, write_excel
import os

def separate_wol_rxns(in_models, wol_model_names, dont_remove_wol_rxns=False):
    """Modifies the in_models in place."""
    wol_gene_prefs = ('Wolbachia', 'AAW')
    wol_models, unique_wol_ids = [], []
    for m, wol_name in zip(in_models, wol_model_names):
        nem_mod_rxns, nem_deleted_rxns = [], []
        wm = Model(wol_name)
        for rxn in m.reactions:
            if not rxn.gene_names: continue
            wol_gene_names, nem_gene_names = [], []
            for g_name in rxn.gene_names.split(';'):
                g_name = g_name.strip()
                if g_name.startswith(wol_gene_prefs):
                    wol_gene_names.append(g_name)
                else:
                    nem_gene_names.append(g_name)
            wol_gene_names = ';'.join(wol_gene_names)
            nem_gene_names = ';'.join(nem_gene_names)
            if wol_gene_names:
                r_id = rxn.id
                if r_id.endswith('_M'):
                    r_id = r_id[:-2]
                w_id = r_id + '_W'
                if w_id not in wm.reactions:
                    w_rxn = Reaction(w_id)
                    setup_wol_reaction(rxn, w_rxn, wol_gene_names)
                    wm.add_reaction(w_rxn)
                    # add attributes.
                if nem_gene_names: # rxn found in nematode and wolbachia. remove wolbachia gene names.
                    nem_mod_rxns.append(rxn.id)
                    rxn.gene_names = nem_gene_names
                elif rxn.id not in nem_dont_delete: # rxn found only in wolbachia. remove it completely.
                    nem_deleted_rxns.append(rxn.id)
                    if dont_remove_wol_rxns == False:
                        rxn.delete()
            else: # rxn had no nematode or wolbachia gene names.
                pass
        wol_models.append(wm)
        if dont_remove_wol_rxns:
            unique_wol_ids.append(nem_deleted_rxns)
    return wol_models, unique_wol_ids
def setup_wol_reaction(n_rxn, w_rxn, wol_gene_names):
    mtb_dict = {}
    n_cmp = ''
    for mtb, coef in n_rxn.metabolites.items():
        if not n_cmp:
            n_cmp = mtb.id[0]
        else:
            if mtb.id[0] != n_cmp:
                print n_rxn.id, n_rxn.reaction
        w_mtb_id = 'W%s' % mtb.id[1:]
        w_mtb = Metabolite(w_mtb_id, name=mtb.name, compartment='w')
        mtb_dict[w_mtb] = coef
    w_rxn.gene_names = wol_gene_names
    w_rxn.bounds = n_rxn.bounds
    w_rxn.name = n_rxn.name
    w_rxn.objective_coefficient = n_rxn.objective_coefficient
    w_rxn.subsystem = n_rxn.subsystem
    w_rxn.enzyme_commission = n_rxn.enzyme_commission
    w_rxn.protein_names = n_rxn.protein_names
    w_rxn.confidence_notes = n_rxn.confidence_notes
    w_rxn.reaction_notes = n_rxn.reaction_notes
    w_rxn.gene_reaction_rule = n_rxn.gene_reaction_rule



# # #  Inputs
files_dir = '/mnt/hgfs/win_projects/brugia_project'
in_model_files = ['model_o_vol_3.5.xlsx', 'model_b_mal_3.5.xlsx']
wol_model_names = ['model_wOv_1', 'model_wBm_1']
out_nem_model_names = ['model_o_vol_4', 'model_b_mal_4']
nem_dont_delete = set(['R00104', 'R00161'])
# # #  Run-time options
save_wol_models = False
save_nem_models = True
test_removed_wol_rxns = True # Check the effect of removing the wolbachia-only reactions

# # #  Run steps
in_model_files = [os.path.join(files_dir, m_file) for m_file in in_model_files]
in_models = [read_excel(m_file, verbose=False) for m_file in in_model_files]
fvas = []
for m in in_models:
    optimize_minimal_flux(m)
    fvas.append(flux_variability_analysis(m))
wol_models, unique_wol_ids = separate_wol_rxns(in_models, wol_model_names, test_removed_wol_rxns)

if test_removed_wol_rxns:
    for m, fva, r_ids in zip(in_models, fvas, unique_wol_ids):
        print('Evaluating %i reactions from %s. Initial flux: %.1f.' % (len(r_ids), m, m.solution.f))
        diffs = test_changed_constraints(r_ids, m, fva, bounds_deltas=None, cumulative=False)
        # Helpful if tried killing each of the reactions, then adding free imports of each biomass category back in, to see which is most affected.
        for d in diffs:
            print d

if save_wol_models:
    out_wol_model_files = [os.path.join(files_dir, '%s-wip.xlsx'%m_name) for m_name in wol_model_names]
    for wm, m_file in zip(wol_models, out_wol_model_files):
        write_excel(wm, m_file)
if save_nem_models:
    out_nem_model_files = [os.path.join(files_dir, '%s-wip.xlsx'%m_name) for m_name in out_nem_model_names]
    for m, m_file in zip(in_models, out_nem_model_files):
        write_excel(m, m_file)