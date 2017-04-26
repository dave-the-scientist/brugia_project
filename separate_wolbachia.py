from cobra import Model, Reaction, Metabolite
from read_excel import read_excel, write_excel
import os

def separate_wol_rxns(in_models, wol_model_names):
    """Modifies the in_models in place."""
    wol_gene_prefs = ('Wolbachia', 'AAW')
    wol_models = []
    for m, wol_name in zip(in_models, wol_model_names):
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
                if r_id not in wm.reactions:
                    w_id = r_id + '_W'
                    w_rxn = Reaction(w_id)
                    # add attributes.
                if nem_gene_names: # rxn found in nematode and wolbachia. remove wolbachia gene names.
                    rxn.gene_names = nem_gene_names
                else: # rxn found only in wolbachia. remove it completely.
                    rxn.delete()
            else: # rxn had no nematode or wolbachia gene names.
                pass
        wol_models.append(wm)
    return wol_models

# # #  Inputs
files_dir = '/mnt/hgfs/win_projects/brugia_project'
in_model_files = ['model_o_vol_3.5.xlsx', 'model_b_mal_3.5.xlsx']
wol_model_names = ['model_wOv_1', 'model_wBm_1']
out_nem_model_names = ['model_o_vol_4', 'model_b_mal_4']
# # #  Run-time options
save_wol_models = False
save_nem_models = False

# # #  Run steps
in_model_files = [os.path.join(files_dir, m_file) for m_file in in_model_files]
in_models = [read_excel(m_file, verbose=False) for m_file in in_model_files]
wol_models = separate_wol_rxns(in_models, wol_model_names)

if save_wol_models:
    out_wol_model_files = [os.path.join(files_dir, '%s-wip.xlsx'%m_name) for m_name in wol_model_names]
    for wm, m_file in zip(wol_models, out_wol_model_files):
        write_excel(wm, m_file)
if save_nem_models:
    out_nem_model_files = [os.path.join(files_dir, '%s-wip.xlsx'%m_name) for m_name in out_nem_model_names]
    for m, m_file in zip(wol_models, out_nem_model_files):
        write_excel(m, m_file)
