from read_excel import read_excel, write_excel
import os

def merge_models(worm_m, wol_m, out_name, out_file):
    rxn_pref_remove = ('TEST_',)
    worm_m += wol_m
    rxns_to_remove = []
    for rxn in worm_m.reactions:
        if rxn.id.startswith(rxn_pref_remove):
            rxns_to_remove.append(rxn.id)
    worm_m.remove_reactions(rxns_to_remove)
    worm_m.id = out_name
    write_excel(worm_m, out_file)
    print('%s saved to %s' % (worm_m, out_file))



# # #  Inputs
files_dir = '/mnt/hgfs/win_projects/brugia_project'
worm_model_files = ['model_b_mal_4-sep.xlsx', 'model_o_vol_4-sep.xlsx']
wol_model_files = ['model_wBm_4-sep.xlsx', 'model_wOv_4-sep.xlsx']
out_model_names = ['model_b_mal_4.5', 'model_o_vol_4.5']
# # #  Run-time options


# # #  Run steps
worm_model_files = [os.path.join(files_dir, m_file) for m_file in worm_model_files]
wol_model_files = [os.path.join(files_dir, m_file) for m_file in wol_model_files]
worm_models = [read_excel(m_file, verbose=False) for m_file in worm_model_files]
wol_models = [read_excel(m_file, verbose=False) for m_file in wol_model_files]

for worm_m, wol_m, out_name in zip(worm_models, wol_models, out_model_names):
    out_file = os.path.join(files_dir, out_name+'-wip.xlsx')
    merge_models(worm_m, wol_m, out_name, out_file)
