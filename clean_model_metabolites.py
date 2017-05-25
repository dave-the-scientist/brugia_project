"""Run on the models_V3.5. Adds mtb compartments, ensures names are transferred, and renames the GXXXXX ids.
"""
import os
from read_excel import read_excel, write_excel

files_dir = '/mnt/hgfs/win_projects/brugia_project'
in_model_names = ['model_o_vol_3.5', 'model_b_mal_3.5']

in_model_files = [os.path.join(files_dir, m_file+'.xlsx') for m_file in in_model_names]
out_model_files = [os.path.join(files_dir, m_file+'-wip.xlsx') for m_file in in_model_names]

in_models = [read_excel(m_file, verbose=False) for m_file in in_model_files]

for m in in_models:
    for mtb in m.metabolites:
        if mtb.id == 'FA_storage_mix': continue
        elif mtb.id[0] == 'C':
            mtb.compartment = 'c'
            if mtb.id[1] == 'C': continue
        elif mtb.id[0] == 'M':
            mtb.compartment = 'm'
            if mtb.id[1] == 'C' or mtb.name: continue
            c_id = 'C' + mtb.id[1:]
            if c_id in m.metabolites and m.metabolites.get_by_id(c_id).name:
                mtb.name = m.metabolites.get_by_id(c_id).name
        elif mtb.id[0] == 'G':
            pass


for m, out_file in zip(in_models, out_model_files):
    write_excel(m, out_file)
