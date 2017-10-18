import os
from cobra.flux_analysis import single_reaction_deletion, double_reaction_deletion
from model_tools import id_bottleneck_metabolites
from read_excel import read_excel

def get_rxns_to_delete(model):
    rxn_to_genes = {}
    for rxn in model.reactions:
        if not rxn.gene_names or not rxn.id.startswith(('R', 'ACYLCOA')):
            continue
        rxn_to_genes[rxn.id] = rxn.gene_names.split(';')
        #if len(rxn_to_genes) >= 10: break
    return rxn_to_genes

def do_deletions(model, reactions, obj_fraction=0.0):
    fraction_epsilon = 0.0001
    rxn_deficiencies = []
    orig_f = float(model.optimize().f)
    s_rates, s_stats = single_reaction_deletion(model, reactions)
    print('%i reactions knocked out; original objective %.1f' % (len(s_stats), orig_f))
    for r_id, new_f in s_rates.items():
        if abs(new_f) < fraction_epsilon:
            new_f = 0.0
        stat = s_stats[r_id]
        if new_f/orig_f <= obj_fraction+fraction_epsilon:
            if stat == 'optimal':
                deficiencies = find_model_deficiencies(model, orig_f, new_f, r_id)
            else:
                deficiencies = 'infeasible'
            rxn_deficiencies.append( ('%s (%.1f)' % (r_id, new_f/orig_f*100), deficiencies) )
    rxn_deficiencies.sort()
    rxn_deficiencies.sort(key=lambda d: d[1])
    print('%i reactions with significant impact:' % len(rxn_deficiencies))
    for r_info, defic in rxn_deficiencies:
        print('%s\n\t%s' % (r_info, defic))
    #d_rates, d_stats = double_reaction_deletion(model, reactions)

def find_model_deficiencies(model, orig_f, new_f, r_id):
    deficiencies = []
    ob = model.reactions.get_by_id(r_id).bounds
    model.reactions.get_by_id(r_id).bounds = (0,0)
    diffs = id_bottleneck_metabolites(model, new_f, 'BIOMASS', threshold=1.0)
    for recovered_f, mtb_id in diffs:
        def_str = '%s (%.1f)' % (mtb_id, recovered_f/orig_f*100)
        sub_defs = []
        for sub_f, sub_mtb_id in id_bottleneck_metabolites(model, new_f, mtb_id.upper(), threshold=1.0):
            sub_defs.append('%s[%.1f]' % (sub_mtb_id, sub_f/orig_f*100))
        if sub_defs:
            def_str += ': %s' % ', '.join(sub_defs)
        deficiencies.append(def_str)
    model.reactions.get_by_id(r_id).bounds = ob
    if not deficiencies:
        return 'unrecoverable'
    else:
        return ', '.join(deficiencies)

# # #  Options
files_dir = '/mnt/hgfs/win_projects/brugia_project'
model_file = 'model_b_mal_4.5-wip.xlsx'
objective_threshold_fraction = 0.1 # Considered significant if resulting objective flux is less than 0.1 (10%) of the original.

# # #  Run steps
model_path = os.path.join(files_dir, model_file)
model = read_excel(model_path, verbose=False)
rxn_to_genes = get_rxns_to_delete(model)
do_deletions(model, list(rxn_to_genes.keys()), objective_threshold_fraction)
