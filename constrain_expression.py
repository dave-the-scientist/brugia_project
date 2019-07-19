"""
Assumptions:
- Gene names in the expression file are of the form "Bm7142".
  - They have some prefix followed by digits. They do not end in a letter, and they are not subdivided like "Bm7142.3".
"""

import os, cobra, pandas
from read_excel import read_excel, write_excel

def extract_gene_reactions(model, worm_gene_prefix, wol_gene_prefix):
    """Returns 2 dicts of the form {'SeqName1':set(reaction1, reaction3), ...}"""
    worm_genes, wol_genes = {}, {}
    with_names = 0
    for rxn in model.reactions:
        if not rxn.gene_names:
            continue
        with_names += 1
        for name in rxn.gene_names.split(';'):
            name = name.strip()
            if '.' in name:
                name = name.partition('.')[0]
            while not name[-1].isdigit():
                name = name[:-1]
            if name.startswith(worm_gene_prefix):
                worm_genes.setdefault(name, set()).add(rxn.id)
            elif name.startswith(wol_gene_prefix):
                wol_genes.setdefault(name, set()).add(rxn.id)
            else:
                print 'Unknown gene name found: %s' % name
    print('Of %i reactions in %s, %i had gene names.' % (len(model.reactions), model, with_names))
    print('%i unique genes from the worm, %i from the Wolbachia.' % (len(worm_genes), len(wol_genes)))
    return worm_genes, wol_genes

def parse_expression_file(filename, sheetname, gene_rxns, conditions, replicate_inds):
    """Returns a dict {'rxn1':[0.3, 0.99, 1.0,... condition_n], ...}.
    """
    rxn_consts = {}
    seq_name_key = 'Sequence Name'
    frame = pandas.read_excel(filename, sheetname)
    if len(frame.columns) != len(set(frame.columns)):
        print('Error: at least one column header was not unique in sheet %s.' % sheetname)
        exit()
    cond_keys = [[cond+ind for ind in replicate_inds if cond+ind in frame.columns] for cond in conditions]
    for ck, cond in zip(cond_keys, conditions):
        if len(ck) == 0:
            print('Error: no replicates found for condition "%s" in sheet %s.' % (cond, sheetname))
            exit()
    for i in frame.index:
        row = frame.ix[i]
        seq_name = row[seq_name_key]
        if seq_name not in gene_rxns:
            continue
        avgs = [sum(row[k] for k in ck)/float(len(ck)) for ck in cond_keys]
        # modify consts if I want more/less severe penalties
        if any(avgs):
            consts = [a/max(avgs) for a in avgs]
        else:
            consts = [0.0 for cond in conditions]
        for rxn in gene_rxns[seq_name]:
            rxn_consts.setdefault(rxn, []).append(consts)
    rxn_constraints = {}
    for rxn, const_list in rxn_consts.items():
        if len(const_list) == 1:
            rxn_constraints[rxn] = dict(zip(conditions, const_list[0]))
        else:
            # This is not the best way. Should deal with AND vs OR for genes; this is a compromise.
            consts = [sum(stage_consts)/float(len(const_list)) for stage_consts in zip(*const_list)]
            rxn_constraints[rxn] = dict(zip(conditions, consts))
    print('Parsed constraints for %i reactions from %s.' % (len(rxn_constraints), sheetname))
    return rxn_constraints

def save_stage_models(model, worm_constraints, wol_constraints, conditions, model_out_str, out_dir, ignore_zeros=True):
    """If ignore_zeros is True, an expression level of zero will not modify a reaction. This is done because it can be easy to miss a transcript in an experiment."""
    for cond in conditions:
        model_name = model_out_str % cond
        m = model.copy()
        m.id = model_name
        for rxn, consts in worm_constraints.items():
            if ignore_zeros and consts[cond] == 0:
                continue
            bounds = m.reactions.get_by_id(rxn).bounds
            bounds = (bounds[0]*consts[cond], bounds[1]*consts[cond])
            m.reactions.get_by_id(rxn).bounds = bounds
        for rxn, consts in wol_constraints.items():
            if ignore_zeros and consts[cond] == 0:
                continue
            bounds = m.reactions.get_by_id(rxn).bounds
            bounds = (bounds[0]*consts[cond], bounds[1]*consts[cond])
            m.reactions.get_by_id(rxn).bounds = bounds
        outfile = os.path.join(out_dir, model_name+'.xlsx')
        write_excel(m, outfile)
        print('Wrote new model to %s' % outfile)


# # #  Inputs
files_dir = '/mnt/hgfs/win_projects/brugia_project'
model_file = 'model_bm_5.xlsx'
expression_file = 'All_Stages_Brugia_Wolbachia_FPKMs.xlsx'
worm_sheet = 'Brugia_FPKMs'
wol_sheet = 'Wolbachia_FPKMs'
# # #  Options
model_out_str = 'model_bm_5_%s'
conditions = ['L3', 'L3D6', 'L3D9', 'L4', 'F30', 'M30', 'F42', 'M42', 'F120', 'M120']
worm_gene_prefix, wol_gene_prefix = 'Bm', 'AAW'
replicate_inds = ['a', 'b', 'c']


# # #  Run steps
model_path = os.path.join(files_dir, model_file)
exp_path = os.path.join(files_dir, expression_file)
model = read_excel(model_path, verbose=False)
worm_genes, wol_genes = extract_gene_reactions(model, worm_gene_prefix, wol_gene_prefix)

worm_constraints = parse_expression_file(exp_path, worm_sheet, worm_genes, conditions, replicate_inds)
wol_constraints = parse_expression_file(exp_path, wol_sheet, wol_genes, conditions, replicate_inds)

save_stage_models(model, worm_constraints, wol_constraints, conditions, model_out_str, files_dir)
