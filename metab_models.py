#xlrd is a module to read excel sheets.
#as is pandas.io.excel
"""Useful atts:
m.reactions or .metabolites: .get_by_id()
rxn.reactants or .products
rxn.get_coefficient('C00689'): -1 for consumed, +1 for produced
rxn.metabolites: {'C00689': -1, ...}
mtb.summary(): rates mtb is being produced and used in the current FBA.
"""
from read_excel import read_excel

def basic_stats(model):
    print('The model %s consists of %i reactions, %i metabolites, and %i genes.' % (model, len(model.reactions), len(model.metabolites), len(model.genes)))
def find_orphans(model):
    mtbs = {}
    for rxn in model.reactions:
        r_id = rxn.id
        for mtb, coef in rxn.metabolites.items():
            if coef > 0:
                # produced.
                pass


ov_model_file = 'model_o_vol.xlsx'
bm_model_file = 'model_b_mal.xlsx'

ov = read_excel(ov_model_file)
bm = read_excel(bm_model_file)

basic_stats(ov)
basic_stats(bm)

find_orphans(bm)
