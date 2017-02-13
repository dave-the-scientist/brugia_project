#xlrd is a module to read excel sheets.
#as is pandas.io.excel
"""Useful atts:
m.reactions or .metabolites: .get_by_id()
rxn.reactants or .products
rxn.get_coefficient('C00689'): -1 for consumed, +1 for produced
rxn.metabolites: {'C00689': -1, ...}
mtb.summary(): rates mtb is being produced and used in the current FBA.
model.metabolites.C00042 is the mtb object for succinate.
model.reactions.R02146.x is the flux through the reaction in the current FBA state.
"""
from read_excel import read_excel

# # #  Options
verbose = True # True or False.

def basic_stats(model):
    descrip_str = '\n\nBasic stats for model: %s' % model
    print('%s\n%s' % (descrip_str, '='*len(descrip_str.strip())))
    print('Contains: %i reactions, %i metabolites, and %i genes.' % (len(model.reactions), len(model.metabolites), len(model.genes)))
    topology_analysis(model)
    end_str = '\nEnd of basic stats for model %s' % model
    print('%s\n%s' % (end_str, '-'*len(end_str.strip())))
def topology_analysis(model):
    not_gen, not_used, blocked, sprfls = get_metab_outliers(model)
    print('\nBlocked reactions: %i allowing no flux.' % len(blocked))
    if len(sprfls) == 0:
        sprfls_mtb_str = '0.'
    else:
        sprfls_mtb_str = '%i described but not part of any reaction.' % len(sprfls)
    print('\nSuperfluous metabolites: %s' % sprfls_mtb_str)
    if len(not_gen) == 0:
        unavail_metab_str = '0. All metabolites produced by reactions.'
    else:
        unavail_metab_str = '%i metabolites used in reactions, but not produced by any.' % len(not_gen)
        if verbose:
            unavail_metab_str += '\n\t%s' % (' '.join(not_gen))
    print('\nUnavailable metabolites: %s' % unavail_metab_str)
    if len(not_used) == 0:
        unused_mtb_str = '0. All metabolites used in reactions.'
    else:
        unused_mtb_str = '%i metabolites produced in reactions, but not used by any.' % len(not_used)
        if verbose:
            unused_mtb_str += '\n\t%s' % (' '.join(not_used))
    print('\nUnused metabolites: %s' % unused_mtb_str)
def get_metab_outliers(model):
    mtbs, blocked = get_metabolite_reactions(model)
    sprfls = list(set(m.id for m in model.metabolites) - set(mtbs.keys()))
    not_gen, not_used = [], []
    for mtb, data in mtbs.items():
        if len(data['product']) == 0:
            not_gen.append(mtb)
        if len(data['reactant']) == 0:
            not_used.append(mtb)
    not_gen.sort(); not_used.sort(); blocked.sort(); sprfls.sort()
    return not_gen, not_used, blocked, sprfls
def get_metabolite_reactions(model):
    mtbs, blocked = {}, []
    for rxn in model.reactions:
        r_id = rxn.id
        forward = True if rxn.upper_bound > 0 else False
        reverse = True if rxn.lower_bound < 0 else False
        if forward == reverse == False:
            blocked.append(r_id)
            continue
        for mtb, coef in rxn.metabolites.items():
            m_id = mtb.id
            prd, rct = False, False
            if forward == reverse == True:
                prd, rct = True, True
            elif coef > 0: # is a product of this reaction.
                if forward:
                    prd = True
                elif reverse:
                    rct = True
            elif coef < 0: # is a reactant of this reaction.
                if reverse:
                    prd = True
                elif forward:
                    rct = True
            else:
                print('\nWarning: in reaction %s, metabolite %i has unknown coefficient "%s"' % (r_id, m_id, coef))
                continue
            if prd:
                mtbs.setdefault(m_id, {'product':[], 'reactant':[]})['product'].append(r_id)
            if rct:
                mtbs.setdefault(m_id, {'product':[], 'reactant':[]})['reactant'].append(r_id)
    return mtbs, blocked

def compare_models(m1, m2):
    descrip_str = '\n\nComparing models: %s and %s' % (m1, m2)
    print('%s\n%s' % (descrip_str, '='*len(descrip_str.strip())))
    compare_reactions(m1, m2)
    compare_metabolites(m1, m2)
    end_str = '\nEnd of comparison between %s and %s' % (m1, m2)
    print('%s\n%s' % (end_str, '-'*len(end_str.strip())))
def compare_reactions(m1, m2):
    m1_rxns = set(r.id for r in m1.reactions)
    m2_rxns = set(r.id for r in m2.reactions)
    common_rxns = m1_rxns & m2_rxns
    print('\nCommon reactions: %i.' % len(common_rxns))
    m1_unique = list(m1_rxns - m2_rxns)
    m1_uniq_str = '%i' % len(m1_unique)
    if verbose:
        m1_unique.sort()
        m1_uniq_str += ':\n\t%s' % (' '.join(m1_unique))
    print('\nReactions unique to %s: %s.' % (m1, m1_uniq_str))
    m2_unique = list(m2_rxns - m1_rxns)
    m2_uniq_str = '%i' % len(m2_unique)
    if verbose:
        m2_unique.sort()
        m2_uniq_str += ':\n\t%s' % (' '.join(m2_unique))
    print('\nReactions unique to %s: %s.' % (m2, m2_uniq_str))
def compare_metabolites(m1, m2):
    not_gen1, not_used1, _, _ = get_metab_outliers(m1)
    not_gen1, not_used1 = set(not_gen1), set(not_used1)
    not_gen2, not_used2, _, _ = get_metab_outliers(m2)
    not_gen2, not_used2 = set(not_gen2), set(not_used2)
    print_comparison(m1, m2, not_gen1, not_gen2, 'unavailable metabolites')
    print_comparison(m1, m2, not_used1, not_used2, 'unused metabolites')
def print_comparison(m1, m2, set1, set2, term):
    common = list(set1 & set2)
    unique1 = list(set1 - set2)
    unique2 = list(set2 - set2)
    if len(common) == 0:
        common_str = '0. No %s shared between the models' % term
    else:
        common_str = '%i %s shared between the models' % (len(common), term)
        if verbose:
            common.sort()
            common_str += ':\n\t%s' % (' '.join(common))
    print('\nCommon %s: %s.' % (term, common_str))
    unique_str1 = '%i' % len(unique1)
    if verbose and unique1:
        unique1.sort()
        unique_str1 += ':\n\t%s' % (' '.join(unique1))
    print('\n%s unique to %s: %s.' % (term.capitalize(), m1, unique_str1))
    unique_str2 = '%i' % len(unique2)
    if verbose and unique2:
        unique2.sort()
        unique_str2 += ':\n\t%s' % (' '.join(unique2))
    print('\n%s unique to %s: %s.' % (term.capitalize(), m2, unique_str2))


verbose = False
ov_model_file = 'model_o_vol.xlsx'
bm_model_file = 'model_b_mal.xlsx'

ov = read_excel(ov_model_file)
bm = read_excel(bm_model_file)

basic_stats(ov)
basic_stats(bm)
compare_models(ov, bm)

ov.optimize()
ov.summary()

bm.optimize()
bm.summary()
