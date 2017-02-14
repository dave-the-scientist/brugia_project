#xlrd is a module to read excel sheets.
#as is pandas.io.excel
"""Useful atts:
m.reactions or .metabolites: .get_by_id()
rxn.reactants or .products
rxn.get_coefficient('C00689'): -1 for consumed, +1 for produced
rxn.metabolites: {'C00689': -1, ...}
mtb.reactions: all of the reactions that produce or consume mtb.
mtb.summary(): rates mtb is being produced and used in the current FBA.
model.metabolites.C00042 is the mtb object for succinate.
model.reactions.R02146.x is the flux through the reaction in the current FBA state.
model.reactions.R02146.reaction is a human-readable description.
"""
from math import sqrt
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

def compare_objective_functions(rxn1, rxn2):
    m1, m2 = rxn1.model, rxn2.model
    coefs1, coefs2 = rxn1.metabolites, rxn2.metabolites
    mtbs1 = set(c.id for c in coefs1)
    mtbs2 = set(c.id for c in coefs2)
    common = list(mtbs1 & mtbs2)
    common.sort(key=lambda c: rxn1.get_coefficient(c))
    print('\nObjective metabolites in common: %i.' % len(common))
    buff = []
    for c in common:
        mtb1 = m1.metabolites.get_by_id(c)
        mtb2 = m2.metabolites.get_by_id(c)
        coef1, coef2 = -rxn1.get_coefficient(c), -rxn2.get_coefficient(c)
        mtb1.reactions
        # flux1 = sum(for r in ^, r.x * r.metabolites[mtb1])
        #flux1, flux2 = ?
        # shadow price = mtb.y; probably useful?
        # print the producing reactions too.
        buff.append('%s  (%.3f | %.3f)' % (c, coef1, coef2))
    print('(coefficient 1 | coefficient 2) => flux 1 | flux2\n%s' % '\n'.join(buff))


def visualize_all_reactions(model, flux_range, colour_range, to_file=None):
    data = [(r.id, r.x) for r in model.reactions]
    output = kegg_search_color_pathway(data, flux_range[0], flux_range[1], colour_range[0], colour_range[1], colour_range[2])
    if not to_file:
        print('\nValues:\n%s' % output)
    else:
        with open(to_file, 'w') as f:
            f.write(output)
        print('\nVisualization values for %s stored at %s' % (model, to_file))
def kegg_search_color_pathway_old(data, min_val, max_val, min_col, max_col):
    """Designed to be visualized using the Search&Color Pathway tool from KEGG, searching against the 'rn' database.
    data = [('rxn1_id', 42.3), ('rxn2_id', -332.17), ...], min_col = '#00ccee'"""
    val_delta = float(max_val - min_val)
    min_rgb = (int(min_col[1:3], 16), int(min_col[3:5], 16), int(min_col[5:7], 16))
    max_rgb = (int(max_col[1:3], 16), int(max_col[3:5], 16), int(max_col[5:7], 16))
    rgb_delta = [hc-lc for hc,lc in zip(max_rgb, min_rgb)]
    buff = []
    for r_id, val in data:
        percent = (val - min_val) / val_delta
        buff.append('%s %s' % (r_id, calc_colour(percent, min_rgb, rgb_delta)))
    return '\n'.join(buff)
def kegg_search_color_pathway(data, min_val, max_val, min_col, zero_col, max_col):
    """Designed to be visualized using the Search&Color Pathway tool from KEGG, searching against the 'rn' database.
    data = [('rxn1_id', 42.3), ('rxn2_id', -332.17), ...], min_col = '#00ccee'"""
    val_delta = float(max_val - min_val)
    min_rgb = (int(min_col[1:3], 16), int(min_col[3:5], 16), int(min_col[5:7], 16))
    max_rgb = (int(max_col[1:3], 16), int(max_col[3:5], 16), int(max_col[5:7], 16))
    zero_rgb = (int(zero_col[1:3], 16), int(zero_col[3:5], 16), int(zero_col[5:7], 16))
    pos_rgb_delta = [c-zc for c,zc in zip(max_rgb, zero_rgb)]
    neg_rgb_delta = [c-zc for c,zc in zip(min_rgb, zero_rgb)]
    buff = []
    for r_id, val in data:
        if val == 0.0:
            col = zero_col
        elif val < 0:
            percent = calc_colour_percent(val, min_val, 'sqrt')
            col = calc_colour(percent, zero_rgb, neg_rgb_delta)
        else:
            percent = calc_colour_percent(val, max_val, 'sqrt')
            col = calc_colour(percent, zero_rgb, pos_rgb_delta)
        buff.append('%s %s' % (r_id, col))
    return '\n'.join(buff)
def calc_colour_percent(val, end_val, xform=None):
    p = float(abs(val)) / abs(end_val)
    if xform == 'sqrt':
        return sqrt(p)
    else:
        return p
def calc_colour(percent, min_rgb, rgb_delta):
    col_vals = tuple(int(round(dc*percent+lc, 0)) for lc,dc in zip(min_rgb, rgb_delta))
    return '#%.2x%.2x%.2x' % (col_vals)


# # #  Options
verbose = False
min_flux = -1000
max_flux = 1000
min_colour = '#ff0000'
zero_colour = '#99ffff' # '#b3ffff' is similar, a little lighter.
max_colour = '#00ff00'
colour_range = (min_colour, zero_colour, max_colour)

test_rxns = ['R01061', 'R01512', 'R00024', 'R01523', 'R01056', 'R01641', 'R01845', 'R01829', 'R01067']
test_data = [(r, min_flux+i*(max_flux-min_flux)/(len(test_rxns)-1)) for i, r in enumerate(test_rxns)]

ov_model_file = 'model_o_vol.xlsx'
bm_model_file = 'model_b_mal.xlsx'
ov_viz_file = 'model_o_vol_flux.txt'
bm_viz_file = 'model_b_mal_flux.txt'

ov = read_excel(ov_model_file)
bm = read_excel(bm_model_file)

#basic_stats(ov)
#basic_stats(bm)
#compare_models(ov, bm)

ov.optimize()
bm.optimize()

ov.summary()
bm.summary()

#print test_data
#print kegg_search_color_pathway2(test_data, min_flux, max_flux, min_colour, zero_colour, max_colour)
visualize_all_reactions(ov, (min_flux, max_flux), colour_range, to_file=ov_viz_file)
visualize_all_reactions(bm, (min_flux, max_flux), colour_range, to_file=bm_viz_file)

compare_objective_functions(ov.reactions.get_by_id('BIOMASS'), bm.reactions.get_by_id('BIOMASS'))

# #  TO modify bounds using RNAseq, just go rpkm/max_rpkm*(upper bound); likewise for lower.
