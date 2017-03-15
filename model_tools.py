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
import os, cobra
from read_excel import read_excel, write_excel


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
    for c_id in common:
        mtb1 = m1.metabolites.get_by_id(c_id)
        mtb2 = m2.metabolites.get_by_id(c_id)
        coef1, coef2 = -coefs1[mtb1], -coefs2[mtb2]
        fluxes1 = [r.x*r.metabolites[mtb1] for r in mtb1.reactions]
        flux1 = sum(f for f in fluxes1 if f > 0)
        fluxes2 = [r.x*r.metabolites[mtb2] for r in mtb2.reactions]
        flux2 = sum(f for f in fluxes2 if f > 0)
        buff.append('%s (%.3f | %.3f) [%.1f | %.1f] => %.1f | %.1f' % (c_id, coef1, coef2, mtb1.y, mtb2.y, flux1, flux2))
    print('(coefficient 1 | coefficient 2) [shadow 1 | shadow 2] => flux 1 | flux2')
    print('%s' % '\n'.join(buff))
def analyze_shadows(model, num):
    obj_c_ids = set()
    for r in model.objective:
        for c in r.metabolites:
            obj_c_ids.add(c.id)
    mtbs_shs = []
    for c in model.metabolites:
        if c.id in obj_c_ids: continue
        mtbs_shs.append(c)
    mtbs_shs.sort(key=lambda c: c.y)

    m = model.copy()
    orig_f = m.optimize().f

    print('\nShadow prices for %s.' % model)

    neg_shs = ['\t%s %.2f'%(c.id, c.y) for c in mtbs_shs[:num]]
    print('Negative (high value):\n%s' % ('\n'.join(neg_shs)))

    #print('High value [ID ShadowPrice ObjDiff ImportFlux]:')
    #for mtb in mtbs_shs[:num]:
    #    new_f, imp_flx = single_metabolite_import_impact(m, mtb)
    #    print('\t%s %i %.1f %.1f' % (mtb.id, mtb.y, new_f-orig_f, imp_flx))

    pos_shs = ['\t%s %.2f'%(c.id, c.y) for c in mtbs_shs[:-num-1:-1]]
    print('Positive (negative value):\n%s' % ('\n'.join(pos_shs)))

def assess_metabolites_impact(model, mtbs):
    # for each never-produced mtb, should also test an export for each never-used mtb. In another function.
    m = model.copy()
    orig_f = model.solution.f
    diffs = []
    for mtb in mtbs:
        new_f, rxn_f = single_metabolite_import_impact(m, mtb)
        diff_f = new_f - orig_f
        if diff_f > 1.0:
            diffs.append((diff_f, rxn_f, mtb.id, mtb.name))
    diffs.sort(key=lambda d: -d[0])
    for d in diffs:
        print('Delta: %i Flux: %i (%s) %s' % d)
    print('%i reactions improved the objective function.' % len(diffs))

def single_metabolite_import_impact(model, mtb):
    """Adds a reversible reaction to import and runs pFBA. Returns the new
    objective flux and the flux through the import reaction. Removes the import
    reaction before return.
    """
    rxn_id = 'TEST_IMPORT_REACTION'
    if rxn_id in model.reactions:
        i = 1
        while '%s_%i' % (rxn_id, i) in model.reactions:
            i += 1
        rxn_id = '%s_%i' % (rxn_id, i)
    rxn = cobra.Reaction(rxn_id)
    rxn.lower_bound = -1000
    rxn.upper_bound = 1000
    rxn.objective_coeficient = 0
    if isinstance(mtb, str):
        mtb = model.metabolites.get_by_id(mtb)
    rxn.add_metabolites({mtb: 1.0})
    model.add_reaction(rxn)
    obj_f = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model).f
    rxn_f = rxn.x
    model.remove_reactions([rxn])
    return obj_f, rxn_f


def visualize_fba_reactions(model, flux_range, colour_range, to_file_str=None):
    data = [(r.id, r.x) for r in model.reactions]
    data.sort(key=lambda r: r[0])
    output = kegg_search_color_pathway_fba(data, flux_range[0], flux_range[1], colour_range[0], colour_range[1], colour_range[2])
    if not to_file_str:
        print('\nValues:\n%s' % output)
    else:
        to_file = to_file_str % 'fba'
        with open(to_file, 'w') as f:
            f.write(output)
        print('\nFBA flux values for %s stored at %s' % (model, to_file))
def kegg_search_color_pathway_fba(data, min_val, max_val, min_col, zero_col, max_col):
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
    elif xform == 'inv': # Used to effectively reverse the gradient direction.
        return 1.0 - p
    else:
        return p
def calc_colour(percent, min_rgb, rgb_delta):
    return '#%.2x%.2x%.2x' % (calc_colour_rgb(percent, min_rgb, rgb_delta))
def calc_colour_rgb(percent, min_rgb, rgb_delta):
    return tuple(int(round(dc*percent+lc, 0)) for lc,dc in zip(min_rgb, rgb_delta))

def visualize_fva_reactions(fva_rxns, flux_range, colour_range, to_file_str):
    data = [(r, vals['minimum'], vals['maximum']) for r,vals in fva_rxns.items()]
    data.sort(key=lambda r: r[0])
    output = kegg_search_color_pathway_fva(data, flux_range[0], flux_range[1], colour_range[0], colour_range[1], colour_range[2])
    to_file = to_file_str % 'fva'
    with open(to_file, 'w') as f:
        f.write(output)
    print('\nFVA values written to %s' % (to_file))
def kegg_search_color_pathway_fva(data, min_val, max_val, min_col, zero_col, max_col):
    """This uses the same colour gradient as the FBA variant, but it also maps
    flux variance. The average of the min and max is set as the 'value' for the
    colour gradient, which is then modified by the variance, towards black. So
    the darker the colour, the higher the variance through the reaction."""
    val_delta = float(max_val - min_val)
    min_rgb = (int(min_col[1:3], 16), int(min_col[3:5], 16), int(min_col[5:7], 16))
    max_rgb = (int(max_col[1:3], 16), int(max_col[3:5], 16), int(max_col[5:7], 16))
    zero_rgb = (int(zero_col[1:3], 16), int(zero_col[3:5], 16), int(zero_col[5:7], 16))
    black_rgb = (0, 0, 0)
    pos_rgb_delta = [c-zc for c,zc in zip(max_rgb, zero_rgb)]
    neg_rgb_delta = [c-zc for c,zc in zip(min_rgb, zero_rgb)]
    buff = []
    for r_id, min_var, max_var in data:
        val = (min_var + max_var) / 2.0
        if val == 0.0:
            col_rgb = zero_rgb
        elif val < 0:
            percent = calc_colour_percent(val, min_val, xform='sqrt')
            col_rgb = calc_colour_rgb(percent, zero_rgb, neg_rgb_delta)
        else:
            percent = calc_colour_percent(val, max_val, xform='sqrt')
            col_rgb = calc_colour_rgb(percent, zero_rgb, pos_rgb_delta)
        var_prcnt = calc_colour_percent(max_var-min_var, val_delta, xform='inv')
        var_col = calc_colour(var_prcnt, black_rgb, col_rgb)
        buff.append('%s %s' % (r_id, var_col))
    return '\n'.join(buff)

def pathway_analysis(models, fvas, pathways):
    range_str = '%i to %i'
    desc_str = '\nFluxes: %s' % (' | '.join('%s %.1f'%(str(m), m.solution.f) for m in models))
    print(desc_str)
    print('-' * len(desc_str.strip()))
    for path_desc, rxn_ids in pathways: # for each pathway:
        max_name_len, val_lens = 0, [0]*len(fvas)
        buff = []
        for rxn_data in rxn_ids: # for each reaction in the pathway:
            name, r_ids = rxn_data[:2]
            rxn_coef = 1.0
            if len(rxn_data) == 3:
                rxn_coef *= rxn_data[2]
            if len(name) > max_name_len:
                max_name_len = len(name)
            vals = [name+':']
            for i, fva in enumerate(fvas): # for each model:
                for r_id in r_ids: # one of these IDs should be in that model.
                    data = fva.get(r_id, None)
                    if data != None:
                        min_max = (int(round(data['minimum'])), int(round(data['maximum'])))
                        if rxn_coef < 0:
                            min_max = (int(round(min_max[1]*rxn_coef)), int(round(min_max[0]*rxn_coef)))
                        val_str = range_str % min_max
                        break
                else: # none of the reaction IDs were found in the model.
                    val_str = 'N/A'
                if len(val_str) > val_lens[i]:
                    val_lens[i] = len(val_str)
                vals.append(val_str)
            buff.append(tuple(vals))
        value_str = '%%-%is %s' % (max_name_len+1, ' | '.join('%%%is' % l for l in val_lens))
        pathway_values = '\n'.join(value_str % v for v in buff)
        print('\t%s:\n%s\n' % (path_desc, pathway_values))

def compare_reaction_mtbs(model1, model2, missing_mtb, outfile):
    rxn_mtbs = set()
    exchanges, transports = set(), set()
    for rxn in model1.reactions:
        rs = tuple(sorted(m.id.strip()[1:] for m in rxn.reactants))
        ps = tuple(sorted(m.id.strip()[1:] for m in rxn.products))
        if any(m[0].isalpha() for m in rs+ps): # Using a custom metabolite.
            continue
        if rs == ps:
            transports.add(rs)
        elif not rs:
            exchanges.add(ps)
        elif not ps:
            exchanges.add(rs)
        else:
            rxn_mtbs.add((rs, ps))
    print('\nParsed %i reactions, %i transports, and %i exchanges.' % (len(rxn_mtbs), len(transports), len(exchanges)))
    missed_exs = set()
    bad_ids = set()
    missing_mtb_ids = set()
    for rxn in model2.reactions:
        rs = tuple(sorted(m.id.strip()[1:] for m in rxn.reactants))
        ps = tuple(sorted(m.id.strip()[1:] for m in rxn.products))
        if not rs or not ps:
            ms = rs if rs else ps
            if ms not in exchanges:
                missed_exs.add(ms)
        else:
            if (rs, ps) in rxn_mtbs or (ps, rs) in rxn_mtbs:
                pass
            else:
                coef = check_rxn_missing_mtb(missing_mtb, rs, ps, rxn_mtbs)
                if coef:
                    missing_mtb_ids.add(rxn.id)
                    rxn.add_metabolites({missing_mtb: coef})
                else:
                    bad_ids.add(rxn.id)
    print('%s is missing %i exchanges present in %s, with %i inconsistencies.' % (model2, len(missed_exs), model1, len(bad_ids)))
    write_excel(model2, outfile)
    print('%i reactions were missing %s; modified and saved as %s' % (len(missing_mtb_ids), missing_mtb, outfile))
def check_rxn_missing_mtb(mtb, rs, ps, rxn_mtbs):
    rs_h = tuple(sorted(rs + (mtb[1:],)))
    ps_h = tuple(sorted(ps + (mtb[1:],)))
    if (rs_h, ps) in rxn_mtbs or (ps, rs_h) in rxn_mtbs:
        return -1.0 # mtb should be a reactant
    elif (rs, ps_h) in rxn_mtbs or (ps_h, rs) in rxn_mtbs:
        return 1.0 # mtb should be a product
    else:
        return 0 # neither reaction modification found



# # #  Parameters
min_flux = -1000
max_flux = 1000
min_colour = '#ff0000'
zero_colour = '#ffff00'
max_colour = '#00ff00'
colour_range = (min_colour, zero_colour, max_colour)

# # #  Run-time options
files_dir = '/mnt/hgfs/win_projects/brugia_project'
model_files = ['model_o_vol_2.5.xlsx', 'model_b_mal_2.5.xlsx', 'model_o_vol_3.xlsx']
mtb_cmp_str = '%s_wip.xlsx'
verbose = True
topology_analysis = False
fba_analysis = False
fva_analysis = True
save_visualizations = False
do_reaction_mtb_comparison = None # 'C00080', or None

# # #  Pathway analysis
pathways_for_analysis = [
    ('Fundamental imports', [
        ('Glucose', ('CARBON_SOURCE','EX00031')),
        ('Oxygen', ('DIFFUSION_2','EX00007'), -1),
        ('CO2', ('DIFFUSION_3','EX00011'), -1),
        ('Water', ('DIFFUSION_1','EX00001'), -1),
        ('Phosphate', ('DIFFUSION_6','EX00009'), -1),
        ('Bicarb', ('DIFFUSION_8','EX00288'))
    ]),
    ('TCA cycle', [
        ('oaa -> cit', ('R00351',), -1),
        ('cit -> icit 1a', ('R01325',)),
        ('cit -> icit 1b', ('R01900',), -1),
        ('icit -> akg NADPH', ('R00267',)),
        ('icit -> akg NADH', ('R00709',)),
        ('akg -> succoa', ('R02570',), -1),
        ('succoa -> succ', ('R00405',), -1),
        ('succ -> fum', ('R02164',)),
        ('fum -> mal', ('R01082',), -1),
        ('mal -> oaa', ('R00342',))
    ]),
    ('Oxidative phosphorylation', [
        ('Complex I', ('R02163',), -1),
        ('Complex II', ('R02164',)),
        ('Complex III', ('R02161',)),
        ('Complex IV', ('R00081',)),
        ('ATP synthase', ('R00086',), -1),
        ('C II reverse', ('R01867',)),
        ('d-oro -> UQH2', ('R01868',))
    ])
]
# ('Non-growth', ('NGAM',))

# # #  Program objects
test_rxns = ['R01061', 'R01512', 'R00024', 'R01523', 'R01056', 'R01641', 'R01845', 'R01829', 'R01067']
test_data = [(r, min_flux+i*(max_flux-min_flux)/(len(test_rxns)-1)) for i, r in enumerate(test_rxns)]

# # #  Run steps
cel_m = read_excel(os.path.join(files_dir, 'iCEL1273.xlsx'), verbose=False)
cel_m.reactions.BIO0103.objective_coefficient = 1.0 # # #  TESTING ONLY

model_files = [os.path.join(files_dir, m_file) for m_file in model_files]
models = [read_excel(m_file, verbose=verbose) for m_file in model_files]

if do_reaction_mtb_comparison:
    new_model_files, new_models = [], []
    cel_m = read_excel(os.path.join(files_dir, 'iCEL1273.xlsx'), verbose=False)
    cel_m.reactions.BIO0103.objective_coefficient = 1.0
    for m_file, model in zip(model_files, models):
        mtb_outfile = mtb_cmp_str % m_file.rpartition('.')[0]
        m = model.copy()
        compare_reaction_mtbs(cel_m, m, do_reaction_mtb_comparison, mtb_outfile)
        new_model_files.append(mtb_outfile)
        new_models.append(read_excel(mtb_outfile))
    model_files += new_model_files
    models += new_models

viz_strs = [os.path.split(m_file.rpartition('.')[0]+'_%s.txt')[1] for m_file in model_files]

for m in models:
    cobra.flux_analysis.parsimonious.optimize_minimal_flux(m)
    # rather than m.optimize(); returns optimal FBA with minimum total flux through network.

if topology_analysis:
    for m in models:
        basic_stats(m)
    # # compare_models(m1, m2)

if verbose:
    for m in models:
        print('\n%s summary:' % m)
        m.summary()

if fba_analysis:
    for m, viz_str in zip(models, viz_strs):
        analyze_shadows(m, 5)
        if save_visualizations:
            visualize_fba_reactions(m, (min_flux, max_flux), colour_range, to_file_str=viz_str)
    # # compare_objective_functions(m1.reactions.get_by_id('BIOMASS'), m2.reactions.get_by_id('BIOMASS'))

fvas = []
if fva_analysis:
    for m, viz_str in zip(models, viz_strs):
        m_fva = cobra.flux_analysis.flux_variability_analysis(m)
        fvas.append(m_fva)
        if save_visualizations:
            visualize_fva_reactions(m_fva, (min_flux, max_flux), colour_range, viz_str)
    pathway_analysis(models, fvas, pathways_for_analysis)

#assess_metabolites_impact(models[0], models[0].metabolites)


# Add pyruvate (C00022) to pathway analysis.
#  - Can go into TCA, anaerobic respiration, produce glucose or other carbs, produce fatty acids, produce alanine, etc.
#  - Monitor these.

# loopless_model = cobra.flux_analysis.loopless.construct_loopless_model(model)
#  - optimize() hadn't completed after 40 hours.
#  - Running fva on this model is slow; hadn't completed after 40 hours.
#  - Better idea is do fva on regular model, remove all reactions with no flux, find loopless model, run fva.
# cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
#  - Does FBA, then finds solution that also minimizes total flux in the system. Returns nothing.
# essential_transports = cobra.flux_analysis.essentiality.assess_medium_component_essentiality(model, the_components=[], solver='cglpk')

# Dead-end reactions are quite common in metaNets.
#  - Could iterate over metabolites, adding sink or transport reactions. See which have an impact on the f.
#  - Manually, search kegg for rxns with unavailable reactants. See where they are, if some are the same pathways, etc. Might indicate missing reactions.
#  - Or, for unavailable mtbs, see if any are connected by reactions. then group them together, see if I end up at an unusable mtb.

# A good way to compare models might be build a new network from a solution.
#  - Start at obj fxn, backtrack to transport reactions using those with flux.
#  - Then when comparing two of these nets, identify the portions that take alternate paths to end up at some common intermediate.
# Related: after optimizing, implement method to trace back from the obj fxn to each mtb. Keep path of rxns (use to watch out for cycles).
#  - Should

# Might be a good idea to compare my models to the published human ones (at least for the general metabolism).

# To modify constraints using RNAseq, just go rpkm/max_rpkm*(upper bound); likewise for lower.
#  - However, this doesn't account for different reactions having different speeds (kd?).
#  - Maybe use a less-strict transformation, like sqrt(rpkm/max).
#  - Could also apply the transcriptomic constraints, then go through reactions and try reverting each. See which has the biggest impact on obj fxn; maybe that one should be relaxed.
