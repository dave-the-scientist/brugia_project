import os, cobra
from read_excel import read_excel, write_excel
from model_tools import id_bottleneck_metabolites


def compare_constraints(good_m, test_m, test_fva, test_biomass_impact):
    good_f, orig_f = good_m.optimize().f, test_m.optimize().f
    print('Originally %s yielded %.2f' % (test_m, orig_f))
    for rxn in good_m.reactions:
        r_id = rxn.id
        t_rxn = test_m.reactions.get_by_id(r_id)
        good_bounds, test_bounds = rxn.bounds, t_rxn.bounds
        if good_bounds != test_bounds:
            t_rxn.bounds = good_bounds
            new_f = test_m.optimize().f
            if new_f - orig_f > 1.0:
                print r_id, new_f, good_bounds, test_bounds

                    print_biomass_impact(r_id, test_m, test_fva, good_f)

            t_rxn.bounds = test_bounds
    if test_biomass_impact:
        diffs = id_bottleneck_metabolites(test_m, orig_f, 'BIOMASS', threshold=1.0)
        # producing_r_id = mtb.id.upper() # used to go a level deeper in biomass rxn.


files_dir = '/mnt/hgfs/win_projects/brugia_project'
good_model = 'model_b_mal_4.5-wip.xlsx'
test_model = 'model_b_mal_5_M30.xlsx'

test_biomass_impact = True

good_m = read_excel(os.path.join(files_dir, good_model))
test_m = read_excel(os.path.join(files_dir, test_model))
test_fva = cobra.flux_analysis.flux_variability_analysis(test_m)

compare_constraints(good_m, test_m, test_fva, test_biomass_impact)
