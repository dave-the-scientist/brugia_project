import os, cobra
from read_excel import read_excel, write_excel

def modify_model(model, cel_m, do_deletions=True):
    rxns_to_add, custom_rxns = {}, {}
    for rxn in cel_m.reactions:
        r_id = rxn.id
        if r_id.startswith('RMC'):
            if parse_unadded_keggs(rxn, '_M', model): # True if at least of the kegg ids has not been implemented.
                custom_rxns.setdefault(r_id, []).append(rxn)
        elif r_id.startswith('RM'):
            keggs = parse_unadded_keggs(rxn, '_M', model)
            for kegg in keggs:
                rxns_to_add.setdefault(kegg, []).append(rxn)
        else:
            continue
    cyto_mod, cyto_rm = 0, 0
    missing_mito, rxns_to_delete = {}, []
    for kegg, rxns in rxns_to_add.items():
        if len(rxns) != 1:
            print('Multiple M reactions for id %s' % (kegg))
            continue
        if kegg in model.reactions: # a cyto version already exists.
            cyto_rxn = model.reactions.get_by_id(kegg)
            if copy_cyto_to_mito(kegg, model, cyto_rxn, rxns[0]):
                cyto_mod += 1
                if 'RC' + kegg[1:] not in cel_m.reactions: # cyto version shouldn't exist.
                    if do_deletions:
                        cyto_rxn.delete()
                        cyto_rm += 1
                    else:
                        rxns_to_delete.append(cyto_rxn)
            else:
                pass # cyto version already correct (generally exchange transports).
        else: # cyto version doesn't exist
            missing_mito[kegg] = rxns[0]
            # uncomment to add
            #if not add_mito_reaction(kegg, '_M', model, rxns[0]):
            #    custom_rxns.setdefault(kegg, []).append(rxns[0])
    cust_ids = sorted(cid for cid in custom_rxns)
    print('\n%i custom reactions:\n%s' % (len(cust_ids), ','.join(cust_ids)))
    print('\n%s: %i modified and %i added a dual form; %i missing.\n' % (model, cyto_rm, cyto_mod-cyto_rm, len(missing_mito)))
    return rxns_to_delete

def parse_unadded_keggs(rxn, suffix, model):
    keggs = []
    kegg_str = getattr(rxn, 'kegg_reaction', '')
    if not kegg_str or kegg_str == 'NA':
        return keggs
    for kegg in kegg_str.split(';'):
        if kegg + suffix not in model.reactions:
            keggs.append(kegg)
    return keggs
def copy_cyto_to_mito(kegg, model, cyto_rxn, cel_rxn):
    new_id = kegg + '_M'
    if cyto_rxn.reaction == cel_rxn.reaction:
        return False
    rxn = cobra.Reaction(new_id)
    rxn.name = cyto_rxn.name
    rxn.bounds = cyto_rxn.bounds
    rxn.subsystem = cyto_rxn.subsystem
    rxn.enzyme_commission = cyto_rxn.enzyme_commission
    rxn.gene_names = cyto_rxn.gene_names
    rxn.confidence_notes = cyto_rxn.confidence_notes
    rxn.reaction_notes = cyto_rxn.reaction_notes
    mtb_dict = {}
    for mtb, coef in cyto_rxn.metabolites.items():
        if mtb.id[0] == 'C':
            new_mtb_id = 'M' + mtb.id[1:]
            mtb_dict[cobra.Metabolite(new_mtb_id, name=mtb.name)] = coef
        else:
            print cyto_rxn.reaction
    rxn.add_metabolites(mtb_dict)
    model.add_reaction(rxn)
    return True

def add_mito_reaction(kegg, suffix, model, cel_rxn):
    for mtb in cel_rxn.metabolites:
        if mtb.id.strip()[1] == 'C': # using custom metabolites.
            return False
    new_id = kegg + suffix
    rxn = cobra.Reaction(new_id)
    rxn.name = cel_rxn.name
    rxn.bounds = cel_rxn.bounds
    rxn.subsystem = cel_rxn.subsystem
    rxn.enzyme_commission = cel_rxn.enzyme_commission
    model.add_reaction(rxn)
    rxn.build_reaction_from_string(cel_rxn.reaction)
    return True


files_dir = '/mnt/hgfs/win_projects/brugia_project'
model_files = ['model_o_vol_3.xlsx', 'model_b_mal_3.xlsx']
out_files = ['model_o_vol_3.5.xlsx', 'model_b_mal_3.5.xlsx']
do_deletions = True

cel_m = read_excel(os.path.join(files_dir, 'iCEL1273.xlsx'), verbose=False)
#cel_m.reactions.BIO0101.objective_coefficient = 1.0 # # #  TESTING ONLY
#cobra.flux_analysis.parsimonious.optimize_minimal_flux(cel_m)
models = [read_excel(os.path.join(files_dir, m_file), verbose=False) for m_file in model_files]

for m, out_file in zip(models, out_files):
    if do_deletions:
        modify_model(m, cel_m, do_deletions)
        write_excel(m, os.path.join(files_dir, out_file))
    else:
        print('Obj before modifications: %.1f' % (m.optimize().f))
        rxns_to_delete = modify_model(m, cel_m, do_deletions)
        diffs = []
        cur_f = m.optimize().f
        for rxn in rxns_to_delete:
            rxn.delete()
            diff_f = m.optimize().f - cur_f
            diffs.append((diff_f, rxn))
            cur_f += diff_f
            if cur_f < 1.0:
                break
        diffs.sort()
        for d in diffs:
            print('%.1f by %s' % (d[0], d[1].id))
