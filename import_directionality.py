import pandas
from cobra import Model, Metabolite, Reaction
from read_excel import read_excel, write_excel

elegans_file = 'iCEL1273.xlsx'
model_in_file = 'model_o_vol.xlsx'
model_out_file = 'model_o_vol_2-wip.xlsx'

reversible_arrow = '<==>'
irreversible_arrow = ('-->', '<--')
l_bound, u_bound = -1000, 1000

old_m = read_excel(model_in_file)

pio = pandas.io.excel.ExcelFile(elegans_file)
rxn_frame = pandas.read_excel(elegans_file, 'Reactions')

def use_directionalities(old_m, disagreements):
    unresolved = []
    for r_id, cel_dir, m_dir, cel_products in disagreements:
        updated = False
        if cel_dir == 'reversible':
            rxn = old_m.reactions.get_by_id(r_id)
            rxn.bounds = (l_bound, u_bound)
            updated = True
        elif cel_dir == 'irreversible':
            if ',' not in r_id: # combined reactions are handled manually
                cel_prods = set()
                for m in cel_products:
                    met = m.replace('M','C').replace('E','C')
                    if ' ' in met:
                        met = met.partition(' ')[2]
                    cel_prods.add(met)
                rxn = old_m.reactions.get_by_id(r_id)
                prds = set(m.id for m in rxn.products)
                rcts = set(m.id for m in rxn.reactants)
                if prds == cel_prods or prds|set(['C00080']) == cel_prods:
                    rxn.bounds = (0, u_bound)
                    updated = True
                elif rcts == cel_prods or rcts|set(['C00080']) == cel_prods:
                    rxn.bounds = (l_bound, 0)
                    updated = True
        if not updated:
            unresolved.append((r_id, cel_dir, m_dir, cel_products))
    return unresolved


processed = set()
disagreements, rvsb, irrvsb, broken = [], 0, 0, 0
for i in range(len(rxn_frame)):
    row = rxn_frame.ix[i]
    rxn_str = row['Machine readable']
    r_ids = row['KEGG']
    if not isinstance(r_ids, unicode) or not isinstance(r_ids, unicode):
        continue
    if r_ids in processed:
        continue
    else:
        processed.add(r_ids)
    r_ids = r_ids.split(';')
    if len(rxn_str) == 0 or len(r_ids) == 0:
        continue
    products = set()
    if reversible_arrow in rxn_str:
        cel_direction = 'reversible'
    elif irreversible_arrow[0] in rxn_str:
        cel_direction = 'irreversible'
        products = set(str(p.strip()) for p in rxn_str.partition(irreversible_arrow[0])[2].split('+'))
    elif irreversible_arrow[1] in rxn_str:
        cel_direction = 'irreversible'
        products = set(str(p.strip()) for p in rxn_str.partition(irreversible_arrow[1])[0].split('+'))
    else:
        print('\nError: could not parse %s' % rxn_str)
    old_directions = []
    for r_id in r_ids:
        if r_id in old_m.reactions:
            rxn = old_m.reactions.get_by_id(r_id)
            if rxn.lower_bound < 0:
                if rxn.upper_bound <= 0:
                    old_directions.append('irreversible')
                else:
                    old_directions.append('reversible')
            else:
                if rxn.upper_bound > 0:
                    old_directions.append('irreversible')
                else:
                    old_directions.append('blocked')
        else:
            old_directions.append('missing')
    agreed = False
    if len(set(old_directions)) == 1 and old_directions[0] == 'missing':
        continue # No matching reactions found in old_m
    elif 'blocked' in old_directions or 'missing' in old_directions:
        broken += 1 # Part of the reaction set exists in old_m, but part is missing or blocked.
    elif cel_direction == 'reversible':
        if 'irreversible' not in old_directions:
            agreed = True
        else:
            rvsb += 1
    else: # Should be irreversible
        if 'irreversible' in old_directions:
            agreed = True
        else:
            irrvsb += 1
    if not agreed:
        #print('%s should be %s, but is: %s' % (', '.join(r_ids), cel_direction, ', '.join(old_directions)))
        disagreements.append((', '.join(r_ids), cel_direction, ', '.join(old_directions), products))
disagreements.sort(key=lambda d: d[0])
disagreements.sort(key=lambda d: d[1])
print('\n%i should have been reversible, %i irreversible, and %i were broken' % (rvsb, irrvsb, broken))

unresolved = use_directionalities(old_m, disagreements)
print('\n'.join('%s should be %s, but is: %s' % (d[:3]) for d in unresolved))
print len(unresolved)
write_excel(old_m, model_out_file)
