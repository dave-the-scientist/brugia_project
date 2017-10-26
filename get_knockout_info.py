"""
Notes:
- Brugia protein sequences: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA10729
- wBm protein sequences: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=292805
"""
import os, cPickle, pandas
from molecbio import sequ
from cobra.flux_analysis import single_reaction_deletion, double_reaction_deletion
from model_tools import id_bottleneck_metabolites
from read_excel import read_excel
import xml.etree.ElementTree as ET


def get_rxns_to_delete(model):
    rxn_to_genes = {}
    for rxn in model.reactions:
        if not rxn.gene_names or not rxn.id.startswith(('R', 'ACYLCOA')):
            continue
        rxn_to_genes[rxn.id] = [g.strip() for g in rxn.gene_names.split(';')]
    return rxn_to_genes

def do_deletions(rxn_data, model, rxn_to_genes, do_double_ko=False, obj_fraction=0.0):
    fraction_epsilon = 0.0001
    orig_f = float(model.optimize().f)
    s_rates, s_stats = single_reaction_deletion(model, list(rxn_to_genes.keys()))
    print('Original objective %.1f; %i reactions knocked out.' % (orig_f, len(s_stats)))
    print('Calculating model deficiencies for each knockout...')
    for r_id, new_f in s_rates.items():
        if abs(new_f) < fraction_epsilon:
            new_f = 0.0
        stat = s_stats[r_id]
        if new_f/orig_f <= obj_fraction+fraction_epsilon:
            if stat == 'optimal':
                deficiencies = find_model_deficiencies(model, orig_f, new_f, r_id)
            else:
                deficiencies = 'infeasible'
            rxn_data[r_id] = {'objective':round(new_f/orig_f*100, 1), 'deficiencies':deficiencies, 'genes':rxn_to_genes[r_id]}
    if do_double_ko:
        double_rxn_ids = [r for r in list(rxn_to_genes.keys()) if r not in rxn_data]
        print('Performing double knockouts on %i candidates...' % len(double_rxn_ids))
        double_ko_data = double_reaction_deletion(model, double_rxn_ids[:5], number_of_processes=3)
        d_r1, d_r2, d_rates = double_ko_data['y'], double_ko_data['x'], double_ko_data['data']

def find_model_deficiencies(model, orig_f, new_f, r_id):
    deficiencies = []
    ob = model.reactions.get_by_id(r_id).bounds
    model.reactions.get_by_id(r_id).bounds = (0,0)
    diffs = id_bottleneck_metabolites(model, new_f, 'BIOMASS', threshold=1.0)
    for recovered_f, mtb_id in diffs:
        def_str = '%s (%.1f)' % (mtb_id, recovered_f/orig_f*100)
        sub_defs = []
        for sub_f, sub_mtb_id in id_bottleneck_metabolites(model, new_f, mtb_id.upper(), threshold=1.0):
            sub_defs.append('%s(%.1f)' % (sub_mtb_id, sub_f/orig_f*100))
        if sub_defs:
            def_str += ' [%s]' % ', '.join(sub_defs)
        deficiencies.append(def_str)
    model.reactions.get_by_id(r_id).bounds = ob
    if not deficiencies:
        return 'unrecoverable'
    else:
        return ', '.join(deficiencies)

def process_gene_data(rxn_data):
    gene_data = {}
    for r_id, data in rxn_data.items():
        for gene in set(data['genes']):
            g_entry = generate_gene_entry(data, r_id, gene)
            gene_data.setdefault(gene, []).append(g_entry)
    for gene, entries in gene_data.items():
        rs_per_g = len(entries)
        if rs_per_g > 1:
            for e in entries:
                e['num_reactions'] = rs_per_g
    return gene_data

def generate_gene_entry(r_data, r_id, gene):
    g_data = {}
    if len(set(r_data['genes'])) == 1:
        g_data['other_genes'] = ''
    else:
        g_data['other_genes'] = ','.join(sorted(list(set(r_data['genes']) - set([gene]))))
    g_data['reaction'] = r_id
    g_data['objective'] = r_data['objective']
    g_data['deficiencies'] = r_data['deficiencies']
    g_data['num_reactions'] = 1
    return g_data

# # #  Save/load functions
def save_data_object(data_obj, file_path):
    with open(file_path, 'wb') as f:
        cPickle.dump(data_obj, f, protocol=0)
    print('Saved data to %s' % file_path)
def load_data_object(file_path):
    with open(file_path, 'rb') as f:
        data_obj = cPickle.load(f)
    print('Loaded data from %s' % file_path)
    return data_obj

def save_data_to_excel(gene_data, gene_data_out_file, expression_headings):
    min_column_width = 10
    sheet_name = 'Single knockouts'
    gene_header = 'Gene ID'
    headers_atts = [('# Reactions','num_reactions'), ('Reaction','reaction'), ('Associated genes','other_genes'), ('Objective %','objective'), ('Biomass deficiencies','deficiencies')]
    ortho_headers = ['Human homologs', 'Homolog quality (% identity|% coverage)']
    data = {h[0]:[] for h in headers_atts+expression_headings}
    for h in [gene_header] + ortho_headers:
        data[h] = []
    gene_order = sorted(list(gene_data.keys()))
    gene_order.sort(key=lambda g:gene_data[g][0]['deficiencies'])
    for gene in gene_order:
        for g_data in gene_data[gene]:
            data[gene_header].append(gene)
            for h, att in headers_atts:
                data[h].append(g_data.get(att, 'NOT FOUND'))
            data[ortho_headers[0]].append(g_data['num_human_prots'])
            data[ortho_headers[1]].append( '%.1f|%.1f' % (g_data['human_prot_identity'],g_data['human_prot_coverage']) )
            for h, ls in expression_headings:
                exp_levels = [g_data['expression_levels'].get(l) for l in ls]
                data[h].append('|'.join(exp_levels))
    col_headers = [gene_header] + [h[0] for h in headers_atts] + [i for i in ortho_headers] + [j[0] for j in expression_headings]
    writer = pandas.ExcelWriter(gene_data_out_file, engine='xlsxwriter')
    df = pandas.DataFrame(data)[col_headers] # The [] specifies the order of the columns.
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    worksheet = writer.sheets[sheet_name]
    for i, h in enumerate(col_headers):
        col_width = max(len(h)+1, min_column_width)
        worksheet.set_column(i, i, col_width)
    writer.save()
    print('Data saved to %s' % gene_data_out_file)

# # #  Getting protein names and sequences
def save_prot_names_list(gene_data):
    prot_list_file = 'utility/b_mal_4.5-wip_single_ko_prot_names.txt'
    prot_list = sorted(gene_data.keys())
    with open(prot_list_file, 'w') as f:
        f.write('\n'.join(prot_list))
    print('Saved protein list to %s' % prot_list_file)

def get_prot_name_translations(gene_data, gen_pept_file):
    print('Parsing %s...' % gen_pept_file)
    prot_to_std, found_names = {}, set()
    with open(gen_pept_file, 'r') as f:
        prot_name, std_name = None, None
        for line in f:
            if prot_name == None and line.startswith('VERSION'):
                prot_name = line.strip().split()[1]
            elif prot_name and "/standard_name=" in line:
                std_name = line.partition('=')[2].strip()[1:-1]
                if std_name in gene_data:
                    prot_to_std[prot_name] = std_name
                    found_names.add(std_name)
                prot_name, std_name = None, None
    for gene in gene_data:
        if gene not in found_names:
            prot_to_std['%s.1' % gene] = gene
    return prot_to_std

def save_prot_sequences(gene_data, prot_to_std, prot_sequences_file):
    prots_fasta_file = 'utility/b_malayi_and_wBm_prots.fasta'
    all_seqs = sequ.loadfasta(prots_fasta_file)
    prots, found_genes = [], set()
    for seq in all_seqs:
        gene = prot_to_std.get(seq.name)
        if not gene: continue
        if gene in found_genes:
            print('Error: multiple sequences were found matching "%s".' % seq.name)
            exit()
        prots.append(sequ.Sequence(name=gene, sequence=seq.seq))
        found_genes.add(gene)
    if len(prots) != len(gene_data):
        print('Warning: only found sequences for %i of %i genes. Missing genes:' % (len(prots), len(gene_data)))
        for g in set(gene_data) - found_genes:
            print(g)
        exit()
    sequ.savefasta(prots, prot_sequences_file)
    print('Saved %i sequences to %s' % (len(prots), prot_sequences_file))
    return prots

# # #  Parsing BLAST output
def parse_blast_xml(gene_data, blast_xml_file):
    min_e_val = 1E-30
    tree = ET.parse(blast_xml_file)
    root = tree.getroot()
    iterations = root.find('BlastOutput_iterations')
    for q_hit in iterations:
        gene = q_hit.find('Iteration_query-def').text
        prot_len = float(q_hit.find('Iteration_query-len').text)
        s_hits = q_hit.find('Iteration_hits')
        num_hits, top_hit_id, top_e_val, top_identity, top_coverage = get_good_hits(s_hits, min_e_val)
        top_coverage = round(top_coverage/prot_len*100.0, 1)
        for g_data in gene_data[gene]:
            g_data['num_human_prots'] = num_hits
            g_data['human_prot_id'] = top_hit_id
            g_data['human_prot_identity'] = top_identity
            g_data['human_prot_coverage'] = top_coverage


def get_good_hits(s_hits, min_e_val):
    """Counts based on the 'Hit_def' field in the subject hits, which is the name. Attempts to remove isoforms and predicted proteins from the count.
    """
    best_hit_id, best_e_val, best_ident, best_coverage = None, min_e_val + 1, 0, 0
    hit_names = set()
    for s_hit in s_hits:
        hit_e_val, hit_ident, hit_coverage = min_e_val + 1, 0, 0
        for hsp in s_hit.find('Hit_hsps'):
            e_val = float(hsp.find('Hsp_evalue').text)
            if e_val < hit_e_val:
                hit_e_val = e_val
                hit_ident = round(float(hsp.find('Hsp_identity').text)/float(hsp.find('Hsp_align-len').text)*100, 1)
                hit_coverage = int(hsp.find('Hsp_query-to').text) - int(hsp.find('Hsp_query-from').text)
        if hit_e_val < min_e_val:
            for nm in s_hit.find('Hit_def').text.lower().split(' >gi'):
                if nm.strip().endswith(' [homo sapiens]'):
                    name = nm.strip()[:-15]
                    break
            else:
                print('Warning: no human protein name found in hit "%s", skipping it.' % s_hit.find('Hit_def').text)
                continue # Something weird going on, probably don't want to count it.
            if name.startswith('predicted: '):
                name = name[11:]
            name = name.partition(' isoform ')[0]
            hit_names.add(name)
            if hit_e_val < best_e_val:
                best_hit_id = s_hit.find('Hit_accession').text.strip()
                best_ident = hit_ident
                best_e_val, best_coverage = hit_e_val, hit_coverage
    if not hit_names:
        return 0, None, None, 0, 0
    return len(hit_names), best_hit_id, best_e_val, best_ident, best_coverage

# # #  Getting expression data
def get_expression_data(gene_data, expression_file, sheetnames, conditions):
    for sheetname in sheetnames:
        parse_expression_sheet(gene_data, expression_file, sheetname, conditions)
    null_exp = {c:' -- ' for c in conditions}
    for gene, entries in gene_data.items():
        for e in entries:
            if 'expression_levels' not in e:
                e['expression_levels'] = null_exp

def parse_expression_sheet(gene_data, filename, sheetname, conditions):
    seq_name_key = 'Sequence Name'
    replicate_inds = ['a', 'b', 'c']
    frame = pandas.read_excel(filename, sheetname)
    if len(frame.columns) != len(set(frame.columns)):
        print('Error: at least one column header was not unique in sheet %s.' % sheetname)
        exit()
    cond_keys = [[cond+ind for ind in replicate_inds if cond+ind in frame.columns] for cond in conditions]
    for i in frame.index:
        row = frame.ix[i]
        seq_name = row[seq_name_key]
        if seq_name not in gene_data:
            continue
        avgs = [sum(row[k] for k in ck)/float(len(ck)) for ck in cond_keys]
        exp = {c:'%.1f'%(a/max(avgs)*100.0) for a,c in zip(avgs, conditions)}
        for entry in gene_data[seq_name]:
            entry['expression_levels'] = exp

# # #  Misc functions
def print_deficiencies(rxn_data):
    r_list = sorted(list(rxn_data.keys()))
    r_list.sort(key=lambda r:rxn_data[r]['deficiencies'])
    print('%i reactions with significant impact:' % len(r_list))
    for r_id in r_list:
        print('%s %.1f%% of objective value.' % (r_id, rxn_data[r_id]['objective']))
        print('\t%s' % rxn_data[r_id]['deficiencies'])
        print('\t%s' % ', '.join(rxn_data[r_id]['genes']))


# # #  I/O options
files_dir = '/mnt/hgfs/win_projects/brugia_project'
model_file = 'model_b_mal_4.5-wip.xlsx'
expression_file = 'All_Stages_Brugia_Wolbachia_FPKMs.xlsx'
expression_sheets = ('Brugia_FPKMs', 'Wolbachia_FPKMs')
gen_pept_file = 'utility/b_malayi_genpept.gp'
blast_xml_file = 'utility/model_b_mal_4.5-wip_single_kos_human_blast.xml'
gene_data_out_file = os.path.join(files_dir, 'bm_4.5_single_ko_gene_info.xlsx')
# # #  Intermediate files
prot_sequences_file = 'utility/model_b_mal_4.5-wip_single_ko_prots.fa'
blast_results_file = 'utility/x'
rxn_ko_data_file = 'utility/model_b_mal_4.5-wip_single_kos_rxns.pkl'
gene_ko_data_file = 'utility/model_b_mal_4.5-wip_single_kos_genes.pkl'
# # #  Run options
objective_threshold_fraction = 0.3 # Considered significant if resulting objective function is less than 0.3 (30%) of the original.
do_double_ko = False
expression_conditions = ['L3', 'L3D6', 'L3D9', 'L4', 'F30', 'M30', 'F42', 'M42', 'F120', 'M120']
expression_headings = [('Larval expression (L3|L3D6|L3D9|L4)', ('L3','L3D6','L3D9','L4')), ('Adult female expression (F30|F42|F120)', ('F30','F42','F120')), ('Adult male expression (M30|M42|M120)', ('M30','M42','M120'))]


# # #  Run steps
if not os.path.isfile(rxn_ko_data_file):
    rxn_data = {}
    model_path = os.path.join(files_dir, model_file)
    model = read_excel(model_path, verbose=False)
    rxn_to_genes = get_rxns_to_delete(model)
    do_deletions(rxn_data, model, rxn_to_genes, do_double_ko, objective_threshold_fraction) # Fills out 'objective', 'deficiencies', and 'genes' of reactions in rxn_data.
    save_data_object(rxn_data, rxn_ko_data_file)
else:
    rxn_data = load_data_object(rxn_ko_data_file)

#print_deficiencies(rxn_data)

if not os.path.isfile(gene_ko_data_file):
    gene_data = process_gene_data(rxn_data)
    get_expression_data(gene_data, os.path.join(files_dir, expression_file), expression_sheets, expression_conditions) # Fills out 'expression_levels'
    if not os.path.isfile(prot_sequences_file):
        prot_to_std = get_prot_name_translations(gene_data, gen_pept_file)
        prots = save_prot_sequences(gene_data, prot_to_std, prot_sequences_file)
    else:
        prots = sequ.loadfasta(prot_sequences_file)

    parse_blast_xml(gene_data, blast_xml_file)

    save_data_object(gene_data, gene_ko_data_file)
else:
    gene_data = load_data_object(gene_ko_data_file)

save_data_to_excel(gene_data, gene_data_out_file, expression_headings)
