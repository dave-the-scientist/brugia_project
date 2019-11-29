"""
Notes:
- Brugia protein sequences: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA10729
- wBm protein sequences: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=292805
- BLASTP against Reference proteins (refseq protein) from Human, using BLOSUM45 matrix.
- BLASTP against nr proteins from O. volvulus and wOv, using BLOSUM45 matrix.
  - Caution about the Oncho results; I'm not sure how many protein sequences have been annotated.
- The ChEMBL search results were performed under the "Target Search" tab on their website. Downloaded as a tab-deliminited file.
"""
import os, cPickle, pandas, re
from molecbio import sequ
from cobra.flux_analysis import single_reaction_deletion, double_reaction_deletion
from model_tools import load_model, id_bottleneck_metabolites
import xml.etree.ElementTree as ET


def get_rxns_to_delete(model):
    rxn_to_genes = {}
    for rxn in model.reactions:
        if not rxn.gene_names or not rxn.id.startswith(('R', 'ACYLCOA', 'N00001')):
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
    header_bg = '#DEEDED'
    sheet_name = 'Single knockouts'
    gene_header = 'Gene ID'
    headers_atts = [('# Reactions','num_reactions'), ('Reaction','reaction'), ('Associated genes','other_genes'), ('Objective %','objective'), ('Biomass deficiencies','deficiencies')]
    ortho_headers = ['Human homologs\n(#|% identity|% coverage)', 'O. volvulus homologs\n(#|% identity|% coverage)']
    chembl_headers = ['# ChEMBL hits', 'ChEMBL hits\n(% identity|species)']
    data = {h[0]:[] for h in headers_atts+expression_headings}
    for h in [gene_header] + ortho_headers + chembl_headers:
        data[h] = []
    gene_order = sorted(list(gene_data.keys()))
    gene_order.sort(key=lambda g:gene_data[g][0]['deficiencies'])
    for gene in gene_order:
        for g_data in gene_data[gene]:
            data[gene_header].append(gene)
            for h, att in headers_atts:
                data[h].append(g_data.get(att, 'NOT FOUND'))
            human_hlogs = '%i | %.1f | %.1f' % (g_data['num_human_prots'], g_data['human_prot_identity'],g_data['human_prot_coverage']) if g_data['num_human_prots'] else ' '
            data[ortho_headers[0]].append(human_hlogs)
            oncho_hlogs = '%i | %.1f | %.1f' % (g_data['num_oncho_prots'], g_data['oncho_prot_identity'],g_data['oncho_prot_coverage']) if g_data['num_oncho_prots'] else ' '
            data[ortho_headers[1]].append(oncho_hlogs)
            data[chembl_headers[0]].append(g_data.get('num_chembl_hits', 0))
            data[chembl_headers[1]].append(g_data.get('chembl_hits', ''))
            if '_max_observed_expression' in g_data['expression_levels']:
                max_expression = round(g_data['expression_levels']['_max_observed_expression'], 1)
            else:
                max_expression = " "
            data[expression_headings[0][0]].append(max_expression)
            for h, ls in expression_headings[1:]:
                exp_levels = [g_data['expression_levels'].get(l) for l in ls]
                data[h].append(' | '.join(exp_levels))
    col_headers = [gene_header] + [h[0] for h in headers_atts] + [i for i in ortho_headers+chembl_headers] + [j[0] for j in expression_headings]
    writer = pandas.ExcelWriter(gene_data_out_file, engine='xlsxwriter')
    df = pandas.DataFrame(data)[col_headers] # The [] specifies the order of the columns.
    df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
    worksheet = writer.sheets[sheet_name]
    header_format = writer.book.add_format({'bold': True, 'text_wrap': True, 'align': 'center', 'valign': 'top', 'bg_color': header_bg, 'border': 1})
    for i, h in enumerate(col_headers):
        col_w = max(len(line.strip()) for line in h.splitlines())
        col_width = max(col_w+1, min_column_width)
        if i in (0, 2, 3, 5, 9):
            col_format = writer.book.add_format({'align': 'left'})
        elif i == 10:
            col_format = writer.book.add_format({'align': 'center'})
        else:
            col_format = writer.book.add_format({'align': 'center'})
        worksheet.set_column(i, i, col_width, col_format)
        worksheet.write(0, i, h, header_format) # Header added manually.
    worksheet.freeze_panes(1, 0) # Freezes header row.
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
    sequ.savefasta(prots, prot_sequences_file, spaces=False, numbers=False)
    print('Saved %i sequences to %s' % (len(prots), prot_sequences_file))
    return prots

# # #  Parsing BLAST output
def parse_blast_xml(gene_data, blast_xml_file, taxon_name, spc_str):
    """taxon_name is used to name the properties saved in gene_data."""
    min_e_val = 1E-30
    property_strs = ['num_%s_prots', '%s_prot_id', '%s_prot_identity', '%s_prot_coverage']
    gi_split_regex = re.compile('\s?>gi\|\S+\|\S+\|\S+\|\s?')
    gene_spc_regex = re.compile('(.+) \[(.+)\]$')
    isoform_regex = re.compile('(.+) (isoform \S+)(.*)$')
    tree = ET.parse(blast_xml_file)
    root = tree.getroot()
    iterations = root.find('BlastOutput_iterations')
    for q_hit in iterations:
        gene = q_hit.find('Iteration_query-def').text
        if gene not in gene_data:
            continue
        prot_len = float(q_hit.find('Iteration_query-len').text)
        s_hits = q_hit.find('Iteration_hits')
        hit_names, top_hit_id, top_e_val, top_identity, top_coverage = get_good_hits(s_hits, min_e_val, spc_str.lower(), gi_split_regex, gene_spc_regex, isoform_regex)
        num_hits = len(hit_names)
        top_coverage = round(top_coverage/prot_len*100.0, 1)
        for g_data in gene_data[gene]:
            for p_str, val in zip(property_strs, [num_hits, top_hit_id, top_identity, top_coverage]):
                g_data[p_str % taxon_name] = val


def get_good_hits(s_hits, min_e_val, spc_str, gi_split_regex, gene_spc_regex, isoform_regex):
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
            name = parse_name_from_hit(s_hit, spc_str, gi_split_regex, gene_spc_regex, isoform_regex)
            if not name:
                continue # A hit was found, but it did not match the spc_str
            hit_names.add(name)
            if hit_e_val < best_e_val:
                best_hit_id = s_hit.find('Hit_accession').text.strip()
                best_ident = hit_ident
                best_e_val, best_coverage = hit_e_val, hit_coverage
    if not hit_names:
        return hit_names, None, None, 0, 0
    return hit_names, best_hit_id, best_e_val, best_ident, best_coverage
def parse_name_from_hit(s_hit, spc_str, gi_split_regex, gene_spc_regex, isoform_regex):
    name = find_gene_from_species(s_hit, spc_str, gi_split_regex, gene_spc_regex)
    if not name:
        return False
    if 'isoform' in name:
        nm, iso, rem = isoform_regex.match(name).groups()
        name = nm + rem
    if name.lower().startswith('predicted: '):
        name = name[11:]
    return name

def find_gene_from_species(s_hit, spc_str, gi_split_regex, gene_spc_regex):
    for hit in gi_split_regex.split( s_hit.find('Hit_def').text ):
        m = gene_spc_regex.match(hit)
        if not m:
            continue
        name, spc = m.groups()
        if spc_str in spc.lower():
            return name
    return False

# # #  Getting expression data
def get_expression_data(gene_data, expression_file, sheetnames, conditions):
    for sheetname in sheetnames:
        parse_expression_sheet(gene_data, expression_file, sheetname, conditions)
    null_exp = {c:'--' for c in conditions}
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
        max_expression = max(avgs)
        exp = {c:'%i'%(round(a/max_expression*100.0, 0) if max_expression else 0) for a,c in zip(avgs, conditions)}
        exp['_max_observed_expression'] = max_expression
        for entry in gene_data[seq_name]:
            entry['expression_levels'] = exp

# # #  Parse ChEMBL search file
def parse_chembl_results(gene_data, chembl_results_file):
    max_e_val = 1E-30
    chembl_data = {}
    total_hits, sig_hits = 0, 0
    with open(chembl_results_file) as f:
        f.readline() # Header
        for line in f:
            if not line.strip():
                continue
            total_hits += 1
            gene, chembl_id, tid, description, uniprot_id, target_type, species, _, _, identity, blast_score, e_value = line.split('\t')
            identity, e_value = float(identity), float(e_value)
            if e_value > max_e_val:
                continue
            sig_hits += 1
            hit_data = {'chembl_id':chembl_id, 'species':species, 'identity':identity, 'e_value':e_value}
            chembl_data.setdefault(gene, []).append(hit_data)
    print('%i of the %i ChEMBL hits were below the E-value threshold of %.1e' % (sig_hits, total_hits, max_e_val))
    for gene, data_list in chembl_data.items():
        if gene not in gene_data:
            continue
        data_list.sort(key=lambda d: d['e_value'])
        chembl_hits = ', '.join('%s (%i | %s)' % (d['chembl_id'], round(d['identity'], 0), d['species']) for d in data_list)
        for g_data in gene_data[gene]:
            g_data['num_chembl_hits'] = len(data_list)
            g_data['chembl_hits'] = chembl_hits


# # #  Misc functions
def print_deficiencies(rxn_data):
    r_list = sorted(list(rxn_data.keys()))
    r_list.sort(key=lambda r:rxn_data[r]['deficiencies'])
    print('%i reactions with significant impact:' % len(r_list))
    for r_id in r_list:
        print('%s %.1f%% of objective value.' % (r_id, rxn_data[r_id]['objective']))
        print('\t%s' % rxn_data[r_id]['deficiencies'])
        print('\t%s' % ', '.join(rxn_data[r_id]['genes']))


# # #  Main paths
files_dir = '/mnt/hgfs/win_projects/brugia_project'
utility_dir = '/home/dave/Desktop/projects/brugia_project/utility'

# # #  Main run options
model_file = 'model_b_mal_4.5-wip.xlsx'
run_str = 'bm_4.5-lo_ox-lo_glu'
wolbachia_ratio = 0.1
objective_threshold_fraction = 0.25 # Considered significant if resulting objective function is less than 0.25 (25%) of the original.
do_double_ko = False
expression_conditions = ['L3', 'L3D6', 'L3D9', 'L4', 'F30', 'M30', 'F42', 'M42', 'F120', 'M120']
expression_headings = [('Max\nexpression',), ('Larval expression\n(L3|L3D6|L3D9|L4)', ('L3','L3D6','L3D9','L4')), ('Adult female expression\n(F30|F42|F120)', ('F30','F42','F120')), ('Adult male expression\n(M30|M42|M120)', ('M30','M42','M120'))]
gene_data_out_file = os.path.join(files_dir, '%s_gene_info.xlsx'%(run_str))

# # #  Required files
expression_file = os.path.join(files_dir, 'All_Stages_Brugia_Wolbachia_FPKMs.xlsx')
expression_sheets = ('Brugia_FPKMs', 'Wolbachia_FPKMs')
gen_pept_file = os.path.join(utility_dir, 'b_malayi_genpept.gp')
human_blast_xml_file = os.path.join(utility_dir, '%s_human_blast.xml'%(run_str))
oncho_blast_xml_file = os.path.join(utility_dir, '%s_oncho_blast.xml'%(run_str))
chembl_results_file = os.path.join(utility_dir, '%s_chembl.txt'%(run_str))

# # #  Intermediate files created
prot_sequences_file = os.path.join(utility_dir, '%s_prots.fa'%(run_str))
rxn_ko_data_file = os.path.join(utility_dir, '%s_rxns.pkl'%(run_str))
gene_ko_data_file = os.path.join(utility_dir, '%s_genes.pkl'%(run_str))


# # #  Run steps
if not os.path.isfile(rxn_ko_data_file):
    rxn_data = {}
    model_path = os.path.join(files_dir, model_file)
    model = load_model(model_path, wolbachia_ratio)
    rxn_to_genes = get_rxns_to_delete(model)
    do_deletions(rxn_data, model, rxn_to_genes, do_double_ko, objective_threshold_fraction) # Fills out 'objective', 'deficiencies', and 'genes' of reactions in rxn_data.
    save_data_object(rxn_data, rxn_ko_data_file)
else:
    rxn_data = load_data_object(rxn_ko_data_file)

#print_deficiencies(rxn_data)

if not os.path.isfile(gene_ko_data_file):
    gene_data = process_gene_data(rxn_data)
    get_expression_data(gene_data, expression_file, expression_sheets, expression_conditions) # Fills out 'expression_levels'
    if not os.path.isfile(prot_sequences_file):
        prot_to_std = get_prot_name_translations(gene_data, gen_pept_file)
        prots = save_prot_sequences(gene_data, prot_to_std, prot_sequences_file)
    else:
        prots = sequ.loadfasta(prot_sequences_file)
    for blast_file in [human_blast_xml_file, oncho_blast_xml_file]:
        if not os.path.isfile(blast_file):
            print('Error: no BLAST results found at %s' % blast_file)
            exit()
    parse_blast_xml(gene_data, human_blast_xml_file, 'human', 'homo sapiens')
    parse_blast_xml(gene_data, oncho_blast_xml_file, 'oncho', 'onchocerca volvulus')
    if not os.path.isfile(chembl_results_file):
        print('Error: no ChEMBL results found at %s' % chembl_results_file)
        exit()
    # parse_chembl_results(gene_data, chembl_results_file) # Where it should be called.
    save_data_object(gene_data, gene_ko_data_file)
else:
    gene_data = load_data_object(gene_ko_data_file)

parse_chembl_results(gene_data, chembl_results_file) # # # Temp place to be called from.

save_data_to_excel(gene_data, gene_data_out_file, expression_headings)
