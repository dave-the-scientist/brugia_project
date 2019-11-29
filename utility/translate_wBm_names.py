"""The Brugia FPKM file did not have a Sequence Name column, which would indicate the locus name for each CDS listed in the Gene column. This script generates that column.
NOTE: Reads the whole genpept file as a string, so don't use with absurdly large files.
"""
import os, pandas

def parse_locus_cds(ent):
    locus, cds = None, None
    for line in ent.splitlines():
        if line.startswith('LOCUS'):
            if locus != None:
                print('\nError: multiple "LOCUS" lines identified in the following entry:\n%s' % ent)
                return False, False
            locus = line.split()[1]
        elif '/locus_tag=' in line:
            if cds != None:
                print('\nError: multiple "/locus_tag" lines identified in the following entry:\n%s' % ent)
                return False, False
            cds = line.partition('=')[2].strip()[1:-1]
    if locus == None:
        print('\nError: could not parse a Locus from the following entry:\n%s' % ent)
        return False, False
    if cds == None:
        print('\nError: could not parse a CDS from the following entry:\n%s' % ent)
        return False, False
    return locus, cds

def parse_genpept(filename):
    cds_trans = {}
    with open(filename) as f:
        f_str = f.read()
    gp_entries = [e.strip() for e in f_str.split('//') if not e.isspace()]
    print('Parsed %i entries.' % len(gp_entries))
    for ent in gp_entries:
        locus, cds = parse_locus_cds(ent)
        if locus == cds == False:
            exit()
        elif cds_trans.get(cds, locus) != locus:
            print('Error: CDS %s was associated with more than 1 locus.' % cds)
            exit()
        cds_trans[cds] = locus
    return cds_trans

def generate_trans_col(cds_trans, filename, sheet, column_header):
    column_data = []
    frame = pandas.read_excel(filename, sheet)
    for i in frame.index:
        row = frame.ix[i]
        cds = row[column_header]
        column_data.append(cds_trans.get(cds, ''))
    return column_data


# # #  Options
expression_file = '/mnt/hgfs/win_projects/brugia_project/All_Stages_Brugia_Wolbachia_FPKMs.xlsx'
sheet = 'Wolbachia_FPKMs'
column_header = 'Gene'
genpept_file = 'wBm_genpept.gp'
column_out_file = 'wBm_cds_gene_translation.txt'

# # #  Run steps
cds_trans = parse_genpept(genpept_file)
column_data = generate_trans_col(cds_trans, expression_file, sheet, column_header)
with open(column_out_file, 'w') as f:
    f.write('\n'.join(column_data))
print('Wrote data for %i rows to %s' % (len(column_data), os.path.realpath(column_out_file)))
