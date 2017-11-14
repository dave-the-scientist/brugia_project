# Downloaded from the openCobra github at https://github.com/opencobra, under the m_model_collection project.
# I have edited it to work with our lab's existing file formats of network models.
from warnings import warn

from six import string_types, iteritems
from math import isnan
import pandas
import os

from cobra import Model, Metabolite, Reaction

# Potential string keys in the excel files
RXN_SHEET_NAMES = {"reaction list", "reactions"}
MET_SHEET_NAMES = {"metabolite list", "metabolites"}
MET_ID_KEYS = {"abbreviation", "abbreviations", "metabolite abbreviation",
    "machine readable id"}
MET_NAME_KEYS = {"metabolite name", "officialname", "metabolite_name", "description",
    "human readable id"}
MET_COMPARTMENT_KEYS = {"compartment"}
MET_FORMULA_KEYS = {"chemical formula", "formula"}
RXN_ID_KEYS = {"abbreviation", "reaction #", "abbrev",
    "reaction id", "reaction abbreviation", "id"}
RXN_NAME_KEYS = {"name", "rxn description", "description"}
RXN_STR_KEYS = {"equation", "reaction formula", "reaction", "machine readable"}
RXN_GPR_KEYS = {"geneassociation", "gpr"}
RXN_LB_KEYS = {"lb", "lower bound", "lower bounds"}
RXN_UB_KEYS = {"ub", "upper bound", "upper bounds"}
RXN_OBJ_KEYS = {"objective"}
RXN_SUBS_KEYS = {"subsystem", "pathway"}
RXN_ECS_KEYS = {"ec number", "enzyme"}
RXN_GENES_KEYS = {"genes"}
RXN_PROTS_KEYS = {"proteins"}
RXN_NOTES_KEYS = {"notes", "comments"}
RXN_CONFD_KEYS = {"confidence score"}
RXN_KEGG_KEYS = {"kegg"}


def escape_str(potential_str):
    if isinstance(potential_str, bytes):
        return potential_str
    if hasattr(potential_str, "encode"):
        return potential_str.encode("ascii", "replace")
    return bytes(potential_str)


def guess_name(potential_names, allowed_names, fail=True):
    """see if any of the potential names are allowed for a category"""
    for potential_name in potential_names:
        if potential_name.lower() in allowed_names:
            return potential_name
    if fail:
        raise ValueError("could not find any of %s in %s" %
                         (str(set(allowed_names)),
                          str(set(potential_names))))
    else:
        return None


def extract(row, keyname, type=str):
    """extract a value which may be missing, for a potentially missing key"""
    if keyname is None or keyname == "skip":
        return type()
    value = row[keyname]
    if isinstance(value, float) and isnan(value):
        return type()
    elif type is str:
        return escape_str(value).strip()
    else:
        return type(value)


def read_excel(
        filename,
        verbose=True,
        rxn_sheet_name=None,
        rxn_sheet_header=0,
        rxn_skip_rows=set(),
        rxn_sheet_converters=None,

        rxn_id_key=None,
        rxn_name_key=None,
        rxn_str_key=None,
        rxn_gpr_key=None,
        rxn_lb_key=None,
        rxn_ub_key=None,
        rxn_fwd_arrow=None,
        rxn_rev_arrow=None,
        rxn_reversible_arrow=None,
        rxn_objective_key=None,
        rxn_subsystem_key=None,
        rxn_ec_num_key=None,
        rxn_genes_key=None,
        rxn_proteins_key=None,
        rxn_confidence_key=None,
        rxn_notes_key=None,
        rxn_kegg_key=None,

        met_sheet_name=None,
        met_sheet_header=0,
        met_id_key=None,
        met_name_key=None,
        met_compartment_key=None,
        met_formula_key=None,

        ):

    # autodetect sheet names
    pio = pandas.io.excel.ExcelFile(filename)
    sheet_names = pio.sheet_names
    if rxn_sheet_name is None:
        # only one sheet means it must be that one
        if len(sheet_names) == 1:
            rxn_sheet_name = sheet_names[0]
        else:
            rxn_sheet_name = guess_name(pio.sheet_names, RXN_SHEET_NAMES)
    if met_sheet_name is None:
        met_sheet_name = guess_name(
            pio.sheet_names, MET_SHEET_NAMES, fail=False)

    # Begin model creation
    try:
        model_id = os.path.splitext(os.path.split(filename)[1])[0]
    except:
        model_id = "imported_model"
    m = Model(model_id)

    # Metabolites - if the sheet is found
    met_renames = {}
    if met_sheet_name is not None and met_sheet_name != "ignore":
        met_frame = pandas.read_excel(filename, met_sheet_name,
                                      header=met_sheet_header)
        # strip spaces from header
        met_frame.columns = [i.strip() for i in met_frame.keys()]
        if met_id_key is None:
            met_id_key = guess_name(met_frame.keys(), MET_ID_KEYS)
        if met_name_key is None:
            met_name_key = guess_name(met_frame.keys(), MET_NAME_KEYS, fail=False)
        if met_compartment_key is None:
            met_compartment_key = guess_name(met_frame.keys(), MET_COMPARTMENT_KEYS, fail=False)
        if met_formula_key is None:
            met_formula_key = guess_name(met_frame.keys(), MET_FORMULA_KEYS, fail=False)
        for i in met_frame.index:
            met_row = met_frame.ix[i]
            met_attributes = {}
            if met_compartment_key is not None:
                compartment = extract(met_row, met_compartment_key)
            if met_formula_key is not None:
                formula = extract(met_row, met_formula_key)
                if formula is not None and formula.lower() != "None":
                    met_attributes["formula"] = formula
            met_id = extract(met_row, met_id_key)
            if len(met_id) == 0:
                continue
            if " " in met_id:
                new_id = met_id.replace(" ", "_")
                met_renames[met_id] = new_id
                if verbose:
                    print("Renamed metabolite '%s' to '%s'" % (met_id, new_id))
                met_id = new_id
            if met_name_key is not None:
                met_name = extract(met_row, met_name_key)
            else:
                met_name = extract(met_row, met_id_key)
            met = Metabolite(met_id, name=met_name, compartment=compartment, **met_attributes)

            try:
                m.add_metabolites(met)
            except ValueError:
                if verbose:
                    print("duplicate metabolite '%s' not added" % met.id)
    elif verbose:
        met_frame = None
        print("metabolite sheet not found")
    # need to rename longer strings first, then shorter strings
    met_rename_list = list(sorted((iteritems(met_renames)),
                                  key=lambda x: len(x[0]), reverse=True))

    # Reactions
    rxn_frame = pandas.read_excel(filename, rxn_sheet_name,
                                  header=rxn_sheet_header,
                                  skiprows=rxn_skip_rows,
                                  converters=rxn_sheet_converters)
    # strip spaces from header
    rxn_frame.columns = [i.strip() for i in rxn_frame.keys()]
    if rxn_id_key is None:
        rxn_id_key = guess_name(rxn_frame.keys(), RXN_ID_KEYS)
    if rxn_str_key is None:
        rxn_str_key = guess_name(rxn_frame.keys(), RXN_STR_KEYS)
    if rxn_name_key is None:
        rxn_name_key = guess_name(rxn_frame.keys(), RXN_NAME_KEYS, fail=False)
        if verbose and rxn_name_key is None:
            print("reaction name column not identified")
    if rxn_gpr_key is None:
        rxn_gpr_key = guess_name(rxn_frame.keys(), RXN_GPR_KEYS, fail=False)
        if verbose and rxn_gpr_key is None:
            print("gene reaction rule column not identified")
    if rxn_lb_key is None:
        rxn_lb_key = guess_name(rxn_frame.keys(), RXN_LB_KEYS, fail=False)
        if verbose and rxn_lb_key is None:
            print("reaction lower bound column not identified")
    if rxn_ub_key is None:
        rxn_ub_key = guess_name(rxn_frame.keys(), RXN_UB_KEYS, fail=False)
        if verbose and rxn_ub_key is None:
            print("reaction upper bound column not identified")
    if rxn_objective_key is None:
        rxn_objective_key = guess_name(rxn_frame.keys(), RXN_OBJ_KEYS, fail=False)
        if verbose and rxn_objective_key is None:
            print("reaction objective function column not identified")
    if rxn_subsystem_key is None:
        rxn_subsystem_key = guess_name(rxn_frame.keys(), RXN_SUBS_KEYS, fail=False)
        if verbose and rxn_subsystem_key is None:
            print("reaction subsystem column not identified")
    if rxn_ec_num_key is None:
        rxn_ec_num_key = guess_name(rxn_frame.keys(), RXN_ECS_KEYS, fail=False)
        if verbose and rxn_ec_num_key is None:
            print("reaction EC number column not identified")
    if rxn_genes_key is None:
        rxn_genes_key = guess_name(rxn_frame.keys(), RXN_GENES_KEYS, fail=False)
        if verbose and rxn_genes_key is None:
            print("reaction gene names column not identified")
    if rxn_proteins_key is None:
        rxn_proteins_key = guess_name(rxn_frame.keys(), RXN_PROTS_KEYS, fail=False)
        if verbose and rxn_proteins_key is None:
            print("reaction protein names column not identified")
    if rxn_confidence_key is None:
        rxn_confidence_key = guess_name(rxn_frame.keys(), RXN_CONFD_KEYS, fail=False)
        if verbose and rxn_confidence_key is None:
            print("reaction confidence score column not identified")
    if rxn_notes_key is None:
        rxn_notes_key = guess_name(rxn_frame.keys(), RXN_NOTES_KEYS, fail=False)
        if verbose and rxn_notes_key is None:
            print("reaction notes column not identified")
    if rxn_kegg_key is None:
        rxn_kegg_key = guess_name(rxn_frame.keys(), RXN_KEGG_KEYS, fail=False)
        if verbose and rxn_kegg_key is None:
            print("KEGG ID column not identified")

    for i in range(len(rxn_frame)):
        row = rxn_frame.ix[i]
        if rxn_id_key == "auto":
            rxn_id = "R%04d" % i
        else:
            rxn_id = extract(row, rxn_id_key)
        rxn_str = extract(row, rxn_str_key)
        # skip empty rows
        if not isinstance(rxn_id, string_types) or \
                not isinstance(rxn_str, string_types) or \
                len(rxn_str) == 0 or len(rxn_id) == 0:
            continue
        # skip duplicated header rows
        if rxn_id == rxn_id_key and rxn_str == rxn_str_key:
            continue
        rxn = Reaction()
        rxn.id = rxn_id
        rxn.name = extract(row, rxn_name_key)
        if rxn.id in m.reactions:
            if verbose:
                print("duplicate reaction '%s' found" % rxn.id)
            # add underscores to the end of duplicate reactions
            while rxn.id in m.reactions:
                rxn.id += "_"
        # fill out descriptive attributes
        if rxn_objective_key is not None:
            rxn.objective_coefficient = extract(row, rxn_objective_key, float)
        if rxn_subsystem_key is not None:
            rxn.subsystem = extract(row, rxn_subsystem_key)
        if rxn_ec_num_key is not None:
            rxn.enzyme_commission = extract(row, rxn_ec_num_key)
        if rxn_genes_key is not None:
            rxn.gene_names = extract(row, rxn_genes_key)
        if rxn_proteins_key is not None:
            rxn.protein_names = extract(row, rxn_proteins_key)
        if rxn_confidence_key is not None:
            rxn.confidence_notes = extract(row, rxn_confidence_key)
        if rxn_notes_key is not None:
            rxn.reaction_notes = extract(row, rxn_notes_key)
        if rxn_kegg_key is not None:
            rxn.kegg_reaction = extract(row, rxn_kegg_key)
        m.add_reaction(rxn)

        # Now build the reaction from the string

        # no need to print new metabolite created if no metaboltie sheet
        verbose_build = verbose and met_frame is not None
        build_kwargs = {"verbose": verbose, "fwd_arrow": rxn_fwd_arrow,
                        "rev_arrow": rxn_rev_arrow,
                        "reversible_arrow": rxn_reversible_arrow}
        try:
            rxn.build_reaction_from_string(rxn_str, **build_kwargs)
        except Exception as e:
            # replace metabolites which have spaces, and try again
            fixed_rxn_str = rxn_str
            for k, v in met_rename_list:
                fixed_rxn_str = fixed_rxn_str.replace(k, v)
            if verbose:
                print("converted '%s' to '%s'" % (rxn_str, fixed_rxn_str))
            try:
                rxn.build_reaction_from_string(fixed_rxn_str, **build_kwargs)
            except Exception as e:
                print("Error parsing %s string '%s'" % (repr(rxn), rxn_str))
                raise e

        # parse gene reaction rule
        gpr = extract(row, rxn_gpr_key)
        if "," in gpr or "+" in gpr:
            # break on ',' (which is or) then '+' (which is and)
            ors = ["( %s  )" % o.replace("+", " and ") if "+" in o else o
                   for o in gpr.split(",")]
            gpr = " or ".join(ors)
        rxn.gene_reaction_rule = gpr
        if rxn_lb_key is not None:
            rxn.lower_bound = float(row[rxn_lb_key])
        if rxn_ub_key is not None:
            rxn.upper_bound = float(row[rxn_ub_key])
        # #  Use default values, with reversibility from rxn.

    # fix upper and lower bounds if they include infinity
    inf = float("inf")
    max_lower_bound = max(abs(i) for i in m.reactions.list_attr("lower_bound")
                          if abs(i) < inf)
    max_upper_bound = max(abs(i) for i in m.reactions.list_attr("upper_bound")
                          if abs(i) < inf)
    inf_replace = max(max_upper_bound, max_lower_bound, 1000.) * 100

    def clip_inf(value):
        if value == inf:
            return inf_replace
        elif value == -inf:
            return -inf_replace
        else:
            return value
    for reaction in m.reactions:
        reaction.lower_bound = clip_inf(reaction.lower_bound)
        reaction.upper_bound = clip_inf(reaction.upper_bound)

    return m

def write_excel(model, outfile):
    rxn_sheet = 'Reaction List'
    rxn_headers = [('Abbreviation','id'), ('Description', 'name'), ('Reaction','reaction'),
        ('GPR', 'gene_reaction_rule'), ('Genes', 'gene_names'), ('Proteins', 'protein_names'),
        ('Subsystem', 'subsystem'), ('Reversible', 'reversibility'),
        ('Lower bound', 'lower_bound'), ('Upper bound', 'upper_bound'),
        ('Objective', 'objective_coefficient'), ('Confidence Score', 'confidence_notes'),
        ('EC Number', 'enzyme_commission'), ('Notes', 'reaction_notes')]
    mtb_sheet = 'Metabolite List'
    mtb_headers = [('Abbreviation', 'id'), ('Description', 'name'), ('Neutral formula', ''),
        ('Charged formula', ''), ('Charge', ''), ('Compartment', 'compartment'), ('KEGG ID', ''),
        ('PubChem ID', ''), ('ChEBI ID', ''), ('InChI string', ''), ('SMILES', ''),
        ('HMDB ID', '')]
    writer = pandas.ExcelWriter(outfile, engine='xlsxwriter')
    rxn_data = generate_sheet_data(model.reactions, rxn_headers)
    rxn_df = generate_dataframe(rxn_sheet, rxn_headers, rxn_data, writer)
    mtb_data = generate_sheet_data(model.metabolites, mtb_headers)
    mtb_df = generate_dataframe(mtb_sheet, mtb_headers, mtb_data, writer)
    writer.save()
def generate_sheet_data(obj_list, headers, id_sort=True):
    """obj_list should be a list of cobra Metabolite or Reaction objects."""
    if id_sort:
        obj_list = sorted(obj_list, key=lambda o: o.id)
    data = {h[0]:[] for h in headers}
    for o in obj_list:
        for h, att in headers:
            data[h].append(getattr(o, att, ''))
    return data
def generate_dataframe(sheet_name, headers, data, writer):
    min_column_width = 6
    col_headers = [h[0] for h in headers]
    df = pandas.DataFrame(data)[col_headers] # The [] specifies the order of the columns.
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    worksheet = writer.sheets[sheet_name]
    for i, h in enumerate(col_headers):
        col_width = max(len(h)+1, min_column_width)
        worksheet.set_column(i, i, col_width)
    worksheet.freeze_panes(1, 0) # Freezes top row.
    return df
