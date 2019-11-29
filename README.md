Brugia project
==============

Unfortunately, these scripts were not designed to be user friendly. A brief summary of each:

* analysis_tools.py
A script made to be run directly. It is used to run a model under many different conditions, and to analyze and visualize the behaviors.
* analyze_metabolomics.py
This script was used to analyze the metabolomics data we received. None of the analyses really panned out.
* clean_model_metabolites.py
A one-time script made to be run directly. This parses the reactions and cleans up the annotations.
* compare_models.py
A one-time script made to be run directly. This was made to identify differences in constrains between life stage-specific models. It was never really finished.
* constrain_expression.py
A script made to be run directly. It takes a model file and RNA-seq data, and constrains the model reactions in a fairly simple manner. The maximum expression for one gene across all life stages is used to normalize all other life stages; that proportion is applied to the bounds of that reaction.
* get_knockout_info.py
A script made to be run directly. This script performs single- and double-knockouts, cross-references them against ChEMBL hits, life-stage specific RNA-seq expression data, and the human and Onchocerca volvulus proteins, then outputs all of the information into an Excel file. It expects that the relevant analyses have already been done.
* import_directionality.py
A one-time script made to be run directly. This compares the model against the iCEL1273 model, and imports the reaction directionality for shared reactions.
* import_mito.py
A one-time script made to be run directly. This compares the model reactions against those in the iCEL1273 model, and aids in developing the mitochondiral compartment.
* merge_wol_worm.py
A one-time script made to be run directly. This was used in conjunction with separate_wolbachia.py in order to properly compartmentalize the model.
* model_tools.py
A script made to be run directly, it also contains many useful importable functions. It is used to analyze models at a coarse level.
* read_excel.py
A library script made to be imported. This was created as at the time COBRAPY could not properly import models from excel files. I believe the current version of the library has remedied this. The script is also used to save a model back to excel.
* separate_wolbachia.py
A one-time script made to be run directly. This was used in conjunction with merge_wol_worm.py in order to properly compartmentalize the model.

There are also scripts in the utility directory:

* jitter_plots.py
Script to generate a combination jitter / bar plot. It was used to plot drug activity data, but is formulated for general data.
* translate_wBm_names.py
A one-time script used to modify the RNA-seq data we had for Wolbachia, to add the "Sequence name" column from the annotated genome file.
