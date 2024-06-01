###################################
# Script to summarize all modeling
# information into one table
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
###################################

import pandas as pd
import glob
import os
import sys
import numpy as np

sys.path.append('../utils')
from create_summary_table import *
import utils

###########################
# Read files
###########################
modeling_script = '../scripts/modeling.py'
mmcif_file = '../scripts/IM_MuSIC_cluster.cif'

analysis_dir = '../results/'
clustering_dir = os.path.join(analysis_dir,'clustering')
rmf3 = os.path.join(clustering_dir,'cluster.0','cluster_center_model.rmf3')

I = get_input_information(mmcif_file)
input_information = I.get_dictionaries()


R = get_representation(clustering_dir)
representation = R.get_dictionaries()

S = read_modeling_information(modeling_script,
                              analysis_dir,
                              clustering_dir)


sampling = S.get_dictionaries_sampling()


samples = S.get_dictionaries_models()

clustering = S.get_dictionaries_clustering()

#S.update_mmcif_file(mmcif_file)

V = read_validation_information(clustering_dir)
validation = V.get_dictionaries()

#V = read_benchmark_information(clustering_dir)
#benchmark = V.get_dictionaries()

SS = get_software_information(mmcif_file)
software = SS.get_dictionaries()

D = get_data_availability(clustering_dir)
data_availability = D.get_dictionaries()

################################################
# Edit dictionaries
# Entries is dictionaries can be edited to add
# other custom information
################################################
input_information['Input information'] = []
input_information['Input information'].append('Atomic structure prediction from Alphafold2 (AF\_Q96GS4, AF\_Q9BY32)')
input_information['Input information'].append('Comparative model of the SLC38A9-RagA-RagB-Ragulator complex (template PDB ID 6WJ2)')
input_information['Input information'].append('AlphaFold-Multimer pairwise predictions of BORCS6 or ITPA and other components')

print(representation.keys())

print(representation['Rigid body (RB) definitions'])

representation['Rigid body (RB) definitions'] = ['RB1: Lamtor1$_{49-68}$,Lamtor1$_{83-146}$,Lamtor2$_{5-125}$,Lamtor3$_{1-124}$,Lamtor4$_{1-99}$, \n Lamtor5$_{1-91}$,RRAGB$_{1-337}$,RRAGB$_{344-374}$,RRAGC$_{62-369}$,SLC38A9$_{39-78}$, \n SLC38A9$_{80-96}$',
                                                 'RB2: BORCS6$_{199-216}$',
                                                 'RB3: BORCS6$_{221-228}$',
                                                 'RB4: BORCS6$_{258-295}$',
                                                 'RB5: BORCS6$_{319-355}$',
                                                 'RB6: ITPA$_{1-194}$']



representation['Resolution of structured components'] = '1 [R1] residue per bead'
representation['Resolution of disordered regions'] = '1 [R1] residues per bead'


representation['Spatial restraints encoded into scoring function'] = representation.pop('Spatial restraints encoded into scoring function')

representation['Spatial restraints encoded into scoring function'].append('Bayesian binary binding Mode restraint')

sampling['CPU time'] = ['5 hours on 80 processors']


validation['Percent of sequence connectivity restraints satisfied per structure'] = ['99 \%']
validation['Percent of excluded volume restraints satisfied per structure'] = ['99 \%']



print(samples)

###################
print('-------------')
print(clustering)

sampling['Replica exchange temperature range'] =  ['1.0 - 4.0']

software['Modeling scripts'] = ['https://github.com/integrativemodeling/IM\_MuSIC']
software['Structure prediction'] = ['AlphaFold2', 'AlphaFold-Multimer']
software['Visualization and plotting'] = ['UCSF Chimera']

################################################
# Convert ordered dictionaries 
# into lists
################################################
input_information_list = dict_to_list(input_information)
representation_list = dict_to_list(representation)
sampling_list = dict_to_list(sampling)
samples_list = dict_to_list(samples)
clustering_list = dict_to_list(clustering)
validation_list = dict_to_list(validation)
software_list = dict_to_list(software)
data_availability_list = dict_to_list(data_availability)

print(sampling_list)


################################################
# Compile all information
# 
################################################
variable_dict = {'complex': 'RAG-Ragulator MuSIC system',
                 'number':1,
                 'input_information': input_information_list, 
                 'representation': representation_list,
                 'sampling': sampling_list,
                 'samples': samples_list,
                 'clustering':clustering_list,
                 'validation':validation_list,
                 'software':software_list}

################################################
# Generate tex, pdf file
################################################
template = utils.get_template('../utils/SI_template.tex')
utils.compile_pdf_from_template(template, variable_dict, './table_SI_IM.pdf')

exit()
