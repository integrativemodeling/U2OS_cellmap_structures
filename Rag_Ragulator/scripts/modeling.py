'''
Script for integrative modeling of MuSIC2
cluster using a topology file
'''
# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.dof
import IMP.pmi.mmcif

#import IMP.pmi.restraints.RestraintBase

#import IMP.bbm
#import IMP.bbm.restraints

import os
import sys
    
num_frames = 60000 
num_mc_steps = 10 

if '--mmcif' in sys.argv:
    num_frames = 5

## Define Input Files  ####


top_dir = "../" 
datadirectory= f"{top_dir}/data/"

# Topology File
topology_file = f"{top_dir}/scripts/topology.txt"


# Initialize model
m=IMP.Model()

#Read in the topology file
topology= IMP.pmi.topology.TopologyReader(topology_file,
                                          pdb_dir = f'{top_dir}/data/pdb',
                                          fasta_dir = f'{top_dir}/data/fasta')

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m,resolutions=[1, 5])

#Each state can be specified by a topology file
s = bs.add_state(topology)


# Generate mmcif file
    
if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput()
    po.system.title = ('Integrative structure of the Rag-Ragulator MuSIC assembly')
    bs.system.add_protocol_output(po)

#Building the system with default rotation and translation parameters
root_hier, dof = bs.execute_macro()
                                  
#Setting up restraints

outputobjects = []

# Connectivity keeps things connected along the backbone (ignores if inside
# same rigid body)
mols=[]
crs=[]
moldict = bs.get_molecules()[0]

for molname in moldict:
    for mol in moldict[molname]:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.add_to_model()
        outputobjects.append(cr)
        crs.append(cr)
        mols.append(mol)

    
#Excluded volume restraints
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=5)
ev.add_to_model()
outputobjects.append(ev)

#External Barrier Restraints
#ex_barrier = IMP.pmi.restraints.RestraintBase.ExternalBarrier(root_hier,radius=)
#ex_barrier.add_to_model()
#outputobjects.append(ex_barrier)

'''
#AlphaFold Restraints

params_199_228 = IMP.bbm.restraints.GaussianVMFMixtureParameters.from_file('/wynton/home/sali/ignacia/CCMI/RAG_Ragulator/AF2/BORCS6_C_Terminus_LAMTOR2_Q9Y2Q5/param/param_A_199_217_params.pkl')
prb1 = IMP.bbm.restraints.GaussianVMFMixturePairRestraint(root_hier,
                                                         subunits = ((199,217,'BORCS6',0), 'Lamtor2'),
                                                         label='BORCS6_199_217',
                                                         params=params_199_228)
prb1.add_to_model()
prb1.set_weight(1.0)
outputobjects.append(prb1)



params_258_295 = IMP.bbm.restraints.GaussianVMFMixtureParameters.from_file('/wynton/home/sali/ignacia/CCMI/RAG_Ragulator/AF2/BORCS6_C_Terminus_LAMTOR2_Q9Y2Q5/param/param_A_258_295_params.pkl')
prb2 = IMP.bbm.restraints.GaussianVMFMixturePairRestraint(root_hier,
                                                         subunits = ((258,295,'BORCS6',0), 'Lamtor2'),
                                                         label='BORCS6_258_295',
                                                         params=params_258_295)
prb2.add_to_model()
prb2.set_weight(1.0)
outputobjects.append(prb2)

params_319_355 = IMP.bbm.restraints.GaussianVMFMixtureParameters.from_file('/wynton/home/sali/ignacia/CCMI/RAG_Ragulator/AF2/BORCS6_C_Terminus_LAMTOR2_Q9Y2Q5/param/param_A_319_355_params.pkl')
prb3 = IMP.bbm.restraints.GaussianVMFMixturePairRestraint(root_hier,
                                                         subunits = ((319,355,'BORCS6',0), (258,295,'BORCS6',0)),
                                                         label='BORCS6_319_355',
                                                         params=params_319_355)
prb3.add_to_model()
prb3.set_weight(1.0)
outputobjects.append(prb3)

params_ITPA = IMP.bbm.restraints.GaussianVMFMixtureParameters.from_file('/wynton/home/sali/ignacia/CCMI/RAG_Ragulator/AF2/ITPA_LAMTOR5_6WJ2_E/param/param_ITPA_params.pkl')
prb4 = IMP.bbm.restraints.GaussianVMFMixturePairRestraint(root_hier,
                                                         subunits = ('ITPA', 'Lamtor5'),
                                                         label='ITPA',
                                                         params=params_ITPA)
prb4.add_to_model()
prb4.set_weight(1.0)
outputobjects.append(prb4)

'''
# Save ref rmf 
output = IMP.pmi.output.Output()
output.init_rmf(f"ini_ranked_0.rmf3", [root_hier])
output.write_rmf(f"ini_ranked_0.rmf3")
output.close_rmf(f"ini_ranked_0.rmf3")


#####################################################
#                      SAMPLING                     #
#####################################################

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(root_hier)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(100)

#--------------------------
# Set MC Sampling Parameters
#---------------------------

mc1=IMP.pmi.macros.ReplicaExchange(m,
                                   root_hier=root_hier,
                                   monte_carlo_sample_objects=dof.get_movers(),
                                   output_objects=outputobjects,
                                   number_of_best_scoring_models=0,
                                   monte_carlo_steps=num_mc_steps,
                                   number_of_frames=num_frames,
                                   global_output_directory="output/")
                                   
mc1.execute_macro()

##############################
# Generate mmcif
##############################  
if '--mmcif' in sys.argv:
    import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    
    fname = os.path.join(top_dir,"data/XLs_all_2020.csv")

    # Add publication
    po.system.title = "Implications of a multiscale structure of the yeast Nuclear Pore Complex"
    po.system.citations.append(ihm.Citation(
        pmid='NA',
        title="Implications of a multiscale structure of the yeast Nuclear Pore Complex",
        journal="Molecular Cell", 
        year=2023,
        volume='NA',
        doi='NA',
        page_range='NA',
        authors=['Akey CA', 'Echeverria I', 'Ouch C',
                 'Nudelman I', 'Shi Y', 'Wang J', 'Weiss TM', 'Shi Y',
                 'Chait BT', 'Sali A','Fernandez-Martinez J', 'Rout MP']))
    
    s = po.system
    print("restraint datasets:", [r.dataset for r in s.restraints])
    # Datasets for XL-MS restraint
    for r in s.restraints:
        print('----', r)
        if isinstance(r, ihm.restraint.CrossLinkRestraint):
            r.linker = ihm.cross_linkers.dsso
            print("XL-MS dataset at:", r.dataset.location.path)
            print("Details:", r.dataset.location.details)
          
    # Correct number of output models to account for multiple runs
    protocol = s.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 2500000

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 9999 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=200000,
                            num_models_end=9999))
    
    mg = ihm.model.ModelGroup(name="Cluster 0")
    po.system.state_groups[-1][-1].append(mg)
    e = ihm.model.Ensemble(model_group=mg,
                       num_models=11,
                       post_process=analysis.steps[-1],
                       name="Cluster 0")
    po.system.ensembles.append(e)
    
    # Add the model from RMF - centroid
    rh = RMF.open_rmf_file_read_only('../results/clustering/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [root_hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    m.update()
    del rh
    model = po.add_model(e.model_group)

    # Add ensemble members
    #models = []
    #rmf_ensemble = '../results/cluster0_random_sel.rmf3'
    #rh = RMF.open_rmf_file_read_only(rmf_ensemble)
    #IMP.rmf.link_hierarchies(rh, [hier])
    #for frame_number in range(rh.get_number_of_frames()):        
    #    IMP.rmf.load_frame(rh, RMF.FrameID(frame_number))
    #    mdl.update()
    #    models.append(hier)
    #    model = po.add_model(e.model_group)

    #model_group = ihm.model.ModelGroup(models, name="Cluster 0")
    #state = ihm.model.State([model_group])
    #s.state_groups.append(ihm.model.StateGroup([state]))
    
    #del rh
    
    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'../results/clustering/cluster.0/LPD_{name}.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
    
    # Add uniprot of proteins
    lpep = ihm.LPeptideAlphabet()
    d = 'Segments used for modeling'

    sd_RRAGC = [ihm.reference.SeqDif(181, lpep['D'], lpep['N'], details=d)]
    
    # name : (uniprot id, mutations, [[db_begin, db_end, entity_begin, entity_end]]
    Uniprot={'Lamtor1.0': ('Q6IAA8',[],[49,158,49,158]),
             'Lamtor2.0': ('Q9Y2Q5',[],[1,125, 1,125]),
             'Lamtor3.0': ('Q9UHA4',[],[1,124,1,124]),
             'Lamtor4.0': ('Q0VGL1',[],[1,99,1,99]),
             'Lamtor5.0': ('O43504',[],[1,91,1,91]),
             'RRAGB.0': ('Q5VZM2',[],[1,374,1,374]),
             'RRAGC.0': ('Q9HB90',sd_RRAGC,[62,369,62,369]),
             'SLC38A9.0': ('Q8NBW4',[],[39,96,39,96]),
             'BORCS6.0': ('Q96GS4',[],[199,355,199,355]),
             'ITPA.0': ('Q9BY32',[],[1,194,1,194]),}
    
    for prot, (entry, sd, limits) in Uniprot.items():
        print(prot, entry, limits)
        ref = ihm.reference.UniProtSequence.from_accession(entry)
        ref.alignments.append(ihm.reference.Alignment(
            db_begin=limits[0], db_end=limits[1], entity_begin=limits[2], entity_end=limits[3],seq_dif=sd))
            
        po.asym_units[prot].entity.references.append(ref)

    # Replace local links with DOIs
    repos = []
    #for subdir, zipname in make_archive.ARCHIVES.items():
    #    print('subdir', subdir)
    #    repos.append(ihm.location.Repository(
    #        doi="10.5281/zenodo.3836213", root="../%s" % subdir,
    #        url="https://zenodo.org/record/3836213/files/%s.zip" % zipname,
    #        top_directory=os.path.basename(subdir)))
    
    po.system.update_locations_in_repositories(repos)

    with open('IM_MuSIC_cluster.cif', 'w') as fh:
        ihm.dumper.write(fh, [po.system])
