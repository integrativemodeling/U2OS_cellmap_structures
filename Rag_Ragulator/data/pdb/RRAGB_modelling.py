# Homology modeling by the automodel class
#
# Demonstrates how to refine only a part of the model.(Here RRAGB sequence corresponding to 6wj2_F)
#
# You may want to use the more exhaustive "loop" modeling routines instead.
#
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

# Override the 'select_atoms' routine in the 'automodel' class:
# (To build an all-hydrogen model, derive from allhmodel rather than automodel
# here.)
class MyModel(automodel):

    def special_patches(self,aln):
        self.rename_segments(segment_ids=('A','B','C','D','E','F','G','H'),renumber_residues=[49,1,1,1,1,1,62,39])
    def select_atoms(self):
        
        # chain F (required for multi-chain models) and loop modelling for few residues from chain B, C, D, E:
         return selection(self.residue_range('1:F', '374:F'), self.residue_range('61:B', '64:B'), self.residue_range('124:B','125:B'), self.residue_range('1:C','2:C'),self.residue_range('122:C','124:C'), self.residue_range('1:D','2:D'), self.residue_range('97:D','99:D'), self.residue_range('90:E','91:E'))


env = environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
# selected atoms do not feel the neighborhood
env.edat.nonbonded_sel_atoms = 2

# Be sure to use 'MyModel' rather than 'automodel' here!
a = MyModel(env,
            alnfile  = 'multi_chain1.ali',     # alignment filename
            knowns   = '6wj2',              # codes of the templates
            sequence = 'RRAGB',              # code of the target
            assess_methods=(assess.DOPE,assess.GA341,assess.normalized_dope))

a.starting_model= 1                # index of the first model
a.ending_model  = 10               # index of the last model
                                   # (determines how many models to calculate)
a.make()                           # do homology modeling

