# Example of how the environmental variables that you need to set up to run
# the PREFMD2 hybridization stage should look like.
# Modify the paths according to the locations where the various programs are
# installed on your system and add the variables to your own .bashrc file.

# HHsuite.
# HHsuite databases.
export HHSUITE_SEQ_DB="$HOME/data/hhsuite/databases/uc30/uc30"
export HHSUITE_PDB_DB="$HOME/data/hhsuite/databases/pdb70/pdb70"
# Assuming that the HHsuite executables are located in the
# $HOME/apps/hhsuite/bin directory.
export PATH="$HOME/apps/hhsuite/bin:$PATH"

# Rosetta.
export ROSETTA_HOME="$HOME/apps/rosetta/rosetta_src_2020.08.61146_bundle"

# Other executables.
# Assuming that the TMalign and GNU parallel executables are located in the
# $HOME/other/bin directory.
export PATH="$HOME/other/bin:$PATH"
