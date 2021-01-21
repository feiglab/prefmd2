# Example of how the environmental variables that you need to set up to run
# PREFMD2 should look like.
# Modify the paths according to the locations where the various programs are
# installed on your system and add the variables to your own .bashrc file.

# prefMD2.
export PREFMD2_HOME="$HOME/apps/prefmd2"

# MMSTB.
export MMTSBDIR="/apps/mmtsb"
export PATH="$MMTSBDIR/bin:$MMTSBDIR/perl:$PATH"

# CHARMM.
export CHARMMEXEC="/apps/charmm/bin/charmm"
export CHARMMDATA="$MMTSBDIR/data/charmm"

# locPREFMD.
export LOCPREFMD="$HOME/apps/locprefmd"
export MOLPROBITY="/blue/apps/MolProbity"

# Scoring tools.
export RWPLUS_HOME="$HOME/apps/rwplus"
# Assuming that the mdconv and TMscore executables are located in the
# $HOME/bin directory.
export PATH="$HOME/bin:$PATH"

# Scwrl4.
# Assuming that the scwrl4 executable is located in the $HOME/apps/scwrl4
# directory.
export PATH="$HOME/apps/scwrl4:$PATH"
