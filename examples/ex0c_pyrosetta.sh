#!/usr/bin/env bash

pixi run -e pyrosetta python << 'EOF'
from pathlib import Path
import pyrosetta
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

pyrosetta.init('-relax:default_repeats 1')

xml_string = '''
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="r15" weights="ref2015.wts"/>
        <ScoreFunction name="r15_cst" weights="ref2015_cst.wts"/>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
        <True name="FullPose"/>
    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
        <ResfileCommandOperation name="relax_task" command="NATAA" residue_selector="FullPose"/>
    </TASKOPERATIONS>
    <MOVERS>
        <AddConstraints name="add_csts">
            <CoordinateConstraintGenerator name="ca_cst" sd="1" ca_only="1"/>
        </AddConstraints>
        <FastRelax name="relax" task_operations="relax_task" scorefxn="r15_cst"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover_name="add_csts"/>
        <Add mover_name="relax"/>
    </PROTOCOLS>
    <OUTPUT scorefxn="r15"/>
</ROSETTASCRIPTS>
'''
xml = XmlObjects.create_from_string(xml_string).get_mover('ParsedProtocol')

pose = pyrosetta.pose_from_pdb('inputs/1a53_clean.pdb')
xml.apply(pose)

Path('outputs').mkdir(parents=True, exist_ok=True)
pose.dump_pdb('outputs/ex0c_1a53_rlx.pdb')
EOF

# --- Description --- #
cat << 'EOF' > /dev/null
Here, we run a pyrosetta script directly from the pyrosetta environment in this
pixi workspace.

The example above has the pyrosetta script embedded as a heredoc, but writing a
separate pyrosetta script and then executing the file from the pyrosetta
environment will work as well.

The [pixi run -e pyrosetta] prefix activates the pyrosetta environment in this
workspace before running the following command. In this case, we are running a
python script from within this pyrosetta environment.

Here, we take the input pdb at [inputs/1a53_clean.pdb] and apply two movers. The
first mover adds constraints to the CA atoms, and the second mover applies the
FastRelax mover with those constraints. The output is then saved to
[outputs/ex0c_1a53_rlx.pdb].

For more information on how to run Rosetta and Pyrosetta, refer to the Rosetta
and PyRosetta documentation.
EOF
