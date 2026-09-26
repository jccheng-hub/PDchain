#!/usr/bin/env bash

# This file provides examples on how to run RFdiffusion, LigandMPNN, and
# PyRosetta directly.
# Refer to their respective githubs and documentation for more info.

# RFdiffusion example
pixi run -e rfdiffusion rfdiffusion \
    diffuser.T=15 \
    'contigmap.contigs=[A1-50/10-20/A66-247]' \
    inference.input_pdb="inputs/1a53_clean.pdb" \
    inference.output_prefix="outputs/ex0a_1a53_remodel" \
    inference.num_designs=3 

# LigandMPNN example
pixi run -e ligandmpnn ligandmpnn \
    --pdb_path "inputs/1a53_clean.pdb" \
    --out_folder "outputs" \
    --model_type protein_mpnn \
    --batch_size 10 \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1 \
    --temperature 0.1 \
    --repack_everything 0 \
    --ligand_mpnn_use_side_chain_context 1

# PyRosetta example
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
