#!/usr/bin/env bash

pixi run -e ligandmpnn ligandmpnn \
    --pdb_path "inputs/1a53_clean.pdb" \
    --out_folder "ex0b_ligandmpnn" \
    --model_type protein_mpnn \
    --batch_size 10 \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1 \
    --temperature 0.1 \
    --repack_everything 0 \
    --ligand_mpnn_use_side_chain_context 1

# --- Description --- #
cat << 'EOF' > /dev/null
The above is an example of running LigandMPNN directly from this pixi workspace.

[pixi run -e ligandmpnn] tells pixi to activate the [ligandmpnn] environment
before running the command that follows, which is [ligandmpnn --pdb_path ...].

The command [ligandmpnn] is linked to the [run.py] script from the LigandMPNN
github page. See the original LigandMPNN github for more details:
https://github.com/dauparas/LigandMPNN.git
EOF
