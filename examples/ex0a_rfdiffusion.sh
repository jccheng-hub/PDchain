#!/usr/bin/env bash

pixi run -e rfdiffusion rfdiffusion \
    diffuser.T=15 \
    'contigmap.contigs=[A1-50/10-20/A66-247]' \
    inference.input_pdb="inputs/1a53_clean.pdb" \
    inference.output_prefix="outputs/ex0a_1a53_remodel" \
    inference.num_designs=3 

# --- Description --- #
cat << 'EOF' > /dev/null
The example above runs RFdiffusion directly from within this pixi workspace.

Essentially, the command [rfdiffusion] is linked directly to the
[run_inference.py] script from the original RFdiffusion repository.

The prefix [pixi run -e rfdiffusion] just tells pixi to activate the rfdiffusion
environment before running the command that follows, which in this case is the
command [rfdiffusion diffuser.T=15 ...].

For more details on how to run RFdiffusion, visit the official github:
https://github.com/RosettaCommons/RFdiffusion.git
EOF
