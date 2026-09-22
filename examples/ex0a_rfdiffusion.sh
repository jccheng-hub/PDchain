#!/usr/bin/env bash

pixi run -e rfdiffusion rfdiffusion \
    'contigmap.contigs=[150-150]' \
    inference.output_prefix=outputs_ex0a/uncond \
    inference.num_designs=3 \
    diffuser.T=15

# --- Description --- #
cat << 'EOF' > /dev/null
The example above runs RFdiffusion directly from within this pixi workspace.

Essentially, the command [rfdiffusion] is linked directly to the
[run_inference.py] script from the original RFdiffusion repository.

The prefix [pixi run -e rfdiffusion] just tells pixi to activate the rfdiffusion
environment before running the command that follows, which in this case is the
command [rfdiffusion 'contigmap.contigs=[150-150]' ...].

For more details on how to run RFdiffusion, visit the official github:
https://github.com/RosettaCommons/RFdiffusion.git
EOF
