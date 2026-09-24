#!/usr/bin/env bash

# Generate control
pixi run rfd_chain NONE --skip_mpnn \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --ca_stdev 1 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 1 \
    --outprefix outputs/ex2_1brs_control

# Generate new designs
pixi run rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --ppi_mode --partial --timesteps 1 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 3 \
    --outprefix outputs/ex2a_1brs_partial

# --- Description --- #
cat << 'EOF' > /dev/null
Here, we diversify an existing protein binder with partial diffusion to generate
more binders similar to the original. In partial diffusion, input backbone atoms
are provided as starting points for noising and denoising, essentially allowing
RFdiffusion to tweak the backbone coordinates before subsequent sequence design.

We first throw the input PDB through the entire design pipeline without any
diffusion or design. The `NONE` keyword at where the contigs is supposed to be
tells `rfd_chain` to skip RFdiffusion. The `--skip_mpnn` option tells
`rfd_chain` to skip MPNN sequence design. This leaves just the Rosetta
refinement with Idealize and FastRelax. This ultimately provides us with a
reference point from which we can compare subsequent designs.

In the second block, we generate new designs by diversifying the existing one
with partial diffusion.

`--ppi_mode` turns on RFdiffusion options presumably helpful for protein binder
design. These were `noise_scale_ca=0.5` and `noise_scale_frame=0.5`. See the
RFdiffusion github documentation for more details.

`--partial` turns on partial diffusion. Note that the output diffused structure
will always match the input structure in length with partial diffusion.

`--fixedres B1-110` is necessary here to prevent MPNN from sequence designing
the target protein.

`--select_met min:ddg` specifies the selection metric during iterative rounds
of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN
confidence scores will be accepted after refinement, but this flag will make it
so that acceptance/rejection depends solely on the specific metric. Here, the
`min:` prefix specifies that we want lower values of `ddg`. If you want to
maximize some metric value instead, you would use the `max:` prefix (e.g.
`--select_met max:protein_mpnn_score`)
EOF
