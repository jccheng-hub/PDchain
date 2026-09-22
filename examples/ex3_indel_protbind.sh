#!/usr/bin/env bash

# Generate control
rfd_chain NONE --skip_mpnn \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --ca_stdev 1 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 1 \
    --outprefix outputs_ex3/1brs_control

# Generate new designs
rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --ppi_mode \
    --fixedres B1-110 \
    --ss_to_contigs --vary_linkers 1 --ss_trim 2 --inpaint_seq \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 20 \
    --outprefix outputs_ex3/1brs_indel_protbind

# --- Description --- #
cat << 'EOF' > /dev/null
Here, we diversify an existing protein binder to generate more binders similar
to the original.

We first throw the input PDB through the entire design pipeline without any
diffusion or design. The `NONE` keyword at where the contigs is supposed to be
tells `rfd_chain` to skip RFdiffusion. The `--skip_mpnn` option tells
`rfd_chain` to skip MPNN sequence design. This leaves just the Rosetta
refinement with Idealize and FastRelax. This ultimately provides us with a
reference point from which we can compare subsequent designs.

In the second block, we generate new designs by diversifying the existing one
with partial diffusion. The contigs `89-89/0 B1-110` tells RFdiffusion to
keep residues B1-110 while diffusing the 89 remaining residues from the input
PDB, which in this case is all residues of chain A. For partial diffusion to
work, this number has to exactly match the number of remaining residues.

`--ppi_mode` turns on RFdiffusion options presumably helpful for protein binder
design. These were `noise_scale_ca=0.5` and `noise_scale_frame=0.5`. See the
RFdiffusion github documentation for more details.

`--fixedres B1-110` is necessary here to prevent MPNN from sequence designing
the target protein.

`--partial --timesteps 2 --ss_trim 100` applies partial diffusion with 2
timesteps. By default, `rfd_chain` prevents secondary structure from being
rediffused, but this can be adjusted with `--ss_trim [int]`. By specifying a
high value like 100, you ensure that all non-fixed residues will be rediffused
regardless of secondary structure. If you want to keep more secondary structure
from the original binder, you can leave `--ss_trim` with a lower integer (e.g.
`--ss_trim 2`)

`--select_met min:ddg` specifies the selection metric during iterative rounds
of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN
confidence scores will be accepted after refinement, but this flag will make it
so that acceptance/rejection depends solely on the specific metric. Here, the
`min:` prefix specifies that we want lower values of `ddg`. If you want to
maximize some metric value instead, you would use the `max:` prefix (e.g.
`--select_met max:protein_mpnn_score`)
EOF
