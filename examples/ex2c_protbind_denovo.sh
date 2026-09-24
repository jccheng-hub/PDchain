#!/usr/bin/env bash

# Generate control
pixi run rfd_chain NONE --skip_mpnn \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 1 \
    --outprefix outputs/ex2_1brs_control

# Generate new designs
pixi run rfd_chain 86-95/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --ppi_hotspots A73 A75 \
    --ppi_mode \
    --fixedres B1-110 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 3 \
    --outprefix outputs/ex2c_1brs_denovo

# --- Description --- #
cat << 'EOF' > /dev/null
Here, we generate de novo protein binders.

The input contigs `80-90/0 A17-138` specifies that we want a binder 80-90
residues long to the target protein. Here, only the residues A17-138 from the
input PDB `inputs/5o45_clean.pdb` is provided to RFdiffusion, and that is the
target RFdiffusion will attempt to generate binders to.

`--ppi_hotspots A37 A39 A98 A100` specifies the hotspot residues for the
binding interface. These were surface-exposed residues on the side of a beta
sheet in the input PDB.

`--ppi_mode` turns on RFdiffusion options presumably helpful for protein binder
design. These were `noise_scale_ca=0.5` and `noise_scale_frame=0.5`. See the
RFdiffusion github documentation for more details.

`--fixedres A17-138` is necessary here to prevent MPNN from redesigning the
residues on the target chain.

The de novo generation of proteins with RFdiffusion can yield unreasonable
backbones. While `rfd_chain` will fill in the sequence with MPNN and refine
with Rosetta, the outputs will not all be high quality. It is recommended to
generate a large pool of designs before subsequent filtering and redesign.
EOF
