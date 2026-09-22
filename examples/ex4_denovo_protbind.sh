#!/usr/bin/env bash

rfd_chain 80-90/0 A17-138 \
    --inpdb inputs/5o45_clean.pdb \
    --ppi_hotspots A37 A39 A98 A100 \
    --ppi_mode \
    --fixedres A17-138 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 20 \
    --outprefix outputs_ex3/5o45_protbind

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
