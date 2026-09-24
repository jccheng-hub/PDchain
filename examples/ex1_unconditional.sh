#!/usr/bin/env bash

pixi run rfd_chain 140-160 \
    --idealize --relax --ca_stdev 1 \
    --model_type protein_mpnn \
    --design_cycles 3 \
    --numdes 5 \
    --outprefix outputs/ex1_rand

# --- Description --- #
cat << 'EOF' > /dev/null
Here, we run `rfd_chain` to randomly generate monomeric proteins with
RFdiffusion, sequence design with ProteinMPNN, and refine the structure with
Rosetta FastRelax.

Run `rfd_chain --help` to display the manual on the command line.

The input argument `140-160` specifies the contigs string. Here we tell
RFdiffusion that we want to generate a backbone 140 to 160 residues long. See
the RFdiffusion github documentation for more information on contigs.

`--idealize` specifies that we want to idealize the backbone after RFdiffusion.
This is done using Rosetta.

`--model_type protein_mpnn` specifies that we want to sequence design with the
ProteinMPNN model.

`--relax --ca_stdev 1` specifies that we want to apply FastRelax with
constraints to the CA atoms in the backbone. A standard deviation of 1 is
applied here. See Rosetta documentation on constraints for more information.

`--design_cycles 2` specifies that we want a total of two rounds of MPNN
sequence design and FastRelax. Here, we have two rounds, so it will proceed as
ProteinMPNN -> FastRelax -> ProteinMPNN -> FastRelax. Iterative rounds of MPNN-
FastRelax typically improves scores from MPNN and Rosetta.

`--numdes 3` specifies that we want a total of three designs.

`--outprefix outputs_ex1/uncond` specifies the path prefix of the output
designs. If the output directory does not exist, it will be made automatically.
EOF
