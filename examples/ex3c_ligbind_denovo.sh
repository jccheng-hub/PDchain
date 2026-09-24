#!/usr/bin/env bash

# Generate control
pixi run rfd_chain NONE --skip_mpnn \
    --inpdb inputs/5rgf_clean.pdb \
    --fixedres A50 A127 \
    --ligname 6NT \
    --idealize --relax \
    --ca_stdev 1 \
    --fix_stdev 0.5 --lig_stdev 0.5 --ap_stdev 1.0 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 1 \
    --outprefix outputs/ex3_5rgf_control

# Generate inverse rotamers
[[ $(echo outputs/ex3c_invrots/*.pdb | wc -w) -gt 10 ]] || {
    pixi run gen_invrots inputs/5rgf_clean.pdb A50:N3 A127:N3 X1 \
        --numrots 10 \
        --outdir outputs/ex3c_invrots
}

# Generate designs
pixi run rfd_chain \
    --indir outputs/ex3c_invrots \
    --fixedres A50 B127 \
    --fixbbres A49-51 B126-128 \
    --ligname 6NT \
    --mintot 180 --addtot 40 --mingap 20 --addgap 80 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --fix_stdev 0.5 --lig_stdev 0.5 --ap_stdev 1.0 \
    --design_cycles 5 \
    --select_met min:ddg \
    --numdes 2000 \
    --outprefix outputs/ex3c_5rgf_denovo

# --- Description --- #
cat << 'EOF' > /dev/null
Here, we take an existing ligand-binding protein, complexed with its ligand, and
diversity it to generate similar ligand-binding proteins.

We first apply partial diffusion to noise the starting backbone before
subsequent MPNN-FastRelax. See RFdiffusion github for more details on partial
diffusion.

`--inpdb inputs/1a53_clean.pdb` specifies the input PDB path.

`--fixedres A52 A158` specifies the residues to fix during design. In this case,
the residues 52 and 158 on chain A will be fixed. These were chosen in this
example because they made direct hydrogen-bonding interactions with the ligand.

`--ligname IGP` specifies the three-letter code of the ligand in the input PDB.

`--partial --timesteps 3 --ss_trim 100` applies partial diffusion with 3
timesteps. By default, `rfd_chain` prevents secondary structure from being
rediffused, but this can be adjusted with `--ss_trim [int]`. By specifying a
high value like 100, you ensure that all non-fixed residues will be rediffused
regardless of secondary structure.

`--lig_stdev 1 --fix_stdev 1 --ap_stdev 1` adds a series of constraints during
Rosetta refinement. `--lig_stdev 1` and `--fix_stdev 1` specifies coordinate
constraints with a standard deviation of 1 for ligand atoms and fixed residue
atoms. In this case, the ligand IGP and the residues A52 and A158 will be
constrained to their starting coordinates. `--ap_stdev 1` adds a series of
AtomPair constraints between the ligand atoms and the fixed residue atoms, also
with a standard deviation of 1. See Rosetta documentation on constraints for
more information.

`--select_met min:ddg` specifies the selection metric during iterative rounds
of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN
confidence scores will be accepted after refinement, but this flag will make it
so that acceptance/rejection depends solely on the specific metric. Here, the
`min:` prefix specifies that we want lower values of `ddg`. If you want to
maximize some metric value instead, you would use the `max:` prefix (e.g.
`--select_met max:protein_mpnn_score`)
EOF
