#!/usr/bin/env bash

# Control
pixi run -w PDchain rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 1 \
    --select_met min:ddg \
    --fixedres B1-110 \
    --skip_mpnn \
    --outprefix outputs/ex2_1brs_control

# Native backbone
pixi run -w PDchain rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --model_type protein_mpnn \
    --fixedres B1-110 \
    --outprefix outputs/ex2a_1brs_natbb

# Partial diffusion
pixi run -w PDchain rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --model_type protein_mpnn \
    --fixedres B1-110 \
    --model_ckpt Complex_base \
    --partial --timesteps 1 \
    --outprefix outputs/ex2b_1brs_partial

# Indel diffusion
pixi run -w PDchain rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --model_type protein_mpnn \
    --fixedres B1-110 \
    --model_ckpt Complex_base \
    --ss_to_contigs --vary_linkers 1 --ss_trim 2-3 \
    --outprefix outputs/ex2c_1brs_indel

# De novo
pixi run -w PDchain rfd_chain 86-95/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --model_type protein_mpnn \
    --fixedres B1-110 \
    --model_ckpt Complex_base \
    --ppi_hotspots B27 B38 B54-59 B82-85 B101-104 \
    --outprefix outputs/ex2d_1brs_denovo
