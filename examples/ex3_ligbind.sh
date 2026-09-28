#!/usr/bin/env bash

# Control
pixi run -w PDchain rfd_chain SKIP \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax \
    --model_type soluble_mpnn \
    --ca_stdev 1 \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --design_cycles 3 \
    --numdes 1 \
    --fixedres A50 A127 --ligname 6NT \
    --select_met min:ddg \
    --natbias 10 \
    --outprefix outputs/ex3_5rgf_control

# Native Backbone
pixi run -w PDchain rfd_chain SKIP \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax \
    --model_type soluble_mpnn \
    --ca_stdev 1 \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres A50 A127 --ligname 6NT \
    --select_met min:ddg \
    --outprefix outputs/ex3a_5rgf_natbb

# Partial Diffusion
pixi run -w PDchain rfd_chain \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax \
    --model_type soluble_mpnn \
    --ca_stdev 1 \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres A50 A127 --ligname 6NT \
    --select_met min:ddg \
    --model_ckpt ActiveSite \
    --partial --timesteps 1 \
    --outprefix outputs/ex3b_5rgf_partial

# Indel
pixi run -w PDchain rfd_chain \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax \
    --model_type soluble_mpnn \
    --ca_stdev 1 \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres A50 A127 --ligname 6NT \
    --select_met min:ddg \
    --model_ckpt ActiveSite \
    --ss_to_contigs --ss_trim 2-3 --vary_linkers 1 \
    --outprefix outputs/ex3c_5rgf_indel

# Generate Inverse Rotamers
pixi run -w PDchain gen_invrots inputs/5rgf_clean.pdb A50:N3 A127:N3 X1 \
    --numrots 10 \
    --parallel 1 \
    --outdir outputs/ex3d_invrots

# De Novo
pixi run -w PDchain rfd_chain \
    --indir outputs/ex3d_invrots \
    --idealize --relax \
    --model_type soluble_mpnn \
    --ca_stdev 1 \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres A50 B127 --ligname 6NT \
    --select_met min:ddg \
    --model_ckpt ActiveSite \
    --fixbbres A49-51 B126-128 \
    --mintot 180 --addtot 40 \
    --outprefix outputs/ex3d_5rgf_denovo
