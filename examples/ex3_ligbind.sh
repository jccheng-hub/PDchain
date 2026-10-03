#!/usr/bin/env bash

# Control
rfd_chain SKIP \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 1 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --model_type soluble_mpnn \
    --skip_mpnn \
    --outprefix outputs/ex3_5rgf_control

# Native backbone
rfd_chain SKIP \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --model_type soluble_mpnn --sc_context \
    --outprefix outputs/ex3a_5rgf_natbb

# Partial diffusion
rfd_chain \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --model_type soluble_mpnn --sc_context \
    --model_ckpt ActiveSite \
    --partial --timesteps 1 \
    --outprefix outputs/ex3b_5rgf_partial

# Indel diffusion
rfd_chain \
    --inpdb inputs/5rgf_clean.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --model_type soluble_mpnn --sc_context \
    --model_ckpt ActiveSite \
    --ss_to_contigs --ss_trim 2-3 --vary_linkers 1 \
    --outprefix outputs/ex3c_5rgf_indel

# Inverse rotamer generation
gen_invrots inputs/5rgf_clean.pdb A50:N3 A127:N3 X1 \
    --numrots 10 \
    --parallel 1 \
    --outdir outputs/ex3d_invrots

# De novo
rfd_chain \
    --indir outputs/ex3d_invrots \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 B127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --model_type soluble_mpnn --sc_context \
    --model_ckpt ActiveSite \
    --fixbbres A49-51 B126-128 --tot_range 180-220 \
    --outprefix outputs/ex3d_5rgf_denovo
