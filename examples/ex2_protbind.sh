# Control
pixi run -w PDchain rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 1 \
    --fixedres B1-110 \
    --select_met min:ddg \
    --natbias 10 \
    --outprefix outputs/ex2_1brs_control

# Native Backbone
pixi run -w PDchain rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres B1-110 \
    --select_met min:ddg \
    --outprefix outputs/ex2a_1brs_natbb

# Partial Diffusion
pixi run -w PDchain rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres B1-110 \
    --select_met min:ddg \
    --model_ckpt Complex_base \
    --partial --timesteps 1 \
    --outprefix outputs/ex2b_1brs_partial

# Indel
pixi run -w PDchain rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres B1-110 \
    --select_met min:ddg \
    --model_ckpt Complex_base \
    --ss_to_contigs --vary_linkers 1 --ss_trim 1-3 \
    --outprefix outputs/ex2c_1brs_indel

# De Novo
pixi run -w PDchain rfd_chain 86-95/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --fixedres B1-110 \
    --select_met min:ddg \
    --model_ckpt Complex_base \
    --ppi_hotspots B73 B75 \
    --outprefix outputs/ex2d_1brs_denovo
