# Control
pixi run -w PDchain rfd_chain NONE \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn --natbias 10 \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 1 \
    --outprefix outputs/ex2_1brs_control

# Partial Diffusion
pixi run -w PDchain rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --partial --timesteps 1 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 3 \
    --outprefix outputs/ex2a_1brs_partial

# Indel
pixi run -w PDchain rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --ss_to_contigs --vary_linkers 1 --ss_trim 1-3 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 3 \
    --outprefix outputs/ex2b_1brs_indel

# De Novo
pixi run -w PDchain rfd_chain 86-95/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --ppi_hotspots A73 A75 \
    --fixedres B1-110 \
    --idealize --relax \
    --model_type protein_mpnn \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 3 \
    --outprefix outputs/ex2c_1brs_denovo
