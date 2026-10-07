load ./ex4_5rgf_denovo_0134_gen3_aywpxz_gen0_orides.pdb
load ./ex4_5rgf_denovo_0134_gen3_aywpxz_gen1_xoytye.pdb
load ./ex4_5rgf_denovo_0134_gen3_aywpxz_gen2_vdnthl.pdb
load ./ex4_5rgf_denovo_0134_gen3_aywpxz_gen3_ijzzqn.pdb
load ./ex4_5rgf_denovo_0134_gen3_aywpxz_gen4_blsucl.pdb
load ./ex4_5rgf_denovo_0134_gen3_aywpxz_gen5_dmyczu.pdb
load ex4_5rgf_denovo_0134_gen3_aywpxz_gen6_vbgmll.pdb
join_states lineage, all, -2
disable not lineage
spectrum count, rainbow, lineage and name CA
set movie_fps, 5
set movie_loop, 0
