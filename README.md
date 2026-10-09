# PDchain
Command-line tools for chaining protein design software (RFdiffusion, LigandMPNN, and PyRosetta) using a pixi workspace. Can iteratively optimize designs based on desired metrics with an *in silico* continuous evolution system.

Supported platforms: osx-arm64, linux-64, linux-aarch64.

## Table of Contents
- [Installation](#installation)
- [Running PDchain Commands](#running-pdchain-commands)
- [Unconditional Monomer Generation](#1---unconditional-monomer-generation)
- [Protein Binder Design](#2---protein-binder-design)
- [Ligand Binder / Enzyme Design](#3---ligand-binder--enzyme-design)
- [Continuous Evolution of Designs](#4---continuous-evolution-of-designs)
- [Fold Validation Methods](#fold-validation-methods)
- [Directly Running RFdiffusion, LigandMPNN and PyRosetta](#directly-running-rfdiffusion-ligandmpnn-and-pyrosetta)
- [Built-In Metrics](#built-in-metrics)
## Installation

Run the following command to install the newest version of pixi.
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
Restart your terminal or source your shell's rc file (e.g. `source ~/.bashrc`) to complete the installation. You can check if pixi is registered with `command -V pixi`.

See https://pixi.prefix.dev/latest/installation/ for more information.

> [!IMPORTANT]
This repository's installation script will automatically pull PyRosetta as a dependency. While the original code in this repository is open-source, PyRosetta is not free for commercial use (free for academic, non-profit, and government institutions). Please ensure you are not violating PyRosetta's terms of service by having the appropriate license before running this repository's installation script.

Once you have pixi installed, run the following to download the repository and navigate into it.
```bash
git clone --recurse-submodules https://github.com/jccheng-hub/PDchain.git && cd PDchain
```

Run the following for a CPU-only installation. This installation is more consistent across a variety of systems, but it doesn't leverage GPU acceleration for RFdiffusion.
```bash
cp box/tomls/cpu/pixi.* . && pixi run install && pixi workspace register --name PDchain
```

If you have an Apple Silicon Mac or a CUDA-compatible Linux/WSL and you want to leverage GPU acceleration for RFdiffusion, run the following instead.
```bash
cp box/tomls/gpu/pixi.* . && pixi run install && pixi workspace register --name PDchain
```

The `pixi workspace reigster --name PDchain` command will register `PDchain` as the workspace name. This allows you to use the `-w PDchain` option for various pixi commands. More on this in the [Running PDchain Commands](#running-pdchain-commands) section.

> [!NOTE]
The CUDA version listed in the pixi.toml file for the Linux GPU installation of RFdiffusion is 11.8. Some newer GPUs from NVIDIA are incompatible with CUDA 11.8 and requires newer versions of CUDA. If this is the case, then the installation may fail despite having a CUDA-compatible GPU. I'm not entirely sure how to get around this (though I would be surprised if there wasn't a way around this), but you can attempt to fix the pixi.toml file to run on a newer version of CUDA. This may require additional adjustments like changing the python/pytorch/dgl versions, shifting to installation from pypi instead of conda-forge, changing the available channels, etc.

Installation will take some time (~10-20 minutes) as it will download all of the weights for RFdiffusion and LigandMPNN in addition to PyRosetta.

## Running PDchain Commands

You can enter a *subshell* with PDchain's default environment with *any* of the following lines:
```bash
pixi shell                                  # Must be in PDchain subdirectory
pixi shell -w PDchain                       # Works anywhere if PDchain is a registered workspace name
pixi shell -m /path/to/PDchain/pixi.toml    # Works anywhere
```
You can always exit out of your PDchain subshell with the command `exit`.

If you want to incorporate PDchain's default environment in your *current shell*, run the following instead:
```bash
source pdchain.sh    # pdchain.sh should have spawned in PDchain's root directory upon installation
```

Generally, booting up a subshell is better for interactive use, and incorporating into current shell is better for scripting use.

Once inside the PDchain default environment, you can now directly access a variety of command-line tools in the PDchain default environment. For example:
```bash
pdchain-info            # Report available tools on PDchain
rfd_chain --help        # Help documentation for rfd_chain
esmfold_relax --help    # Help documentation for esmfold_relax
```

If you want to run PDchain commands without having to change/switch environments whatsoever, you can run *any* of the following lines:
```bash
pixi run rfd_chain --help                                  # Only works in PDchain subdirectory
pixi run -w PDchain rfd_chain --help                       # Works anywhere if PDchain is a registered workspace name
pixi run -m /path/to/PDchain/pixi.toml rfd_chain --help    # Works anywhere
```

> [!TIP]
If you don't mind using the free and open source version of PyMOL, it is available in PDchain's default environment. You can simply run `pymol` (if you're in the environment) or `pixi run -w PDchain pymol` (if you're outside the environment).

## 1 - Unconditional Monomer Generation
> [!IMPORTANT]
All examples in this README are meant to be executed in PDchain's default environment from the `examples` directory. In other words, run `source pdchain.sh` (or `pixi shell -w PDchain`) and `cd examples` before attempting to run these example blocks of code.

```bash
rfd_chain 140-160 \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --model_type protein_mpnn \
    --outprefix outputs/ex1_rand
```

Here, we run `rfd_chain` to randomly generate monomeric proteins with RFdiffusion, sequence design with ProteinMPNN, and refine the structure with Rosetta FastRelax.

Run `rfd_chain --help` to display the manual on the command line.

The input argument `140-160` specifies the contigs string. Here we tell RFdiffusion that we want to generate a backbone 140 to 160 residues long. See the RFdiffusion github documentation for more information on contigs.

`--idealize` specifies that we want to idealize the backbone after RFdiffusion. This is done using Rosetta.

`--model_type protein_mpnn` specifies that we want to sequence design with the ProteinMPNN model.

`--relax --ca_stdev 1` specifies that we want to apply FastRelax with constraints to the CA atoms in the backbone. A standard deviation of 1 is applied here. See Rosetta documentation on constraints for more information.

`--design_cycles 3` specifies that we want a total of three rounds of MPNN sequence design and FastRelax. The output of each round is only accepted if it improves Rosetta score and/or MPNN score.

`--numdes 3` specifies that we want a total of three designs.

`--outprefix outputs/ex1_rand` specifies the path prefix of the output designs.

## 2 - Protein Binder Design

The following command Rosetta refines the input protein-protein complex, which in this case is the barnase-barstar complex, without touching the backbone or sequence. This generates a control (no design done) that we can use as a reference for the actual design runs later in this section.

```bash
rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 1 \
    --select_met min:ddg \
    --fixedres B1-110 \
    --skip_mpnn \
    --outprefix outputs/ex2_1brs_control
```

The `SKIP` keyword at where the contigs is supposed to be tells `rfd_chain` to skip RFdiffusion, and `--skip_mpnn` skips the MPNN sequence design step. The result is a computational control refined by Rosetta.

### 2a - Protein Binder Redesign with Native Backbone

The following command will diversify the binder sequence. RFdiffusion is not applied here (just MPNN-FastRelax).

```bash
rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres B1-110 \
    --model_type protein_mpnn \
    --outprefix outputs/ex2a_1brs_natbb
```

`--fixedres B1-110` prevents MPNN from sequence designing the target protein (the barnase on chain B).

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN confidence scores will be accepted after refinement, but this flag will make it so that acceptance/rejection depends solely on the specified metric. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section on [Built-In Metrics](#built-in-metrics) for available metrics.

### 2b - Protein Binder Redesign with Partial Diffusion
The following command will noise/denoise the backbone of the binder (the barstar on chain A) with partial diffusion before applying cycles of MPNN-FastRelax.

```bash
rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres B1-110 \
    --model_type protein_mpnn \
    --model_ckpt Complex_base \
    --partial --timesteps 3 \
    --outprefix outputs/ex2b_1brs_partial
```

`--model_ckpt Complex_base` specifies the model checkpoint used for RFdiffusion. `Complex_base` is recommended for protein binder design.

`--partial` turns on partial diffusion. Note that the output diffused structure will always match the input structure in length with partial diffusion.

`--timesteps 3` specifies the number of timesteps for RFdiffusion. When using partial diffusion, this number is allowed to dip below 15. The higher number of timesteps, the greater the deviation from the original input structure.

### 2c - Protein Binder Redesign with Indels

The following command will diversify the binder by rediffusing loop regions while allowing for insertions and deletions.

```bash
rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres B1-110 \
    --model_type protein_mpnn \
    --model_ckpt Complex_base \
    --ss_to_contigs --ss_trim 1 --vary_linkers 1 \
    --outprefix outputs/ex2c_1brs_indel
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. Loop residues will be masked from RFdiffusion.

`--ss_trim 1` will allow helix/sheet residues neighboring loop residues to be treated as loop residues (which are then masked from RFdiffusion). Sheets are trimmed by the specified value, and helices are trimmed by twice the specified value.

`--vary_linkers 1` will allow the loop regions to vary in length by 1 residue. So if a loop region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

### 2d - Protein Binder De Novo Design

The following command will generate de novo protein binders.

```bash
rfd_chain 86-95/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres B1-110 \
    --model_type protein_mpnn \
    --model_ckpt Complex_base \
    --ppi_hotspots B27 B38 B54-59 B82-85 B101-104 \
    --outprefix outputs/ex2d_1brs_denovo
```

Here, we provided the contigs `86-95/0 B1-110` to specify that we want to generate a backbone 86-95 residues long while preserving our target (chain B residues 1-110).

`--ppi_hotspots B27 B38 B54-59 B82-85 B101-104` provides hotspots for RFdiffusion, and it will attempt to generate a backbone near those specified residues.

De novo design often requires generating thousands of designs and computationally screening through them to get something "reasonable". If we compare the designs to the original control, we are likely to see designs that actually perform worse on many desirable metrics. For example, if we were to check the ddg of the outputs compared to control, (e.g. with `grep -H '^ddg ' outputs/ex2*.pdb`), we are likely to see the control outperform most if not all of the de novo designs. While barnase-barstar is an incredibly tight complex (so the bar [HA!] is pretty high here), large scale generation and screening will still often be necessary to have high confidence in your designs.

Because we are unlikely to get excellent designs right away, it is often necessary to take an agreeable-but-not-excellent de novo design and diversify around that with indels and partial diffusion to get something better. This is the main reason why the *in silico* continuous evolution system in PDchain was built. See the section on [Contnuous Evolution of Designs](#4---continuous-evolution-of-designs) for more details.

## 3 - Ligand Binder / Enzyme Design
The following example takes a protein complexed with a ligand (in this case, a Kemp eliminase complexed with its transition-state analog) and refines it with Rosetta (without any sequence or backbone design). This command is meant to generate a computational control to serve as a reference point for later design runs.

```bash
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
```

Setting the contigs to `SKIP` skips RFdiffusion entirely, and `--skip_mpnn` skips the MPNN sequence design step.


### 3a - Ligand Binder / Enzyme Redesign with Native Backbone

The following command will redesign the sequence using the native backbone (no RFdiffusion done here).

```bash
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
```

`--fixedres A50 A127` fixes the the residues 50 and 127 on chain A. In this example, they correspond to an Asp and Gln that hydrogen-bond to the ligand 6NT.

`--ligname 6NT` specifies the ligand from the input structure you want to keep. Without specifying the ligand names, `rfd_chain` will ignore it entirely.

`--lig_stdev 0.5` applies Rosetta coordinate constraints to the ligand. This one was done with a standard deviation of 0.5, which makes it stronger than the CA coordinate constraints specified with `--ca_stdev 1` in this example.

`--ap_stdev 0.5` applies Rosetta distance (AtomPair) constraints between the ligand and fixed residues. This option only works when both `--fixedres` and `--ligname` are specified. These constraints will preserve the relative geometry between the ligand and the fixed residues (the Asp on A50 and the Gln on A127 in this case).

`--model_type soluble_mpnn` specifies that the use of the SolubleMPNN model, which was not trained to be ligand-aware. When you specify a `--model_type` that isn't `ligand_mpnn` but you include a ligand, the model you specified will be applied first, then the residues within 8 angstroms of the ligand will be redesigned with LigandMPNN. This ensures that at least the residues surrounding the ligand are being redesigned in a ligand-conscious manner. If you want to override LigandMPNN's design with Rosetta's FastDesign for those residues within 8 angstroms, add the flag `--rosetta_lig_nbr`.

`--sc_context` allows LigandMPNN to use fixed residues as additional ligand atoms. Presumably this makes the sequence design of LigandMPNN sensitive to the existing rotamer.

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN confidence scores will be accepted after refinement, but this flag will make it so that acceptance/rejection depends solely on the specified metric. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section on [Built-in Metrics](#built-in-metrics) for more metrics.


### 3b - Ligand Binder / Enzyme Redesign with Partial Diffusion

The following example will noise/denoise the input backbone with partial diffusion before MPNN-FastRelax.

```bash
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
    --partial --timesteps 3 \
    --outprefix outputs/ex3b_5rgf_partial
```

`--partial` enables partial diffusion, and `--timesteps 3` sets the number of noising/denoising steps to 3. Increase the value for `--timesteps` if you want more diversity from the original input.

### 3c - Ligand Binder / Enzyme Redesign with Indels

The following example will preserve helix and sheet residues while allowing loop regions to rediffuse with varying lengths, allowing for insertions and deletions.

```bash
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
    --ss_to_contigs --ss_trim 1 --vary_linkers 1 \
    --outprefix outputs/ex3c_5rgf_indel
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. Loop residues will be masked from RFdiffusion.

`--ss_trim 1` will allow helix/sheet residues neighboring loop residues to be treated as loop residues (which are masked from RFdiffusion). Sheets are trimmed by the specified value, and helices are trimmed by twice the specified value.

`--vary_linkers 1` will allow the loop regions to vary in length by 1 residue. So if a loop region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

### 3d - Ligand Binder / Enzyme De Novo Design

> [!NOTE]
This repository integrates the original RFdiffusion, not RFdiffusion2 or RFdiffusion3. While RFdiffusion2 and RFdiffusion3 supports ligand binder / enzyme design without the need to specify starting backbone coordinates or residue indices, the original RFdiffusion does not. In other words, using RFdiffusion for this design task is not the most elegant strategy given the existence of its successors. Still, it is entirely possible to scaffold active sites with the original RFdiffusion, and this example will show how it can be done in the context of this pixi workspace.

The first thing we need is a set of inverse rotamers. The following example shows how we can generate them.

```bash
gen_invrots inputs/5rgf_clean.pdb A50:N3 A127:N3 X1 \
    --numrots 10 \
    --parallel 1 \
    --outdir outputs/ex3d_invrots
```

The command `gen_invrots` will randomly generate inverse rotamers from the input PDB (given as the first argument). The second argument and onward will tell `gen_invrots` exactly what residues to keep. In this example, we are keeping residues A50 (catalytic Asp), A127 (catalytic Gln), and X1 (ligand).

If a residue specified by `gen_invrots` is an amino acid, you have the option of placing those residues on secondary structures. This follows the format of `[residue]:[H/E/N][int]`, where H, E, and N represent helix, sheet, and any (either helix or sheet), and the integer that follows specifies the number of residues to add at each end of that residue. So for example, our `A50:N3` will scaffold residue A50 on either a 7-residue beta-strand or a 7-residue helix with the residue of interest (Asp on A50) being right in the middle.

`--numrots 10` specify the number of inverse rotamers you want to generate. Once the number of generated rotamers exceed this value, `gen_invrots` will stop generating more.

`--parallel 1` specifies the number of parallel jobs you want to run simultaneously to generate inverse rotamers. Here, we only spawn one job at a time, but increasing this value will make generating inverse rotamers faster. This should not exceed the number of CPU cores you have on your system.

`--outdir outputs/ex3d_invrots` specify where we want to store the inverse rotamers. This directory will serve as the input for the actual design run.

Once you have generated your inverse rotamers, you can take a peak at them in pymol to see if they are reasonable with `pymol $(echo outputs/ex3d_invrots/*.pdb | cut -d' ' -f1-10)`. Notice that the chain letter code is now different. While A50 is still on A50, A127 is now on B127, and X1 (the ligand) is now on C1. This is because `gen_invrots` has to split each input residue onto its own chain to prevent potential overlap in both residue number and residue chain when the additional residues for secondary structures are added. As a result, we need to be make sure we specify the new chain+residue ids in the subsequent design command.

The following command will generate a de novo binder to the target ligand using the inverse rotamers we previously generated.

```bash
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
```

`--indir outputs/ex3d_invrots` specifies that we want the input to be a directory, and that input directory is where our previously generated inverse rotamers are located `outputs/ex3d_invrots`. When you specify an input directory (with `--indir`) as opposed to an input pdb (with `--inpdb`), a random .pdb file from that input directory is chosen as the input pdb.

As mentioned previously, what used to be A127 is now B127 because of `gen_invrots` placing each input fixed residue on its own chain.

`--fixbbres A49-51 B126-128` specifies what residues we want fixed during RFdiffusion. Given that our main residues are on A50 and B127, this example here retains an extra residue on both ends for both A50 and B127. Note that we could have gone up to the range of A47-53 and B124-130 because of how we generated those inverse rotamers (7-residue peptide scaffolding each key residue).

`--tot_range 180-220` provides the overall length of the desired monomer. Here, we are requesting a monomer 180-220 residues long. This option in conjunction with `--fixbbres` will trigger `rfd_chain` to run `gen_contigs` to generate a contigs string based on these constraints. You still have the option to specify your own contigs as the first argument if you want more specific control.

Note that the order in which the residues are specified for `--fixbbres` matters. The command `gen_contigs`, which `rfd_chain` automatically accesses when a contigs string is not provided, will add residues in between the fixed regions *in the specified order* as well as additional residues at the N and C termini. In this case, the output will follow the order of `[diffused]-[A49-51]-[diffused]-[B126-128]-[diffused]` in which `[diffused]` represent the backbone regions filled in by RFdiffusion. Had we specified `--fixbbres B126-128 A49-51` instead, the order between the two fixed regions would have been flipped (i.e. `[diffused]-[B126-128]-[diffused]-[A49-51]-[diffused]`).

`--model_ckpt ActiveSite` tells RFdiffusion to use the ActiveSite checkpoint, which is tailored for scaffolding small motifs. This is more necessary here as we are essentially scaffolding small peptides floating in space. This isn't as necessary for the partial diffusion or indel examples above, but it was kept in those examples too for consistency.

> [!NOTE]
RFdiffusion generates backbones around ligands with an auxillary potential (see RFdiffusion documentation for more details). However, even with this potential, RFdiffusion can produce backbones that clash with the input ligand. The command `rfd_chain` was therefore programmed to exclude backbones with excess clashes to ligand atoms. As a result, you may not see the same amount of design outputs as intended. For example, we have `--numdes 3` in the command above, but you may end up with just `ex3d_5grf_denovo_0002.pdb` in the outputs without the first and third designs because they failed the clash checker. If you want to guarantee the exact number of outputs as specified by `--numdes`, add the option `--persistent` to make RFdiffusion try again if it fails. If you want to turn off the clash checker entirely, you can set the threshold to an excessively high number like `--clash_cut 9999` so that every backbone RFdiffusion generates will always pass regardless of ligand clashes. This is not recommended, however, as MPNN may get an unreasonable input and Rosetta may have a hard time resolving those clashes.

As mentioned in the example with de novo protein binder design, it is unlikely to get "reasonable" de novo designs from a small-scale computational run. Getting promising design candidates often require generating thousands of designs and screening through them. Sometimes, even the best designs from a large-scale batch might not meet all of your desired criteria. For example, they might be globular with ample secondary structure composition, but the binding pocket is completely buried/exposed. In these situations, it might be better to redesign these candidates to optimize for those desired properties rather than to repeatedly fish for new designs that meet all of your criteria all at once. See the section on [Continuous Evolution of Designs](#4---continuous-evolution-of-designs) to see how one can evolve designs based on desirable metrics.

## 4 - Continuous Evolution of Designs
Not every design output will possess properties you are looking for; this is especially true for de novo designs. Often times, you end up with designs that check some boxes but not others. In these situations, it may be worth attempting optimization of these designs by using them as starting points for cycles of diversification and selection.

We did some version of this in the protein binding and ligand binding examples in which we used the `--select_met min:ddg` option to tune the MPNN-FastRelax protocol built into `rfd_chain`. As MPNN-FastRelax is applied through multiple cycles, it will only accept the new output if it improves the specified score. In this case, we instructed `rfd_chain` to select for lower values of `ddg`. This means that, at the end of each cycle of MPNN-FastRelax, the output is only accepted if it has lower `ddg` than the input.

The main benefit to `--select_met` is that it can direct MPNN-FastRelax to select for sequences (from its current backbone) more tailored to your design goal, but it is only concerned about a singular metric. In addition, `rfd_chain` itself has no knowledge of how any of your other designs compare to the current one you are optimizing.

The command `evo_rfd_chain` was written to address these limitations. It will take an input batch of designs, generate a pool of designs by diversifying those inputs, then select the "best" in that pool based on specified metrics for further diversification. The output design will then get added into the same pool of designs before the next round of selection and diversification begins.

The remainder of this section will walk through a complete example of ligand binder / enzyme design and subsequent design evolution. Because every command after the first will require the outputs of the previous command as inputs, this repository also provides all outputs from my own run of these commands. These outputs are stored in the `inputs` directory, and they also act as inputs for the subsequent commands. If you want to emulate the entire process (which may take a couple overnight computation runs) with *only* your own outputs, you'll want to adjust the input paths of these code blocks accordingly.

### 4a - Generating Initial Design Pool
The following command will generate a set of inverse rotamers (identical to the ligand binder / enzyme design example).
```bash
gen_invrots inputs/5rgf_clean.pdb A50:N3 A127:N3 X1 \
    --numrots 200 \
    --parallel 1 \
    --outdir outputs/4a_invrots
```
A copy of my own outputs is located at `inputs/4a_invrots`.

The following command will take an input directory of inverse rotamers, then generate an initial design pool.
```bash
rfd_chain \
    --indir inputs/4a_invrots \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --numdes 200 \
    --select_met min:ddg \
    --fixedres A50 B127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --model_type soluble_mpnn --sc_context \
    --model_ckpt ActiveSite \
    --fixbbres A49-51 B126-128 --tot_range 160-200 \
    --persistent \
    --outprefix outputs/4a_oripool/ex4_5rgf_denovo
```
A copy of my own outputs is located at `inputs/4a_oripool`.

### 4b - Selecting Candidates for Evolution
The directory `inputs/4a_oripool` contains an initial batch of de novo designs, and many of these will not have desirable properties. For example, the ligand might be too exposed, the scaffold may be too elongated or loopy, etc.

To avoid having to parse through a large batch of designs by hand, you can leverage the `filter_pool` command to select designs based on your desired criteria.

```bash
filter_pool min:fc_compact min:fc_dev min:rog_ala \
    --indir inputs/4a_oripool \
    --perc 5 \
    --outdir outputs/4b_filtered
```

Open up the outputs in pymol for a sanity check.
```bash
pymol outputs/4b_filtered/*.pdb
```

Here, we are instructing `filter_pool` to rank designs in `inputs/4a_oripool` based on three metrics: `fc_compact`, `fc_dev`, and `rog_ala`. In this case, we are aiming for lower values of each with the `min` prefix. The `max` and `val` prefix are also available (run `filter_pool --help` for more details). Lower values of `fc_compact` and `fc_dev` correlate with pocket formation around the ligand, and lower values of `rog_ala` correlate with the compactness and globularity of the overall protein scaffold. See [Built-In Metrics](#built-in-metrics) for more info on metrics.

The choice of metrics for `filter_pool` may not be immediately obvious *a priori*. You might need to experiment with various metrics to see which combination selects for more reasonable scaffolds depending on your active site. Generally, the three used up here is a good start for ligand binders.

Keep in mind that selecting for more metrics doesn't necessarily correlate with better outputs, especially for this small batch of 200 designs. The algorithm works by finding designs that are *least offensive* across all listed metrics, and trying to be inoffensive for a large number of metrics all at once will yield designs that do inadequately on everything.

`--perc 5` specifies that we want to output the top 5% of designs, which will be stored in `outputs/4b_filtered` in this case.

> [!Note]
The command `filter_pool` ranks designs with a relativistic reverse-elimination approach. Essentially, every specified metric given to `filter_pool` will generate a list ordering the designs based on that metric. In this case, we will have ordered lists for for `fc_compact`, `fc_dev`, and `rog_ala`. The algorithm then starts eliminating the worst designs on each metric. During each elimination round, (at most) three designs are eliminated: the design with the worst `fc_compact`, the design with the worst `fc_dev`, and the design with the worst `rog_ala`. This process repeats itself until all designs are eliminated, and the one eliminated last is considered the best. This ranking approach will therefore highly rank designs that are well-rounded based on the specified metric relative to other designs in the pool.

### 4c - Evolving Designs (Major Rediffusion)
The scaffold `inputs/4b_filtered/ex4_5rgf_denovo_0134.pdb` seems promising. The ligand is a bit too buried, but the scaffold wraps around it nicely without being unreasonably loopy. There are a couple others in this batch that might also work, but I will stick with this one for the following section.

Every design that is outputted by `rfd_chain` will contain information that allows them to be evolved by `evo_rfd_chain`. Notably, information such as the original contigs, the fixed residues, the ligand names, the Rosetta constraints, and various selection metrics are all included in the footer of each PDB (try `cat inputs/4b_filtered/ex4_5rgf_denovo_0134.pdb` to display the contents of our selected design to see for yourself).

We can now run our evolution system on our candidate design. Here, we specify metrics that we might be interested in selecting for in the context of this newer scaffold.
```bash
evo_rfd_chain \
    min:fixedres_perc_loop min:rog_ala val:dsasa=0.85 \
    --inpdb inputs/4b_filtered/ex4_5rgf_denovo_0134.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --model_type soluble_mpnn --rosetta_lig_nbr \
    --model_ckpt ActiveSite \
    --select_met min:ddg \
    --poolsize 200 --poolperc 10 \
    --tlim 8h --stop_at_capacity \
    --ss_to_contigs --ss_trim 2 --vary_linkers 1 \
    --helix_cap 20 --reset_perc 50 \
    --outdir outputs/4c_evopool
```

Here, we are telling `evo_rfd_chain` to select designs within `inputs/selepool` based on three metrics: `fixedres_perc_loop`, `rog_ala`, and `dsasa`. We want to minimize `fixedres_perc_loop` and `rog_ala` while keeping `dsasa` close to `0.85`. `fixedres_perc_loop` measures the percent loop composition around fixed residues, `rog_ala` measures the (approximate) radius of gyration of the protein scaffold, and `dsasa` gives us the percent burial of our ligand. See [Built-In Metrics](#built-in-metrics) for details on metrics.

`--inpdb` can specify input designs directly for `evo_rfd_chain`. Alternatively, you can use `--indir` to specify an input directory if you want to evolve all designs from that input directory.

`--poolsize 200` specifies the maximum number of designs in the output directory, where the design evolution takes place. All designs from the input directory will be copied into the output directory before design evolution begins. If the number of designs in the output directory exceed the poolsize, the worst ones based on the specified metrics will be archived.

`--poolperc 10` specifies the percent of your current poolsize that will be randomly selected as candidates for diversification. Once candidates are randomly selected, the best one (according to your specified metrics) on that list will be selected for diversification. For example, if your current pool has 50 designs, `--poolperc 10` means that 5 designs (10% of 50) from your current pool will be randomly selected as candidates for redesign; the top-ranking design of that list of candidates is then selected for redesign. Tune this number higher if you want a stronger bias for top-ranking designs. A value of 100 means that the top design has a 100% chance of being selected. A value of 0 means that the input design is selected randomly from the pool in an unbiased fashion (the candidate list has to have at least one design at minimum).

`--tlim 8h` specifies the time limit of your run. This example has this command running for 8 hours maximum.

`--stop_at_capacity` tells `evo_rfd_chain` to stop running when the pool size reaches maximum capacity (which is 200 in this example due to `--poolsize 200`). Without this flag, `evo_rfd_chain` will continue until the time limit is hit, and the worst-ranking designs will be archived as needed to keep the poolsize equal or under the limit.

Because of `--ss_to_contigs`, `--ss_trim 2`, and `--vary_linkers 1`, we are doing ss-based rediffusion and varying the loop linker lengths to allow for insertion and deletion. However, this on its own is more suited for more minor backbone optimization rather than more drastic backbone resampling.

To promote more drastic scaffold remodeling with `evo_rfd_chain`, we added the options `--helix_cap 20` and `--reset_perc 50`. With `--helix_cap 20`, we force helices longer than 20 to be rediffused by RFdiffusion. In addition, `--reset_perc 50` specifies the percent chance in which all residues, except for the fixed residues and the stretch of ss residues attached to them, will be rediffused by RFdiffusion.

As `evo_rfd_chain` progresses, the files `ranked.txt` and `ranked_metrics.txt` will spawn in the output directory. The `ranked.txt` file provides an ordered ranked list (from best to worst) of the current pool. The `ranked_metrics.txt` file reports the metrics being used for design evolution.

A copy of the resulting pool from my own design evolution run is included at `inputs/4c_evopool`. Note that your own results may vary from my own due to the randomness of diversification.

### 4d - Evolving Designs (Minor Rediffusion)
Once you have an agreeable scaffold, you can continue sampling and optimizing in the context of that scaffold with lighter diversification strategy. The results of my own run (stored at `inputs/4c_evopool`) ranked `ex4_5rgf_denovo_0134_gen3_aywpxz.pdb` in the top 10. It had a reasonably deep pocket without being excessively buried, and we will be using this in the following example.

Note that our evolved scaffold has differed quite a bit from the original.
```bash
pymol inputs/4b_filtered/ex4_5rgf_denovo_0134.pdb inputs/4c_evopool/ex4_5rgf_denovo_0134_gen3_aywpxz.pdb
```

Because we don't actually want to introduce any major changes to our new scaffold, we can tone down our diversification approach. For example, we can restrict our approach to just partial diffusion + MPNN-FastRelax.
```bash
evo_rfd_chain \
    min:ddg min:fixedres_reu \
    min:fin_reu_per_res min:fixedres_perc_loop \
    --inpdb inputs/4c_evopool/ex4_5rgf_denovo_0134_gen3_aywpxz.pdb \
    --idealize --relax --ca_stdev 1 \
    --design_cycles 3 \
    --model_type soluble_mpnn --rosetta_lig_nbr \
    --model_ckpt ActiveSite \
    --select_met min:ddg \
    --poolsize 200 --poolperc 10 \
    --tlim 8h --stop_at_capacity \
    --partial --timesteps 3 \
    --outdir outputs/4d_evopool
```

Not only did we change our diversification approach (with `--partial --timesteps 3` instead of `--ss_to_contigs --ss_trim 2 ...`), we also changed the metrics that we are selecting for. Back when we were using a more drastic diversification strategy that involved rediffusing entire chunks of the protein, our selection strategy was geared towards scaffold quality (e.g. `fixedres_perc_loop`, `rog_ala`) and suitable pocket formation (e.g. `dsasa`). Because our current diversification strategy is far less drastic on the backbone, we can afford to focus more on sequence-specific metrics such as `ddg` and `fixedres_reu`.

A copy of the resulting pool for this partial diffusion evolution run is located at `inputs/4d_evopool`. As usual, the designs are ordered based on their relative rankings in the `ranked.txt` file (better designs at the top of the list). Copies of the top 3 are stored at `inputs/4d_evopool_top`.

## Fold Validation Methods
### ESMFold and OmegaFold
PDchain comes with a few different commands that can facilitate fold validation. The commands `esmfold_relax` and `omegafold_relax` will leverage ESMFold (server) and OmegaFold (local). Note that these two commands will only work with monomeric proteins.

In the following example, we will take the 3 designs in `inputs/4d_evopool_top` and submit them to the ESMFold webserver, download the outputs, align to our original design with PyMOL, then refine with Rosetta.
```bash
esmfold_relax inputs/4d_evopool_top/*.pdb --outdir outputs/esmfold
```

When using `esmfold_relax`, the raw outputs from ESMFold will be stored in the directory specified by `--foldrepo`. By default, `--foldrepo` points to the PDchain root directory (e.g. `--foldrepo /Users/johncheng/Workspaces/PDchain/foldrepo`). Should you try to fold the same sequence again, `esmfold_relax` will pull the raw structure from `--foldrepo` instead of queuing the ESMFold webserver. Only up to 1000 structures are stored in `--foldrepo`, and the oldest ones will be deleted when this limit is exceeded.

If the ESMFold server didn't time out on you, then you can try viewing the outputs in PyMOL.
```bash
pymol inputs/4d_evopool_top/*.pdb outputs/esmfold/*.pdb
```

Not all designed sequences will have their folds validated. Their CA-RMSDs are stored in the ESMFold PDBs, which you can view with `grep`.
```bash
grep -H '^rmsd_ca ' outputs/esmfold/*.pdb
```

The command `omegafold_relax` behaves identically to `esmfold_relax` except it uses PDchain's local installation of OmegaFold instead.

### AlphaFold3 Server
AlphaFold3 does not have a command-line interface (that I am aware of) for automated submission of jobs. However, PDchain contain commands that can facilitate the submission of jobs through the AF3 webserver. This is done with the `af3webtools` command.
```bash
af3webtools prep inputs/4d_evopool_top/*.pdb --outdir outputs
```

This should spawn a json file in the output directory. This json can then be uploaded to [alphafoldserver.com](#alphafoldserver.com) as job drafts, which can then be subsequently submitted so long as you have enough jobs left for the day.

Should you download a batch of AF3-predicted structures from the AF3 server, they would be downloaded as a zip file. The command `af3webtools` will also facilitate the extraction of these models and the structural alignment to your designs.
```bash
af3webtools unzip inputs/folds_2026_10_09_01_43.zip --refdir inputs/4d_evopool_top --outdir outputs
```

`--refdir` specifies the reference directory that contains the designs. For this to work, the name of the designs has to match what was submitted onto the AF3 server. If `af3webtools prep` was used to prepare the AF3 server json files from the same designs, then the names should already be consistent.

This should spawn three directories: `af3_outputs`, `aln_pdbs`, and `aln_pses`. The `af3_outputs` directory stores the raw data from the zip file. The `aln_pdbs` directory stores the PDB files aligned to your designs. The `aln_pses` directory stores the pymol session files that each contain a design with all 5 aligned AF3 models.

## Directly Running RFdiffusion, LigandMPNN, and PyRosetta
While PDchain provides wrapper commands for stringing together RFdiffusion, LigandMPNN, and PyRosetta, you can also directly access these programs individually.

For RFdiffusion and LigandMPNN, the only difference between how you would run them "normally" and how you would run them in PDchain is the script/command you call and how to call it. The main python script of each program is linked to a command that exists only in their respective environments.

For example, the main script from RFdiffusion is `run_inference.py`; in PDchain, this script is linked to the `rfdiffusion` command, which is only executable from the `rfdiffusion` environment.

In other words, rather than running...
```bash
python /path/to/RFdiffusion/scripts/run_inference.py ...
```
... you would instead run ...
```bash
pixi run -w PDchain -e rfdiffusion rfdiffusion ...
```

`pixi run -w PDchain` calls the PDchain workspace, and `-e rfdiffusion` calls the `rfdiffusion` *environment*. The last `rfdiffusion` is the actual *command* linked to `run_inference.py` from the original RFdiffusion repository.

Here is a complete example for RFdiffusion meant to be run from the `examples` directory:
```bash
pixi run -w PDchain -e rfdiffusion rfdiffusion \
    diffuser.T=15 \
    'contigmap.contigs=[A1-50/10-20/A66-247]' \
    inference.input_pdb="inputs/1a53_clean.pdb" \
    inference.output_prefix="outputs/ex0a_1a53_remodel" \
    inference.num_designs=1
```

Here is an example for LigandMPNN:
```bash
pixi run -w PDchain -e ligandmpnn ligandmpnn \
    --pdb_path "inputs/1a53_clean.pdb" \
    --out_folder "outputs" \
    --model_type protein_mpnn \
    --batch_size 1 \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1 \
    --temperature 0.1 \
    --repack_everything 0 \
    --ligand_mpnn_use_side_chain_context 1
```

The same idea applies for PyRosetta. Though since PyRosetta doesn't have a singular executable and instead requires python scripts, there isn't a specialized `pyrosetta` command. The `pyrosetta` environment still exists, of course, so you would instead run something like `pixi run -w PDchain -e pyrosetta python my_pyrosetta_script.py`. The following example has the python script (that imports pyrosetta) in the form of a heredoc, which also works as well.
```bash
pixi run -w PDchain -e pyrosetta python << 'EOF'
from pathlib import Path
import pyrosetta
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

pyrosetta.init('-relax:default_repeats 1')

xml_string = '''
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="r15" weights="ref2015.wts"/>
        <ScoreFunction name="r15_cst" weights="ref2015_cst.wts"/>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
        <True name="FullPose"/>
    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
        <ResfileCommandOperation name="relax_task" command="NATAA" residue_selector="FullPose"/>
    </TASKOPERATIONS>
    <MOVERS>
        <AddConstraints name="add_csts">
            <CoordinateConstraintGenerator name="ca_cst" sd="1" ca_only="1"/>
        </AddConstraints>
        <FastRelax name="relax" task_operations="relax_task" scorefxn="r15_cst"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover_name="add_csts"/>
        <Add mover_name="relax"/>
    </PROTOCOLS>
    <OUTPUT scorefxn="r15"/>
</ROSETTASCRIPTS>
'''
xml = XmlObjects.create_from_string(xml_string).get_mover('ParsedProtocol')

pose = pyrosetta.pose_from_pdb('inputs/1a53_clean.pdb')
xml.apply(pose)

Path('outputs').mkdir(parents=True, exist_ok=True)
pose.dump_pdb('outputs/ex0c_1a53_rlx.pdb')
EOF
```

## Built-In Metrics
There are a number of built-in metrics that are available for use. Note that secondary structure determination from `rfd_chain` is done with the DISICL algorithm (dihedral-only; no h-bonding interactions calculated).
|General Metrics|Description|
|:---:|:---|
|`totres`|total number of residues|
|`fin_reu`|total Rosetta Energy Units (REU); lower is better|
|`fin_reu_per_res`| total REU divided by the number of residues; lower is better|
|`protein_mpnn_score`| ProteinMPNN confidence score; higher is better; only scored with `--model_type protein_mpnn`|
|`soluble_mpnn_score`| SolubleMPNN confidence score; higher is better only scored with `--model_type soluble_mpnn`|
|`ligand_mpnn_score`| LigandMPNN confidence score; higher is better; only scored with `--model_type ligand_mpnn`|
|`perc_helix`|helix composition of protein scaffold|
|`perc_sheet`|sheet composition of protein scaffold|
|`perc_loop`|loop composition of protein scaffold|
|`fixedres_perc_helix`|helix composition of fixed residues +/- 3|
|`fixedres_perc_sheet`|sheet composition of fixed residues +/- 3|
|`fixedres_perc_loop`|loop composition of fixed residues +/- 3|
|`ddg`|ddg metric from Rosetta (complex only); lower is better|
|`dsasa`|dsasa metric from Rosetta (complex only)|
|`cst_rmsd`|root mean squared deviations from specified constraints; lower is better|
|`hbonds_to_lig_[resid]`|number of hbonds to `[resid]` calculated from Rosetta (protein-ligand complex only); follows the pdb numbering format (resi+chain) (e.g. `hbonds_to_lig_1X`)|

The following metrics are calculated with the poly-ala version of the design. These are meant to be sequence-agnostic. The term "focus residues" here refer to ligands and fixed residues.
|Poly-Ala Metrics|Description|
|:---:|:---|
|`clash`|clashes between focus residues and poly-ala scaffold; lower is better|
|`rog_ala`|approximate radius of gyration (heavy atoms only; all heavy atoms treated with equal weight); lower ~ compact and globular|
|`mp_dev`|standard deviation of distances between scaffold midpoint and all scaffold atoms; lower ~ hollow and globluar|
|`mp_fc_dst`|distance between the midpoint of focus residues and the midpoint of the rest of the poly-ala scaffold|
|`fc_compact`|root mean squared distance between focus midpoint and atoms of the 50 closest residues; lower ~ compact backbone around ligands; works best for compact ligands|
|`fc_dev`|standard deviation of distances between focus midpoint and atoms of the 50 closest residues; lower ~ hollow cavity around ligands; works best for compact ligands|

## Acknowledgements
PDchain was built on top of the following works:
* RosettaCommons/RFdiffusion
* YaoYinYing/RFdiffusion (mps- and cpu-compatible RFdiffusion)
* YaoYinYing/SE3Transformer (mps- and cpu-compatible SE3Transformer)
* dauparas/LigandMPNN
* RosettaCommons/rosetta
* facebookresearch/esm
* HeliXonProtein/OmegaFold
* schrodinger/pymol-open-source
* rdkit/rdkit
* matteoferla/rdkit-to-params
