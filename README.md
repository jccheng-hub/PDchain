# PDchain

> [!NOTE]
> This README is still a work in progress!

Command-line tools for chaining protein design software (RFdiffusion, LigandMPNN, and PyRosetta) using a pixi workspace. Can iteratively optimize designs based on desired metrics with an *in silico* continuous evolution system.

Supported platforms: osx-arm64, linux-64, linux-aarch64.

## Table of Contents
- [Installation](#installation)
- [Using the `pdchain` function (shortcut for pixi commands)](#using-the-optional-pdchain-function-shortcut-for-pixi-commands)
- [Usage Examples](#usage-examples)
  - [Unconditional Monomer Generation](#1---unconditional-monomer-generation)
  - [Protein Binder Design](#2---protein-binder-design)
    - [Protein Binder Redesign with Native Backbone](#2a---protein-binder-redesign-with-native-backbone)
    - [Protein Binder Redesign with Partial Diffusion](#2b---protein-binder-redesign-with-partial-diffusion)
    - [Protein Binder Redesign with Indels](#2c---protein-binder-redesign-with-indels)
    - [Protein Binder De Novo Design](#2d---protein-binder-de-novo-design)
  - [Ligand Binder / Enzyme Design](#3---ligand-binder--enzyme-design)
    - [Ligand Binder / Enzyme Redesign with Native Backbone](#3a---ligand-binder--enzyme-redesign-with-native-backbone)
    - [Ligand Binder / Enzyme Redesign with Partial Diffusion](#3b---ligand-binder--enzyme-redesign-with-partial-diffusion)
    - [Ligand Binder / Enzyme Redesign with Indels](#3c---ligand-binder--enzyme-redesign-with-indels)
    - [Ligand Binder / Enzyme De Novo Design](#3d---ligand-binder--enzyme-de-novo-design)

## Installation

Run the following command to install the newest version of pixi.
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
Restart your terminal or source your shell's rc file (e.g. `source ~/.bashrc`) to complete the installation. You can check if pixi is registered with `command -V pixi`.

See https://pixi.prefix.dev/latest/installation/ for more information.

> [!IMPORTANT]
> This repository's installation script will automatically pull PyRosetta as a dependency. While the original code in this repository is open-source, PyRosetta is not free for commercial use (free for academic, non-profit, and government institutions). Please ensure you are not violating PyRosetta's terms of service by having the appropriate license before running this repository's installation script.

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

> [!NOTE]
> The CUDA version listed in the pixi.toml file for the Linux GPU installation of RFdiffusion is 11.8. Some newer GPUs from NVIDIA are incompatible with CUDA 11.8 and requires newer versions of CUDA. If this is the case, then the installation will fail despite having a CUDA-compatible GPU. To get around this, you can attempt to fix the pixi.toml file to run on a newer version of CUDA. This may require additional adjustments like changing the python/pytorch/dgl versions, shifting to installation from pypi instead of conda-forge, and/or changing the available channels.

Installation will take some time (~10-20 minutes) as it will download all of the weights for RFdiffusion and LigandMPNN in addition to PyRosetta.

## Using the optional `pdchain` function (shortcut for pixi commands)

Normally, when you want to run a command from a pixi workspace, you would need to prefix your command with `pixi run`. For example:
```bash
pixi run rfd_chain --help
```

But this only works if you are within the pixi workspace's subdirectories. If you want to access that workspace regardless of your working directory path, you would have to specify either (1) the path to pixi.toml or (2) the workspace name if it has one.
```bash
pixi run -m /path/to/PDchain/pixi.toml rfd_chain --help
pixi run -w PDchain rfd_chain --help
```

This can be a bit verbose, so PDchain will spawn a script called `pdchain.sh` upon finishing the installation process. If you `source` this file, you will register the `pdchain` function, which is essentially a shortcut that allows you to run PDchain commands from anywhere with minimal typing.
```bash
source pdchain.sh
```

If you have registered the `pdchain` function, you can use `pdchain` as a prefix to access commands from within PDchain. For example...
```bash
pdchain rfd_chain --help
```

In addition to acting as a shortcut, you also have access to some additional commands.
```bash
pdchain list      # Lists available PDchain-specific commands
pdchain activate  # Modify current shell for direct access to PDchain's default environment
pdchain shell     # Create a new shell for direct access to PDchain's default environment
```

If you ran `pdchain activate` or `pdchain shell`, this will enable you to run PDchain commands without any prefixes.
```bash
rfd_chain --help
```

Using `pdchain activate` or `pdchain shell` will allow the command line tools in the PDchain environment to temporarily supercede your existing tools. Because this environment has its own `bash`, `sed`, `grep`, `awk`, `coreutils` (`ls`, `cd`, `rm`, `mkdir`, etc.), and `pymol`, your commands will default to these versions upon activating the environment.

Under the hood, `pdchain shell` is just a shortcut for `pixi shell -m /path/to/PDchain/pixi.toml`, and `pdchain activate` is just a shortcut for `eval "$(pixi shell-hook -m /path/to/PDchain/pixi.toml)"`.

If you don't already have a preferred molecular viewing software, PDchain's default environment does come with the open-source version of PyMOL, which you can access in a number of ways as mentioned above, summarized here.
```bash
pixi run -m /path/to/PDchain/pixi.toml pymol  # Run from anywhere
pixi run -w PDchain pymol                     # Run from anywhere (if PDchain is registered as the workspace name)
pixi run pymol                                # Run from PDchain subdirectories
pdchain pymol                                 # Run from anywhere (if `source pdchain.sh` was ran)
pymol                                         # Run from anywhere after `pdchain activate` or `pdchain shell`
```

For more information on the `pdchain` function, run the following.
```bash
pdchain --help
```

## Usage Examples
Example scripts can be found in the `examples` directory. All example scripts are written as if they would be ran with `examples` as the working directory.

### 1 - Unconditional Monomer Generation
```bash
pixi run -w PDchain rfd_chain 140-160 \
    --idealize --relax --ca_stdev 1 \
    --model_type protein_mpnn \
    --design_cycles 3 \
    --numdes 3 \
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

### 2 - Protein Binder Design

The following command Rosetta refines the input protein-protein complex, which in this case is the barnase-barstar complex, without touching the backbone or sequence. This generates a control (no design done) that we can use as a reference for the actual design runs later in this section.

```bash
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
```

The `SKIP` keyword at where the contigs is supposed to be tells `rfd_chain` to skip RFdiffusion.

`--natbias 10` applies a biasing weight of 10 towards the native (input) residues at every position during MPNN sequence design. A weight of 10 basically forces MPNN to recover the input residues at every position, preventing any actual sequence design but still allowing MPNN to rebuild/repack the side chains.

### 2a - Protein Binder Redesign with Native Backbone

The following command will diversify the binder sequence. RFdiffusion is not applied here (just MPNN-FastRelax). The command is identical to the control command except it does not have the `--natbias 10` option.

```bash
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
```

`--fixedres B1-110` is necessary here to prevent MPNN from sequence designing the target protein (the barnase on chain B).

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN confidence scores will be accepted after refinement, but this flag will make it so that acceptance/rejection depends solely on the specified metric. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section **Built-in Metrics (WIP)** for available metrics.

### 2b - Protein Binder Redesign with Partial Diffusion
The following command will diversify the binder (the barstar on chain A) with partial diffusion before applying cycles of MPNN-FastRelax.

```bash
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
```

`--model_ckpt Complex_base` specifies the model checkpoint used for RFdiffusion. `Complex_base` is recommended for protein binder design.

`--partial` turns on partial diffusion. Note that the output diffused structure will always match the input structure in length with partial diffusion.

`--timesteps 1` specifies the number of timesteps for RFdiffusion. When using partial diffusion, this number is allowed to dip below 15. The higher number of timesteps, the greater the deviation from the original input structure.

### 2c - Protein Binder Redesign with Indels

The following command will diversify the binder by rediffusing loop regions while allowing for insertions and deletions.

```bash
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
    --ss_to_contigs --vary_linkers 1 --ss_trim 2-3 \
    --outprefix outputs/ex2c_1brs_indel
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. Loop residues will be masked from RFdiffusion.

`--vary_linkers 1` will allow the loop regions to vary in length by 1 residue. So if a loop region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

`--ss_trim 2-3` will allow 2-3 helix/sheet residues neighboring loop residues to be masked from RFdiffusion. This gives RFdiffusion more wiggle room when filling in those masked regions.

The indel diversification approach is more computationally expensive than partial diffusion because it requires a minimum of 15 timesteps during RFdiffusion, but the additional diversity it provides can be beneficial depending on the design goal.

### 2d - Protein Binder De Novo Design

The following command will generate de novo protein binders.

```bash
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
    --ppi_hotspots B27 B38 B54-59 B82-85 B101-104 \
    --outprefix outputs/ex2d_1brs_denovo
```

Here, we provided the contigs `86-95/0 B1-110` to specify that we want to generate a backbone 86-95 residues long while preserving our target (chain B residues 1-110).

`--ppi_hotspots B27 B38 B54-59 B82-85 B101-104` provides hotspots for RFdiffusion, and it will attempt to generate a backbone near those specified residues.

De novo design often requires generating thousands of designs and computationally screening through them to get something "reasonable". If we compare the designs to the original control, we are likely to see designs that actually perform worse on many desirable metrics. For example, if we were to check the ddg of the outputs compared to control, (e.g. with `grep -H '^ddg ' outputs/ex2*.pdb`), we are likely to see the control outperform most if not all of the de novo designs. While barnase-barstar is an incredibly tight complex (so the bar [HA!] is pretty high here), large scale generation and screening will still often be necessary to have high confidence in your designs.

Because we are unlikely to get excellent designs right away, it is often necessary to take an agreeable-but-not-excellent de novo design and diversify around that with indels and partial diffusion to get something better. This is the main reason why the *in silico* continuous evolution system in PDchain was built. See the section on `evo_rfd_chain` for more details.

### 3 - Ligand Binder / Enzyme Design
The following example takes a protein complexed with a ligand (in this case, a Kemp eliminase complexed with its transition-state analog) and refines it with Rosetta (without any sequence or backbone design). This command is meant to generate a computational control to serve as a reference point for later design runs.

```bash
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
```

Setting the contigs to `SKIP` skips RFdiffusion entirely, and `--natbias 10` forces MPNN to recover the native sequence. Together, these prevent any backbone or sequence design.

### 3a - Ligand Binder / Enzyme Redesign with Native Backbone

The following command will redesign the sequence using the native backbone (no RFdiffusion done here). The command is identical to that of the control in the previous section except without `--natbias 10`.

```bash
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
```

`--fixedres A50 A127` fixes the the residues 50 and 127 on chain A. In this example, they correspond to an Asp and Gln that hydrogen-bond to the ligand 6NT.

`--ligname 6NT` specifies the ligand from the input structure you want to keep. Without specifying the ligand names, `rfd_chain` will ignore it entirely.

`--model_type soluble_mpnn` specifies that the use of the SolubleMPNN model, which was not trained to be ligand-aware. When you specify a `--model_type` that isn't `ligand_mpnn` but you include a ligand, the model you specified will be applied first, then the residues within 8 angstroms of the ligand will be redesigned with LigandMPNN. This ensures that at least the residues surrounding the ligand are being redesigned in a ligand-conscious manner.

`--lig_stdev 0.5` applies Rosetta coordinate constraints to the ligand. This one was done with a standard deviation of 0.5, which makes it stronger than the CA coordinate constraints specified with `--ca_stdev 1` in this example.

`--ap_stdev 0.5` applies Rosetta distance (AtomPair) constraints between the ligand and fixed residues. This option only works when both `--fixedres` and `--ligname` are specified. These constraints will preserve the relative geometry between the ligand and the fixed residues (the Asp on A50 and the Gln on A127 in this case).

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN confidence scores will be accepted after refinement, but this flag will make it so that acceptance/rejection depends solely on the specified metric. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section **Built-in Metrics (WIP)** for available metrics.


### 3b - Ligand Binder / Enzyme Redesign with Partial Diffusion

The following example will noise/denoise the input backbone with partial diffusion before MPNN-FastRelax.

```bash
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
```

`--partial` enables partial diffusion, and `--timesteps 1` sets the number of noising/denoising steps to 1. Increase the value for `--timesteps` if you want more diversity from the original input.

### 3c - Ligand Binder / Enzyme Redesign with Indels

The following example will preserve helix and sheet residues while allowing loop regions to rediffuse with varying lengths, allowing for insertions and deletions.

```
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
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. Loop residues will be masked from RFdiffusion.

`--vary_linkers 1` will allow the loop regions to vary in length by 1 residue. So if a loop region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

`--ss_trim 2-3` will allow 2-3 helix/sheet residues neighboring loop residues to be masked from RFdiffusion. This gives RFdiffusion more wiggle room when filling in those masked regions.

The indel diversification approach is more computationally expensive than partial diffusion because it requires a minimum of 15 timesteps during RFdiffusion, but the additional diversity it provides can be beneficial depending on the design goal.

### 3d - Ligand Binder / Enzyme De Novo Design

> [!NOTE]
> This repository integrates the original RFdiffusion, not RFdiffusion2 or RFdiffusion3. While RFdiffusion2 and RFdiffusion3 supports ligand binder / enzyme design without the need to specify starting backbone coordinates or residue indices, the original RFdiffusion does not. In other words, using RFdiffusion for this design task is not the most elegant strategy given the existence of its successors. Still, it is entirely possible to scaffold active sites with the original RFdiffusion, and this example will show how it can be done in the context of this pixi workspace.

The first thing we need is a set of inverse rotamers. The following example shows how we can generate them.

```bash
pixi run -w PDchain gen_invrots inputs/5rgf_clean.pdb A50:N3 A127:N3 X1 \
    --numrots 10 \
    --parallel 1 \
    --outdir outputs/ex3d_invrots
```

The command `gen_invrots` will randomly generate inverse rotamers from the input PDB (given as the first argument). The second argument and onward will tell `gen_invrots` exactly what residues to keep. In this example, we are keeping residues A50 (catalytic Asp), A127 (catalytic Gln), and X1 (ligand).

If a residue specified by `gen_invrots` is an amino acid, you have the option of placing those residues on secondary structures. This follows the format of `[residue]:[H/E/N][int]`, where H, E, and N represent helix, sheet, and any (either helix or sheet), and the integer that follows specifies the number of residues to add at each end of that residue. So for example, our `A50:N3` will scaffold residue A50 on either a 7-residue beta-strand or a 7-residue helix with the residue of interest (Asp on A50) being right in the middle.

`--numrots 10` specify the number of inverse rotamers you want to generate. Once the number of generated rotamers exceed this value, `gen_invrots` will stop generating more.

`--parallel 1` specifies the number of parallel jobs you want to run simultaneously to generate inverse rotamers. Here, we only spawn one job at a time, but increasing this value will make generating inverse rotamers faster. This should not exceed the number of CPU cores you have on your system.

`--outdir outputs/ex3d_invrots` specify where we want to store the inverse rotamers. This directory will serve as the input for the actual design run.

Once you have generated your inverse rotamers, you can take a peak at them in pymol to see if they are reasonable with `pixi run -w PDchain pymol $(echo outputs/ex3d_invrots/*.pdb | cut -d' ' -f1-10)`. Notice that the chain letter code is now different. While A50 is still on A50, A127 is now on B127, and X1 (the ligand) is now on C1. This is because `gen_invrots` has to split each input residue onto its own chain to prevent potential overlap in both residue number and residue chain when the additional residues for secondary structures are added. As a result, we need to be make sure we specify the new chain+residue ids in the subsequent design command.

The following command will generate a de novo binder to the target ligand using the inverse rotamers we previously generated.

```
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
```

`--indir outputs/ex3d_invrots` specifies that we want the input to be a directory, and that input directory is where our previously generated inverse rotamers are located `outputs/ex3d_invrots`. When you specify an input directory (with `--indir`) as opposed to an input pdb (with `--inpdb`), a random .pdb file from that input directory is chosen as the input pdb.

As mentioned previously, what used to be A127 is now B127 because of `gen_invrots` placing each input fixed residue on its own chain.

`--fixbbres A49-51 B126-128` specifies what residues we want fixed during RFdiffusion. Given that our main residues are on A50 and B127, this example here retains an extra residue on both ends for both A50 and B127. Note that we could have gone up to the range of A47-53 and B124-130 because of how we generated those inverse rotamers (7-residue peptide scaffolding each key residue).

`--mintot 180 --addtot 40` provides the overall length of the desired monomer. Here, we are requesting a monomer 180-220 residues long (minimum 180 with a maximum addition of 40). This option in conjunction with `--fixbbres` will trigger `rfd_chain` to run `gen_contigs` to generate a contigs string based on these constraints. You still have the option to specify your own contigs as the first argument if you want more specific control.

Note that the order in which the residues are specified for `--fixbbres` matters. The command `gen_contigs`, which `rfd_chain` automatically accesses when a contigs string is not provided, will add residues in between the fixed regions *in the specified order* as well as additional residues at the N and C termini. In this case, the output will follow the order of `[diffused]-[A49-51]-[diffused]-[B126-128]-[diffused]` in which `[diffused]` represent the backbone regions filled in by RFdiffusion. Had we specified `--fixbbres B126-128 A49-51` instead, the order between the two fixed regions would have been flipped (i.e. `[diffused]-[B126-128]-[diffused]-[A49-51]-[diffused]`).

`--model_ckpt ActiveSite` tells RFdiffusion to use the ActiveSite checkpoint, which is tailored for scaffolding small motifs. This is more necessary here as we are essentially scaffolding small peptides floating in space. This isn't as necessary for the partial diffusion or indel examples above, but it was kept in those examples too for consistency.

> [!NOTE]
> RFdiffusion generates backbones around ligands with an auxillary potential (see RFdiffusion documentation for more details). However, even with this potential, RFdiffusion can produce backbones that clash with the input ligand. The command `rfd_chain` was therefore programmed to exclude backbones with excess clashes to ligand atoms. As a result, you may not see the same amount of design outputs as intended. For example, we have `--numdes 3` in the command above, but you may end up with just `ex3d_5grf_denovo_0002.pdb` in the outputs without the first and third designs because they failed the clash checker. If you want to guarantee the exact number of outputs as specified by `--numdes`, add the option `--persistent` to make RFdiffusion try again if it fails. If you want to turn off the clash checker entirely, you can set the threshold to an excessively high number like `--clash_cut 9999` so that every backbone RFdiffusion generates will always pass regardless of ligand clashes. This is not recommended, however, as MPNN may get an unreasonable input and Rosetta may have a hard time resolving those clashes.

As mentioned in the example with de novo protein binder design, it is highly unlikely to get any reasonable de novo designs from a small-scale computational run. Getting promising design candidates often require generating thousands of designs and screening through them. Sometimes, the best designs from even a large-scale batch can have undesirable properties, such low secondary structure component or a high radius of gyration. These designs can then be subjected to redesign with partial diffusion or indels for further optimization. See the section on `evo_rfd_chain` to see how one can evolve designs based on desirable metrics.

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
