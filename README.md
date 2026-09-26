# PDchain

> [!NOTE]
> This README is still a work in progress!

Command-line tools for chaining protein design software (RFdiffusion, LigandMPNN, and PyRosetta) using a pixi workspace. Can iteratively optimize designs based on desired metrics with an *in silico* continuous evolution system.

Supported platforms: osx-arm64, linux-64, linux-aarch64.

## Table of Contents
- [Installation](#installation)
- [Usage Examples](#usage-examples)
  - [Unconditional Monomer Generation](#1---unconditional-monomer-generation)
  - [Protein Binder Design](#2---protein-binder-design)
    - [Protein Binder Redesign with Partial Diffusion](#2a---protein-binder-redesign-with-partial-diffusion)
    - [Protein Binder Redesign with Indels](#2b---protein-binder-redesign-with-indels)
    - [Protein Binder De Novo Design](#2c---protein-binder-de-novo-design)

## Installation
If your system doesn't have pixi already, run the following command to install it. 
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
cp box/tomls/cpu.toml pixi.toml && pixi run install && pixi workspace register --name PDchain
```

If you have an Apple Silicon Mac or a CUDA-compatible Linux/WSL and you want to leverage GPU acceleration for RFdiffusion, run the following instead.

```bash
cp box/tomls/gpu.toml pixi.toml && pixi run install && pixi workspace register --name PDchain
```

Installation will take some time (~10-20 minutes) as it will download all of the weights for RFdiffusion and LigandMPNN in addition to PyRosetta.

## Usage Examples
Example scripts can be found in the `examples` directory. All example scripts are written as if they would be ran from with `examples` as the working directory.

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

The following command Rosetta refines the input PDB (barnase-barstar complex), which already has a protein-protein interaction. This generates a *computational control* (no design done) that we can use as a reference for later design protocols.

```bash
pixi run -w PDchain rfd_chain NONE \
    --inpdb inputs/1brs_af3mod0.pdb \
    --idealize --relax \
    --model_type protein_mpnn --natbias 10 \
    --ca_stdev 1 \
    --design_cycles 3 \
    --select_met min:ddg \
    --numdes 1 \
    --outprefix outputs/ex2_1brs_control
```

The `NONE` keyword at where the contigs is supposed to be tells `rfd_chain` to skip RFdiffusion.

`--natbias 10` applies a biasing weight of 10 towards the native (input) residues at every position during MPNN sequence design. A weight of 10 basically forces MPNN to recover the input residues at every position, preventing any actual sequence design but still allowing MPNN to rebuild/repack the side chains.

Rosetta Idealize and FastRelax refinement is then applied after the MPNN sequence "design", and the output will serve a computational control.

### 2a - Protein Binder Redesign with Partial Diffusion
The following command will diversify the binder (the barstar on chain A) with partial diffusion before applying cycles of MPNN sequence design + Rosetta FastRelax.

```bash
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
```

`--partial` turns on partial diffusion. Note that the output diffused structure will always match the input structure in length with partial diffusion.

`--fixedres B1-110` is necessary here to prevent MPNN from sequence designing the target protein (the barnase on chain B).

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score and/or MPNN confidence scores will be accepted after refinement, but this flag will make it so that acceptance/rejection depends solely on the specified metric. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section **ADD SECTION HERE** for available metrics.

### 2b - Protein Binder Redesign with Indels

The following command will diversify the binder by rediffusing loop regions while allowing for insertions and deletions.

```bash
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
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. Loop residues will be masked from RFdiffusion.

`--vary_linkers 1` will allow the loop regions to vary in length by 1 residue. So if a loop region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

`--ss_trim 1-3` will allow 1-3 helix/sheet residues neighboring loop residues to be masked from RFdiffusion. This gives RFdiffusion more wiggle room when filling in those masked regions.

The indel diversification approach is more computationally expensive than partial diffusion because it requires a minimum of 15 timesteps during RFdiffusion, but the additional diversity it provides can be beneficial depending on the design goal.

### 2c - Protein Binder De Novo Design

The following command will generate de novo protein binders.

```bash
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
```

Here, we provided the contigs `86-95/0 B1-110` to specify that we want to generate a backbone 86-95 residues long while preserving our target (chain B residues 1-110).

`--ppi_hotspots A73 A75` provides hotspots for RFdiffusion, and it will attempt to generate a backbone near those specified residues.

De novo design often requires generating thousands of designs and computationally screening through them to get something "reasonable". If we compare the designs to the original control, we are likely to see designs that actually perform worse on many desirable metrics. For example, if we were to check the ddg of the outputs compared to control, (e.g. with `grep -H '^ddg ' outputs/ex2*.pdb`), we are likely to see the control outperform most if not all of the de novo designs. While barnase-barstar is an incredibly tight complex (so the bar [HA!] is pretty high here), large scale generation and screening will still often be necessary to have high confidence in your designs.

Because we are unlikely to get excellent designs right away, it is often necessary (at least from my experience), to take an agreeable-but-not-excellent de novo design and diversify around that with indels and partial diffusion to get something better, and that is the main reason why this *in silico* continuous evolution system here was built. See the section on `evo_rfd_chain` for more details.

## Using the `pdchain` function (shortcut for pixi commands)
Finishing the installation process should spawn a script called `pdchain.sh` in the root directory. If you `source` this file, you will have access to the `pdchain` function, which is essentially a shortcut for pixi commands tailored for PDchain.

```bash
# Register the pdchain function
source /path/to/PDchain/pdchain.sh

# List available CLI tools from PDchain
pdchain list

# Help with pdchain function
pdchain --help

# Run from any directory
pdchain rfd_chain 103-103 --relax --outprefix outputs/sample3

# Activate the PDchain environment to directly use its command line tools from any directory
pdchain activate    # OR pdchain shell
rfd_chain 104-104 --relax --outprefix outputs/sample4
```

Using `pdchain activate` will allow the command line tools in the PDchain environment to temporarily supercede your existing tools. Because this environment has its own `bash`, `sed`, `grep`, `awk`, `coreutils` (`ls`, `cd`, `rm`, `mkdir`, etc.), and `pymol`, your commands will default to these versions upon activating the environment.

Pixi currently doesn't have a way to "deactivate" the environment in your current shell, so you'll have to close the terminal to escape the environment provided by `pdchain activate`. If you want a reversible interactive shell, then use `pdchain shell` instead. This pushes you into a subshell that you can `exit` out of. Only use `pdchain shell` interactively (not within a bash script like `~/.bashrc`, as it seems to cause pixi to be trapped in some infinite loop from my experience).

For those who are familiar with pixi, `pdchain shell` is just a shortcut for `pixi shell -m /path/to/PDchain/pixi.toml`, and `pdchain activate` is just a shortcut for `eval "$(pixi shell-hook -m /path/to/PDchain/pixi.toml)"`.

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
