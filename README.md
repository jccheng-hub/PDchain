# PDchain
An accessible cross-platform pixi workspace with command-line tools for chaining protein design software (RFdiffusion, LigandMPNN, and PyRosetta). Can iteratively optimize designs based on desired metrics with an *in silico* continuous evolution system.

Supported platforms: osx-arm64, linux-64, linux-aarch64

A GPU workspace setup is available for Macs and CUDA-compatible Linux (or Windows Subsystem for Linux) systems. For Linux or WSL systems without CUDA compatibility, CPU-only setups are also available.

## Table of Contents
- [Installation](#installation)
- [Running PDchain Commands](#running-pdchain-commands)
- [Unconditional Monomer Generation](#unconditional-monomer-generation)
- [Protein Binder Design](#protein-binder-design)
  - [Cyclic Peptide Design](#cyclic-peptide-de-novo-design)
- [Ligand Binder / Enzyme Design](#ligand-binder--enzyme-design)
- [Continuous Evolution of Designs](#continuous-evolution-of-designs)
- [Fold Validation Methods](#fold-validation-methods)
- [Directly Running RFdiffusion, LigandMPNN and PyRosetta](#directly-running-rfdiffusion-ligandmpnn-and-pyrosetta)
- [Built-In Metrics](#built-in-metrics)
- [Disclosure on AI usage](#disclosure-on-ai-usage)
- [Acknowledgements](#acknowledgements)

## Installation
> [!IMPORTANT]
This repository's installation script will automatically pull PyRosetta as a dependency. While the original code in this repository is open-source, PyRosetta is not free for commercial use (free for academic, non-profit, and government institutions). Please ensure you are not violating PyRosetta's terms of service by having the appropriate license before running this repository's installation script.

Run the following command to install the newest version of pixi.
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
Restart your terminal or source your shell's rc file (e.g. `source ~/.bashrc`) to complete the installation. You can check if pixi is registered with `command -V pixi`.

See https://pixi.prefix.dev/latest/installation/ for more information.

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

The `pixi workspace register --name PDchain` command will register `PDchain` as the workspace name. This allows you to use the `-w PDchain` option for various pixi commands. More on this in the [Running PDchain Commands](#running-pdchain-commands) section.

Installation will take some time (~10-20 minutes) as it will download all of the weights for RFdiffusion and LigandMPNN in addition to PyRosetta.

## Running PDchain Commands

You can enter a *subshell* with PDchain's default environment with *one* of the following `pixi shell` commands:
```bash
pixi shell                                  # Must be in PDchain subdirectory
pixi shell -w PDchain                       # Works anywhere if PDchain is a registered workspace name
pixi shell -m /path/to/PDchain/pixi.toml    # Works anywhere
```
You can always exit out of your PDchain subshell with the command `exit`.

If you want to incorporate PDchain's default environment in your *current shell*, run *one* of the following instead:
```bash
source pdchain.sh                     # pdchain.sh should have spawned in PDchain's root directory upon installation
eval "$(pixi shell-hook -w PDchain)"  # Works anywhere if PDchain is a registered workspace name
```

Generally, booting up a subshell is better for interactive use, and incorporating into current shell is better for scripting use.

Once inside the PDchain default environment, you can now directly access a variety of command-line tools in the PDchain default environment. For example:
```bash
pdchain-info            # Report available tools on PDchain
rfd_chain --help        # Help documentation for rfd_chain
esmfold_relax --help    # Help documentation for esmfold_relax
```

If you want to run PDchain commands without having to change/switch environments whatsoever, you can prefix your command with `pixi run`; any of the following approaches can work.
```bash
pixi run rfd_chain --help                                  # Only works in PDchain subdirectory
pixi run -w PDchain rfd_chain --help                       # Works anywhere if PDchain is a registered workspace name
pixi run -m /path/to/PDchain/pixi.toml rfd_chain --help    # Works anywhere
```

If you don't mind using the free and open source version of PyMOL, it is available in PDchain's default environment. You can simply run `pymol` (if you're in the environment) or `pixi run -w PDchain pymol` (if you're outside the environment).

### Running Example Commands
All example commands in this README are meant to be executed in PDchain's default environment from the `examples` directory, and the files containing example inputs are provided at `examples/inputs.zip`. Before running the example code blocks, run `pixi shell -w PDchain` (to enter the PDchain default environment), `cd examples` (to change the working directory to `examples`), and `unzip inputs.zip` (to extract the input files).

## Unconditional Monomer Generation
```bash
rfd_chain 140-160 \
    --numdes 3 \
    --outprefix outputs/ex1_rand
```

Here, we run `rfd_chain` to randomly generate monomeric proteins with RFdiffusion, sequence design with ProteinMPNN, and refine the structure with Rosetta FastRelax.

Run `rfd_chain --help` to display the manual on the command line.

The input argument `140-160` specifies the contigs string. Here we tell RFdiffusion that we want to generate a backbone 140 to 160 residues long. See the RFdiffusion github documentation for more information on contigs.

`--numdes 3` specifies that we want a total of three designs.

`--outprefix outputs/ex1_rand` specifies the path prefix of the output designs.

## Protein Binder Design

### Control Refinement

The following command Rosetta refines the input protein-protein complex, which in this case is the barnase-barstar complex, without touching the backbone or sequence. This generates a control (no design done) that we can use as a reference for the actual design runs later in this section.

```bash
rfd_chain SKIP \
    --inpdb inputs/1brs_af3mod0.pdb \
    --skip_mpnn \
    --numdes 1 \
    --select_met min:ddg \
    --outprefix outputs/1brs_control
```

The `SKIP` keyword at where the contigs is supposed to be tells `rfd_chain` to skip RFdiffusion, and `--skip_mpnn` skips the MPNN sequence design step. The result is a computational control refined by Rosetta.

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score will be accepted after refinement, but this flag will make it so that acceptance depends on the specified metric instead. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section on [Built-In Metrics](#built-in-metrics) for available metrics.

Note that in this control case, no MPNN sequence design is actually done, but the iterative rounds of FastRelax (followed by acceptance based on `min:ddg`) will still run.

### Protein Binder Redesign with Partial Diffusion
The following command will noise/denoise the backbone of the binder (the barstar on chain A) with partial diffusion before applying cycles of MPNN-FastRelax.

```bash
rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --numdes 3 \
    --select_met min:ddg \
    --partial --timesteps 3 \
    --outprefix outputs/1brs_partial
```

`--fixedres B1-110` fixes residues 1-110 on chain B (barnase, in this case). This flag prevents the specified residues from being rediffused by RFdiffusion and redesigned by MPNN.

`--partial` turns on partial diffusion. Note that the output diffused structure will always match the input structure in length with partial diffusion.

`--timesteps 3` specifies the number of timesteps for RFdiffusion. When using partial diffusion, this number is allowed to dip below 15. The higher number of timesteps, the greater the deviation from the original input structure.

### Protein Binder Redesign with Linker Regeneration

The following command will diversify the binder by regenerating the linker residues connecting stretches of secondary structure while allowing for insertions and deletions.

```bash
rfd_chain \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --numdes 3 \
    --select_met min:ddg \
    --ss_to_contigs --ss_trim 1 --vary_linkers 1 \
    --outprefix outputs/1brs_linker
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. This invokes linker regeneration. The default linkers are loop residues.

`--ss_trim 1` will allow helix/sheet residues neighboring loop residues to be treated as additional linker residues (which are then masked from RFdiffusion). Sheets are trimmed by the specified value, and helices are trimmed by twice the specified value.

`--vary_linkers 1` will allow the linker regions to vary in length by 1 residue. So if a linker region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

### Protein Binder De Novo Design

The following command will generate de novo protein binders.

```bash
rfd_chain 86-95/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --numdes 3 \
    --select_met min:ddg \
    --ppi_hotspots B27 B38 B54-59 B82-85 B101-104 \
    --outprefix outputs/1brs_denovo
```

Here, we provided the contigs `86-95/0 B1-110` to specify that we want to generate a backbone 86-95 residues long while preserving our target (chain B residues 1-110).

`--ppi_hotspots B27 B38 B54-59 B82-85 B101-104` provides hotspots for RFdiffusion, and it will attempt to generate a backbone near those specified residues. These hotspots were chosen in this example due to them forming contacts with the wild-type binder barstar.

De novo design often requires generating thousands of designs and computationally screening through them to get something "reasonable". If we compare the designs to the original control, we are likely to see designs that actually perform worse on many desirable metrics. For example, if we were to check the ddg of the outputs compared to control, (e.g. with `grep -H '^ddg ' outputs/ex2*.pdb`), we are likely to see the control outperform most if not all of the de novo designs. While barnase-barstar is an incredibly tight complex (so the bar [HA!] is pretty high here), large scale generation and screening will still often be necessary to have high confidence in your designs.

Because we are unlikely to get excellent designs right away, it is often necessary to take an agreeable-but-not-excellent de novo design and diversify around that with indels and partial diffusion to get something better. This is the main reason why the *in silico* continuous evolution system in PDchain was built. See the section on [Continuous Evolution of Designs](#4---continuous-evolution-of-designs) for more details.

### Cyclic Peptide De Novo Design

The following command will generate a *cyclic* peptide to bind a target protein. It's the same idea as the protein binder example but just with the additional option `--cyclize_chain A`.
```bash
rfd_chain 8-12/0 B1-110 \
    --inpdb inputs/1brs_af3mod0.pdb \
    --fixedres B1-110 \
    --numdes 3 \
    --select_met min:ddg \
    --ppi_hotspots B27 B38 B54-59 B82-85 B101-104 \
    --cyclize_chain A \
    --outprefix outputs/1brs_cycpep
```

Note that running PDchain's version of RFdiffusion will always assign the de novo backbone as chain A in the final protein complex. Should your input PDB have the *target* protein on chain A, that chain A of your input PDB will be shifted to chain B while the new chain A is reserved for the de novo backbone. Because the chain you specify with this option refers to the chain of the eventual complex, you will always want to specify `--cyclize_chain A` when generating de novo cyclic peptides  *regardless of the chain ID of the target protein* in the input PDB.

## Ligand Binder / Enzyme Design

### Control Refinement

The following example takes a protein complexed with a ligand (in this case, a Kemp eliminase complexed with its transition-state analog) and refines it with Rosetta (without any sequence or backbone design). This command is meant to generate a computational control to serve as a reference point for later design runs.

```bash
rfd_chain SKIP \
    --inpdb inputs/5rgf_clean.pdb \
    --numdes 1 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --skip_mpnn \
    --outprefix outputs/5rgf_control
```

The `SKIP` keyword at where the contigs is supposed to be tells `rfd_chain` to skip RFdiffusion, and `--skip_mpnn` skips the MPNN sequence design step. The result is a computational control refined by Rosetta.

`--select_met min:ddg` specifies the selection metric during iterative rounds of MPNN-FastRelax. By default, designs with improved Rosetta score will be accepted after refinement, but this flag will make it so that acceptance depends on the specified metric instead. Here, the `min:` prefix specifies that we want lower values of `ddg`. If you want to maximize some metric value instead, you would use the `max:` prefix (e.g. `--select_met max:protein_mpnn_score`). If you want to lean towards some specific values, you would use the `val` prefix followed by `=[desired_value]` (e.g. `--select_met val:dsasa=0.7`). See section on [Built-In Metrics](#built-in-metrics) for available metrics.

Note that in this control case, no MPNN sequence design is actually done, but the iterative rounds of FastRelax (followed by acceptance based on `min:ddg`) will still run.

The combination of the four following flags preserves active site geometry during refinement:
- `--fixedres A50 A127` fixes residues 50 and 127 on chain A.
- `--ligname 6NT` specifies the ligand(s) to preserve from the input structure. 
- `--lig_stdev 0.5` constrains ligand atoms to their starting coordinates (with CoordinateConstraint) with a 0.5 Å standard deviation.
- `--ap_stdev 0.5` constrains distances (with AtomPair) between all ligand atoms and all fixed residue atoms with a 0.5 Å standard deviation. This only works if both `--fixedres` and `--ligname` are invoked.

### Ligand Binder / Enzyme Redesign with Partial Diffusion

The following example will noise/denoise the input backbone with partial diffusion before MPNN-FastRelax.

```bash
rfd_chain \
    --inpdb inputs/5rgf_clean.pdb \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --sc_context \
    --rfd_ckpt ActiveSite \
    --partial --timesteps 3 \
    --outprefix outputs/5rgf_partial
```

`--partial` enables partial diffusion, and `--timesteps 3` sets the number of noising/denoising steps to 3. Increase the value for `--timesteps` if you want more diversity from the original input.

`--rfd_ckpt ActiveSite` specifies the ActiveSite checkpoint for RFdiffusion. This is typically recommended for scaffolding smaller motifs. This is optional here.

`--sc_context` treats side chain atoms of fixed residues as extra ligand atoms. This presumably makes LigandMPNN sensitive to rotamers.

Note that by default, the ProteinMPNN model is used (i.e. the option `--mpnn_model protein_mpnn` is assumed). If ligands are specified while using a non-LigandMPNN model, an additional round of MPNN redesign for residues within 8 Å of each ligand will be conducted with LigandMPNN, ensuring that the final sequence was ligand-aware.

### Ligand Binder / Enzyme Redesign with Linker Regeneration

The following example will diversify the scaffold by regenerating the linker residues connecting stretches of secondary structure while allowing for insertions and deletions.

```bash
rfd_chain \
    --inpdb inputs/5rgf_clean.pdb \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 A127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --sc_context \
    --rfd_ckpt ActiveSite \
    --ss_to_contigs --ss_trim 1 --vary_linkers 1 \
    --outprefix outputs/5rgf_linker
```

`--ss_to_contigs` will generate a contigs string based on the secondary structure of the input PDB. This invokes linker regeneration. The default linkers are loop residues.

`--ss_trim 1` will allow helix/sheet residues neighboring loop residues to be treated as additional linker residues (which are then masked from RFdiffusion). Sheets are trimmed by the specified value, and helices are trimmed by twice the specified value.

`--vary_linkers 1` will allow the linker regions to vary in length by 1 residue. So if a linker region was originally 5 residues long, that region can end up 4-6 residues long in the output design.

### Ligand Binder / Enzyme De Novo Design

> [!NOTE]
This repository integrates the original RFdiffusion, not its successors, due to its cross-platform compatibility. However, RFdiffusion successors (e.g. RFdiffusion3) support ligand binder / enzyme design without the need to specify starting backbone coordinates or residue indices, and the original RFdiffusion does not. In other words, using RFdiffusion for this design task is not the most elegant strategy given the existence of its successors. Still, it is entirely possible to scaffold active sites with the original RFdiffusion, and this example will show how it can be done in the context of this pixi workspace.

The first thing we need is a set of inverse rotamers. The following example shows how we can generate them.

```bash
gen_invrots inputs/5rgf_clean.pdb A50:E2 A127:E2 X1 \
    --numrots 10 \
    --parallel 1 \
    --outdir outputs/invrots
```

The command `gen_invrots` will randomly generate inverse rotamers from the input PDB (given as the first argument). The second argument and onward will tell `gen_invrots` exactly what residues to keep. In this example, we are keeping residues A50 (catalytic Asp), A127 (catalytic Gln), and X1 (ligand).

If a residue specified by `gen_invrots` is an amino acid, you have the option of placing those residues on secondary structures. This follows the format of `[residue]:[H/E/N][int]`, where H, E, and N represent helix, sheet, and any (either helix or sheet), and the integer that follows specifies the number of residues to add at each end of that residue. So for example, `A50:E2` will scaffold residue A50 on a 5-residue beta-strand with the residue of interest (Asp on A50) being right in the middle.

`--numrots 10` specify the number of inverse rotamers you want to generate. Once the number of generated rotamers exceed this value, `gen_invrots` will stop generating more.

`--parallel 1` specifies the number of parallel jobs you want to run simultaneously to generate inverse rotamers. Here, we only spawn one job at a time, but increasing this value will make generating inverse rotamers faster. This should not exceed the number of CPU cores you have on your system.

`--outdir outputs/invrots` specify where we want to store the inverse rotamers. This directory will serve as the input for the actual design run.

Once you have generated your inverse rotamers, you can take a peak at them in pymol to see if they are reasonable with `pymol $(echo outputs/invrots/*.pdb | cut -d' ' -f1-10)`. Notice that the chain letter code is now different. While A50 is still on A50, A127 is now on B127, and X1 (the ligand) is now on C1. This is because `gen_invrots` has to split each input residue onto its own chain to prevent potential overlap in both residue number and residue chain when the additional residues for secondary structures are added. As a result, we need to be make sure we specify the new chain+residue ids in the subsequent design command.

The following command will generate a de novo binder to the target ligand using the inverse rotamers we previously generated.

```bash
rfd_chain \
    --indir outputs/invrots \
    --numdes 3 \
    --select_met min:ddg \
    --fixedres A50 B127 --ligname 6NT \
    --lig_stdev 0.5 --ap_stdev 0.5 \
    --sc_context \
    --rfd_ckpt ActiveSite \
    --fixbbres A49-51 B126-128 --tot_range 180-220 \
    --outprefix outputs/5rgf_denovo
```

`--indir outputs/invrots` specifies that we want the input to be a directory, and that input directory is where our previously generated inverse rotamers are located. When you specify an input directory (with `--indir`) as opposed to an input pdb (with `--inpdb`), a random .pdb file from that input directory is chosen as the input pdb.

As mentioned previously, what used to be A127 is now B127 because of `gen_invrots` placing each input fixed residue on its own chain.

`--fixbbres A49-51 B126-128` specifies what residues we want fixed during RFdiffusion. Given that our main residues are on A50 and B127, this example here retains an extra residue on both ends for both A50 and B127. Note that we could have gone up to the range of A48-52 and B125-129 because of how we generated those inverse rotamers (5-residue peptide per key residue).

`--tot_range 180-220` provides the overall length of the desired monomer. Here, we are requesting a monomer 180-220 residues long. This option in conjunction with `--fixbbres` will trigger `rfd_chain` to run `gen_contigs` to generate a contigs string based on these constraints. You still have the option to specify your own contigs as the first argument if you want more specific control.

Note that the order in which the residues are specified for `--fixbbres` matters. The command `gen_contigs`, which `rfd_chain` automatically accesses when a contigs string is not provided, will add residues in between the fixed regions *in the specified order* as well as additional residues at the N and C termini. In this case, the output will follow the order of `[diffused]-[A49-51]-[diffused]-[B126-128]-[diffused]` in which `[diffused]` represent the backbone regions filled in by RFdiffusion. Had we specified `--fixbbres B126-128 A49-51` instead, the order between the two fixed regions would have been flipped (i.e. `[diffused]-[B126-128]-[diffused]-[A49-51]-[diffused]`).

`--rfd_ckpt ActiveSite` tells RFdiffusion to use the ActiveSite checkpoint, which is tailored for scaffolding small motifs. This is more necessary here as we are essentially scaffolding small peptides floating in space. This isn't as necessary for the partial diffusion or the linker regeneration examples above, but it was kept in those examples too for consistency.

> [!NOTE]
RFdiffusion generates backbones around ligands with an auxillary potential (see RFdiffusion documentation for more details). However, even with this potential, RFdiffusion can produce backbones that clash with the input ligand. The command `rfd_chain` was therefore programmed to exclude backbones with excess clashes to ligand atoms. As a result, you may not see the same amount of design outputs as intended. For example, we have `--numdes 3` in the command above, but you may end up with just `ex3d_5grf_denovo_0002.pdb` in the outputs without the first and third designs because they failed the clash checker. If you want to guarantee the exact number of outputs as specified by `--numdes`, add the option `--persistent` to make RFdiffusion try again if it fails. If you want to turn off the clash checker entirely, you can set the threshold to an excessively high number like `--clash_cut 9999` so that every backbone RFdiffusion generates will always pass regardless of ligand clashes. This is not recommended, however, as MPNN may get an unreasonable input and Rosetta may have a hard time resolving those clashes.

It is unlikely to get "reasonable" de novo designs from a small-scale computational run. Getting promising design candidates often require generating thousands of designs and screening through them. Sometimes, even the best designs from a large-scale batch might not meet all of your desired criteria. For example, they might be globular with ample secondary structure composition, but the binding pocket is completely buried/exposed. In these situations, it might be better to redesign these candidates to optimize for those desired properties rather than to repeatedly fish for new designs that meet all of your criteria all at once. See the section on [Continuous Evolution of Designs](#4---continuous-evolution-of-designs) to see how one can evolve designs based on desirable metrics.

## Continuous Evolution of Designs
> [!NOTE]
> This section is under construction!

Not every design output will possess properties you are looking for; this is especially true for de novo designs. Often times, you end up with designs that check some boxes but not others. In these situations, it may be worth attempting optimization of these designs by using them as starting points for cycles of diversification and selection.

...

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

This should spawn a json file in the output directory. This json can then be uploaded to [alphafoldserver.com](https://alphafoldserver.com) as job drafts, which can then be subsequently submitted so long as you have enough jobs left for the day.

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

... you would instead run...
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
There are a number of built-in metrics that are available for use. Here, they are categorized into sequence-agnostic and sequence-dependent metrics. The sequence-agnostic metrics are calculated with the poly-ala version of the design. Generally, sequence-agnostic metrics are best employed when diffusing in major chunks of the protein backbone (e.g. de novo or major rediffusion), and sequence-dependent metrics are best employed when the backbone undergoes minor changes (e.g. partial diffusion or fixed-backbone redesign).

Note that secondary structure determination from `rfd_chain` is done with the DISICL algorithm (dihedral-only; no h-bonding interactions calculated). Also, the term "focus residues" below refer to ligands and fixed residues.

|Sequence-Agnostic (Poly-Ala) Metrics|Description|
|:---:|:---|
|`totres`|total number of residues|
|`perc_helix`|helix composition of protein scaffold|
|`perc_sheet`|sheet composition of protein scaffold|
|`perc_loop`|loop composition of protein scaffold|
|`fixedres_perc_helix`|helix composition of fixed residues +/- 3|
|`fixedres_perc_sheet`|sheet composition of fixed residues +/- 3|
|`fixedres_perc_loop`|loop composition of fixed residues +/- 3|
|`clash`|clashes between focus residues and poly-ala scaffold; lower is better|
|`rog_ala`|approximate radius of gyration (heavy atoms only; all heavy atoms treated with equal weight); lower ~ compact and globular|
|`rog_ala_A`|approximate radius of gyration (heavy atoms only; all heavy atoms treated with equal weight) for chain A; lower ~ compact and globular|
|`rog_ala_B`|approximate radius of gyration (heavy atoms only; all heavy atoms treated with equal weight) for chain B; lower ~ compact and globular|
|`mp_dev`|standard deviation of distances between scaffold midpoint and all scaffold atoms; lower ~ hollow and globluar|
|`mp_fc_dst`|distance between the midpoint of focus residues and the midpoint of the rest of the poly-ala scaffold|
|`fc_compact`|root mean squared distance between focus midpoint and atoms of the 50 closest residues; lower ~ compact backbone around ligands; works best for compact ligands|
|`fc_dev`|standard deviation of distances between focus midpoint and atoms of the 50 closest residues; lower ~ hollow cavity around ligands; works best for compact ligands|

|Sequence-Dependent Metrics|Description|
|:---:|:---|
|`fin_reu`|total Rosetta Energy Units (REU); lower is better|
|`fin_reu_per_res`| total REU divided by the number of residues; lower is better|
|`protein_mpnn_score`| ProteinMPNN confidence score; higher is better; only scored with `--model_type protein_mpnn`|
|`soluble_mpnn_score`| SolubleMPNN confidence score; higher is better only scored with `--model_type soluble_mpnn`|
|`ligand_mpnn_score`| LigandMPNN confidence score; higher is better; only scored with `--model_type ligand_mpnn`|
|`ddg`|ddg metric from Rosetta (complex only); lower is better|
|`dsasa`|dsasa metric from Rosetta (complex only)|
|`cst_rmsd`|root mean squared deviations from specified constraints; lower is better|
|`netcharge`|net charge metric from Rosetta|
|`hbonds_to_lig_[resid]`|number of hbonds to `[resid]` calculated from Rosetta (protein-ligand complex only); follows the pdb numbering format (resi+chain) (e.g. `hbonds_to_lig_1X`)|

## Disclosure on AI usage
AI was used by the author as a resource for learning a range of programming techniques (many of which were leveraged when making this repository). All of the scripts in PDchain (e.g. scripts in `box/protocols` and `box/tools`) were hand coded; none were AI generated nor modified from AI-generated work.

## Acknowledgements
PDchain was built on top of the following works:
* [RosettaCommons/RFdiffusion](https://github.com/RosettaCommons/RFdiffusion.git)
* [YaoYinYing/RFdiffusion (mps- and cpu-compatible RFdiffusion)](https://github.com/YaoYinYing/RFdiffusion)
* [YaoYinYing/SE3Transformer (mps- and cpu-compatible SE3Transformer)](https://github.com/YaoYinYing/SE3Transformer)
* [dauparas/LigandMPNN](https://github.com/dauparas/LigandMPNN)
* [RosettaCommons/rosetta](https://github.com/RosettaCommons/rosetta)
* [facebookresearch/esm](https://github.com/facebookresearch/esm)
* [HeliXonProtein/OmegaFold](https://github.com/HeliXonProtein/OmegaFold)
* [schrodinger/pymol-open-source](https://github.com/schrodinger/pymol-open-source)
* [rdkit/rdkit](https://github.com/rdkit/rdkit)
* [matteoferla/rdkit-to-params](https://github.com/matteoferla/rdkit-to-params)
