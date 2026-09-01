# PDchain
Cross-platform command line tools for chaining protein design software (RFdiffusion, LigandMPNN, and PyRosetta) using a pixi workspace. Can iteratively optimize designs based on desired metrics with an in-silico continuous evolution system. Compatible with Linux, Mac, and Windows Subsystem for Linux. Supported platforms: osx-arm64, linux-64, linux-aarch64.

## Installation
Make sure pixi is installed on your command line. 

```
# Install pixi if you don't have it already
curl -fsSL https://pixi.sh/install.sh | sh
```

Follow the instructions provided in the terminal upon installation. You may need to restart your terminal and/or add ~/.pixi/bin to your PATH. See https://pixi.prefix.dev/latest/installation/ for more information.


Once you have pixi installed, run the following commands.

```
git clone ...                   # Download repository
cd PDchain                      # Navigate into root directory
pixi run install                # Install all environments

# Linux-aarch64 will be cpu-only, whereas osx-arm64 and linux-64 will make use of mps and cuda, respectively.
# If you want the entire environment to be cpu-only regardless of system, then run the following:
# cp box/tomls/cpu.toml pixi.toml && pixi run install
```

Installation will take some time (~10-20 minutes) as it will download all of the weights for RFdiffusion and LigandMPNN in addition to PyRosetta.

## Example one-off usage
```
# Make sandbox directory for quick tests
mkdir -p sandbox && cd sandbox

# Generate unconditional relaxed monomer 101 residues long
pixi run rfd_chain 101-101 --relax --outprefix outputs/sample1
```

Pixi is designed with workspaces in mind, so these `pixi run` commands will normally only work if you're inside the PDchain directory or its subdirectories. When outside the workplace, you'll need to specify the full path to `pixi.toml` with the `-m` option.

```
# Can run from any directory
pixi run -m /path/to/PDchain/pixi.toml rfd_chain 102-102 --relax --outprefix outputs/sample2
```

## Using the `pdchain` function (shortcut for pixi commands)
Finishing the installation process should spawn a script called `pdchain.sh` in the root directory. If you `source` this file, you will have access to the `pdchain` function, which is essentially a shortcut for pixi commands tailored for PDchain.

```
# Register the pdchain function
source /path/to/PDchain/pdchain.sh

# Help with pdchain function
pdchain --help

# Run pdchain from any directory
pdchain rfd_chain 103-103 --relax --outprefix outputs/sample3

# Activate the PDchain environment to directly use its command line tools from any directory
pdchain activate    # OR pdchain shell
rfd_chain 104-104 --relax --outprefix outputs/sample4
```

Using `pdchain activate` will allow the command line tools in the PDchain environment to temporarily supercede your existing tools. Because this environment has its own `bash`, `sed`, `grep`, `awk`, `coreutils` (`ls`, `cd`, `rm`, `mkdir`, etc.), and `pymol`, your commands will default to these versions upon activating the environment.

Pixi currently doesn't have a way to "deactivate" the environment in your current shell, so you'll have to close the terminal to escape the environment provided by `pdchain activate`. If you want a reversible interactive shell, then use `pdchain shell` instead. This pushes you into a subshell that you can `exit` out of. Only use `pdchain shell` interactively (never within a bash script like `~/.bashrc`, as it seems to cause pixi to be trapped in some infinite loop from my experience).

For those who are familiar with pixi, `pdchain shell` is just a shortcut for `pixi shell -m /path/to/PDchain/pixi.toml`, and `pdchain activate` is just a shortcut for `eval "$(pixi shell-hook -m /path/to/PDchain/pixi.toml)"`.

## Acknowledgements

The work here is made possible by the following software!

* RosettaCommons/RFdiffusion

* YaoYinYing/RFdiffusion (mps- and cpu-compatible RFdiffusion)

* YaoYinYing/SE3Transformer (mps- and cpu-compatible SE3Transformer)

* dauparas/LigandMPNN

* RosettaCommons/rosetta

* facebookresearch/esm

* schrodinger/pymol-open-source

* rdkit/rdkit

* matteoferla/rdkit-to-params

* prefix-dev/pixi
