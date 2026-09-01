#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

usage () { cat << EOF
Usage: bash install.sh TYPE
Installation script for PDchain

Parameters:
    TYPE                    Type of installation (cpu or gpu)

Options:
  --scriptdirs [str]        Script directories in box to add as CLI tools
  --help                    Display this help and exit
EOF
}

# --- Defaults --- # {{{
args=($(echo $0 $@ | sed "s/--.*//"))
opts=($(echo $0 $@ | sed "s/--/\n--/" | sed "1d"))
optchk () {
    [[ " ${opts[*]} " =~ " $1 " ]] && echo 1 || echo 0
}
optarg () {
    sed "s/--/\n--/g" <<< "${opts[*]}" |
    sed -n "1,/^$1 /s/^$1 //p" | sed 's/ \+$//'
}

val_opts=(scriptdirs)
bool_opts=(help)
scriptdirs="protocols tools"

default_vals () {
    [[ ${#val_opts[@]} -ge 1 ]] && {
        echo -e "\nDefaults:"
        for opt in ${val_opts[@]} ; do
            [[ -n ${!opt} ]] && echo "  --${opt} ${!opt}" || echo "  --${opt} NONE"
        done | sort | xargs printf "  %-19s %-19s %-19s %-s\n"
    }
    [[ ${#bool_opts[@]} -ge 1 ]] && {
        echo -e "\nFlags:"
        sed 's/ /\n/g' <<< "${bool_opts[@]}" |
        sed 's/^/--/g' | xargs printf "  %-39s %-s\n"
    }
}

all_opts=()
for opt in ${bool_opts[@]} ; do
    [[ $(optchk --$opt) == 1 ]] && all_opts+=(--$opt)
    declare "$opt=$(optchk --$opt)"
done
for opt in ${val_opts[@]} ; do
    val=$(optarg --${opt})
    [[ -n $val ]] && declare "$opt=$val"
    [[ ${!opt} == REQUIRED ]] && help="1"
    [[ -n ${!opt} ]] && all_opts+=(--$opt ${!opt})
done
[[ $((${#args[@]}-1)) -lt 1 ]] && help="1"
[[ $help == 1 ]] && { usage ; default_vals ; exit ;}
# }}}

set -e
shopt -s nullglob

# Download models
mkdir -p $pixiroot/models
modlinks=(
    http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
    http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
    http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
    http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
    http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
    http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
    http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_002.pt 
    https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_010.pt 
    https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_020.pt 
    https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_030.pt 
    https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_005_25.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_010_25.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_020_25.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_030_25.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/per_residue_label_membrane_mpnn_v_48_020.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/global_label_membrane_mpnn_v_48_020.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_002.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_010.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_020.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_030.pt
    https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_sc_v_32_002_16.pt
)
for link in ${modlinks[@]} ; do
    target="$pixiroot/models/$(basename $link)"
    [[ ! -e $target ]] && wget -nc -P $pixiroot/models $link
done

# Tweak RFdiffusion
rfddir="$pixiroot/box/programs/RFdiffusion"
pkgdir=$(find $pixiroot/.pixi/envs/rfdiffusion -name "site-packages" | head -n1)
[[ -e $pkgdir ]] && {
    echo "Adjusting RFdiffusion environment..."
    rsync -ah --info=progress2 --mkpath $rfddir/examples $pkgdir
    rsync -ah --info=progress2 --mkpath $pixiroot/models $pkgdir
    rsync -ah --info=progress2 --mkpath $rfddir/config $pixiroot/.pixi/envs/rfdiffusion
}
echo

# Add RFdiffusion symlink
rfdlink="$pixiroot/.pixi/envs/rfdiffusion/bin/rfdiffusion"
if [[ $1 == "cpu" ]] ; then
    cp $rfddir/scripts/run_inference_cpu.py $rfddir/scripts/run_inference.py
    echo "$rfddir/scripts/run_inference.py is set up for CPU use only."
elif [[ $1 == "gpu" ]] ; then
    cp $rfddir/scripts/run_inference_gpu.py $rfddir/scripts/run_inference.py
    echo "$rfddir/scripts/run_inference.py is set up for GPU use (mps or cuda)."
fi
ln -sfn $pixiroot/box/programs/RFdiffusion/scripts/run_inference.py $rfdlink &&
echo "Added symlink: $rfdlink"

# Add LigandMPNN symlink
mpnnlink="$pixiroot/.pixi/envs/ligandmpnn/bin/ligandmpnn"
ln -sfn $pixiroot/box/programs/LigandMPNN/run.py $mpnnlink &&
echo "Added symlink: $mpnnlink"

# Install omegafold github
pixi run -m $pixiroot/pixi.toml -e omegafold pip install --no-deps git+https://github.com/HeliXonProtein/OmegaFold.git

# Add box protocols and tools to bin
boxdir="$pixiroot/box"
bindir="$pixiroot/.pixi/envs/default/bin"
cmdtools=()
for dir in $scriptdirs ; do
    scripts=($boxdir/$dir/*.sh)
    for script in ${scripts[@]} ; do
        bn=$(basename ${script%.*})
        ln -sfn $script $bindir/$bn
        cmdtools+=("$bn")
    done
done

# Report CLI tools
echo -e "\nThe following CLI tools were added:"
cmdtools_f=$(
    sed 's/ /\n/g' <<< ${cmdtools[*]} | sort |
    xargs -n4 printf "    %-19s %-19s %-19s %-s\n"
)
echo "$cmdtools_f"

# Print functions
sed -E 's/^ {4}//' << EOF > $pixiroot/pdchain.sh
    _pdchain_usage() { cat << EOF
    usage: pdchain [CMD]
    Shortcut for running tools in the PDchain workspace.
    Replace [CMD] with any PDchain command line tools.

    For example:
    pdchain rfd_chain 100-100 --outprefix outputs/example0 --relax
    pdchain mpros --help

    CMD special options:
        activate                Activates default PDchain environment within the
                                current shell. Allows direct input commands without
                                pdchain prefix.

                                E.g. pdchain activate
                                rfd_chain 100-100 --outprefix outputs/example1

                                Keep in mind that the PDchain environment has GNU
                                coreutils, so commands such as ls, cd, mkdir, etc.
                                will default to them. This is true for bash, sed,
                                grep, awk, and pymol as well.

        shell                   Enters subshell with PDchain environment. Allows
                                direct input commands just like with "activate".
                                Because you're in a subshell, non-exported functions
                                and variables will not be inherited. You can exit
                                the subshell with the "exit" command.

                                WARNING: Do not put this in your ~/.bashrc! From my
                                experience (as of March 2026), it gets trapped in
                                some infinite loop. Might be a bug with pixi.

        list                    List available command line tools from PDchain.

        --help                  Show this message

    If you want to run RFdiffusion, LigandMPNN, etc. directly, you must specify the
    environment name with -e, then specify the command with the same name.
    E.g. pdchain -e rfdiffusion rfdiffusion [...]
    E.g. pdchain -e ligandmpnn ligandmpnn [...]
    EOF
    }
    
    _pdchain_list() {
        echo "List of available command line tools:"
        echo "${cmdtools[*]}" | tr ' ' '\n' | sort | xargs -n4 printf "    %-19s %-19s %-19s %-s\n"
    }

    pdchain() {
        _pdchain_toml="$pixiroot/pixi.toml"
    
        if [ "\$1" = "activate" ]; then
            eval "\$(pixi shell-hook -m "\$_pdchain_toml")"
        elif [ "\$1" = "shell" ]; then
            pixi shell -m "\$_pdchain_toml"
        elif [ "\$1" = "list" ]; then
            _pdchain_list
        elif [ "\$1" = "--help" ] || [ "\$#" -eq 0 ]; then
            _pdchain_usage
        else
            pixi run -m "\$_pdchain_toml" "\$@"
        fi
    }

    if [ -n "\$BASH_VERSION" ] ; then
        export -f pdchain _pdchain_usage _pdchain_list
    fi
EOF

echo "
Run the following command to register pdchain as a function

    source $pixiroot/pdchain.sh
"
