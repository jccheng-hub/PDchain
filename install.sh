#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

usage () { cat << EOF
Usage: bash install.sh [CPU/GPU]
Installation script for PDchain

Options:
  --scriptdirs [str]        Directories in box with scripts to be added as CLI
                            tools. The directories 'protocols' and 'tools' are
                            included by default.
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
scriptdirs=""

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

# Declare type of installation
if [[ $1 =~ ^CPU|GPU$ ]] ; then
    echo "Setting up PDchain for $1 usage..."
else
    echo "Need to specify CPU or GPU." && exit 1
fi

# Install omegafold github
pixi run -m $pixiroot/pixi.toml -e omegafold pip install --no-deps git+https://github.com/HeliXonProtein/OmegaFold.git

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

# Add RFdiffusion symlink
rfddir="$pixiroot/box/programs/RFdiffusion"
rfdlink="$pixiroot/.pixi/envs/rfdiffusion/bin/rfdiffusion"
rfd_py="$pixiroot/.pixi/envs/rfdiffusion/bin/run_inference.py"
ln -sfn $rfd_py $rfdlink &&
echo "Added symlink: $rfdlink"

# Add LigandMPNN symlink
mpnnlink="$pixiroot/.pixi/envs/ligandmpnn/bin/ligandmpnn"
ln -sfn $pixiroot/box/programs/LigandMPNN/run.py $mpnnlink &&
echo "Added symlink: $mpnnlink"

# Add box protocols and tools to bin
boxdir="$pixiroot/box"
bindir="$pixiroot/.pixi/envs/default/bin"
cmdtools=()
for dir in protocols tools $scriptdirs ; do
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
sed -E 's/^ {4}//' << EOF > $bindir/pdchain-info && chmod +x $bindir/pdchain-info
    #!/usr/bin/env bash
    cat << 'EOF'
    If you are running this command, you should already be in the PDchain default
    environment. In addition to the tools listed below, you will also have access to
    the open-source version of PyMOL. This default environment is what allows all
    the command-line tools in PDchain (which are all essentially bash scripts) to
    run. This also means that this command-line is now populated with GNU coreutils,
    GNU awk, GNU sed, etc., as those were used in the making of PDchain commands.

    Example usage, such as protein binder design, enzyme design, and continuous
    design evolution are provided in the README.md, which is best viewed on the
    github: https://github.com/jccheng-hub/PDchain.git

    EOF
    echo "List of available command-line tools:"
    echo "${cmdtools[*]}" | tr ' ' '\n' | sort | xargs -n4 printf "    %-19s %-19s %-19s %-s\n"
EOF
pixi shell-hook -m $pixiroot/pixi.toml > $pixiroot/pdchain.sh

echo "
Run the following command to activate the PDchain default environment:
source $pixiroot/pdchain.sh
"
