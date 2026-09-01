#!/usr/bin/env bash
usage () { cat << EOF
Usage: track_lineage DESIGN
Track the lineage for the specified design (from evo_rfd_chain)

Parameters:
    DESIGN                  Path to design

Options:
  --help                    Display this help and exit
EOF
}

# >>> Defaults >>> {{{
args=($(echo $0 $@ | sed "s/--.*//"))
opts=($(echo $0 $@ | sed "s/--/\n--/" | sed "1d"))
optchk () {
    [[ " ${opts[*]} " =~ " $1 " ]] && echo 1 || echo 0
}
optarg () {
    sed "s/--/\n--/g" <<< "${opts[*]}" |
    sed -n "1,/^$1 /s/^$1 //p" | sed 's/ \+$//'
}

val_opts=()
bool_opts=(help)

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
# <<< Defaults <<< }}}

# Check Pixi
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

design=(${args[1]})
pooldir=$(dirname $design)
archive="$pooldir/archive"

lineage=($(
    {
        ls "$design"
        bn=$(sed -n 's/^inpdb //p' $design)
        parent=$(find $pooldir $archive -maxdepth 1 -name "${bn}.pdb" 2>/dev/null | head -n1)
        while ls "$parent" 2>/dev/null ; do
            bn=$(sed -n 's/^inpdb //p' $parent)
            parent=$(find $pooldir $archive -maxdepth 1 -name "${bn}.pdb" 2>/dev/null | head -n1)
        done
    } | tac
))

loadcmd=$(echo "${lineage[@]}" | sed 's/ /\n/g' | sed 's/^/load /')
first=$(basename $lineage | cut -d. -f1)
pymol -d "
$loadcmd
join_states lineage, all, -2
disable not lineage
spectrum count, rainbow, lineage and name CA
set movie_fps, 5
set movie_loop, 0
"
