#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.' >&2
    exit 1
}

usage () { cat << EOF
Usage: calc_molar_abs PDB
Calculates molar absorption coefficient of the input PDB
Also prints out the number of tryptophans, tyrosines, and cysteines

Parameters:
    PDB                     Input PDB

Options:
  --help                    Display this help and exit
EOF
}

args=($(echo $0 $@ | sed "s/--.*//"))
opts=($(echo $0 $@ | sed "s/--/\n--/" | sed "1d"))
[[ $((${#args[@]}-1)) -lt 1 ]] && help="1"
[[ $help == 1 ]] && { usage ; default_vals ; exit 1 ;}

inpdb=${args[1]}
awk '
    $1=="ATOM" && $3=="CA" && $4~/^(TRP|TYR|CYS)$/ { n[$4]++ }
    END {
        print "n_W " n["TRP"] + 0
        print "n_Y " n["TYR"] + 0
        print "n_C " n["CYS"] + 0
        print "molar_abs " (n["TRP"]+0) * 5500 + (n["TYR"]+0) * 1490 + (n["CYS"]+0) * 125
    }
' $inpdb
