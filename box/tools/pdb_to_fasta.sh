#!/usr/bin/env bash
usage () { cat << EOF
Usage: pdb_to_fasta PDB...
Make fasta from PDBs

Parameters:
    PDB                     Input PDB file

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

pdb_to_fasta () {
    local inpdbs=$@
    awk -v OFS=',' '
    BEGIN {
        olc["ALA"]="A" ; olc["ARG"]="R" ; olc["ASN"]="N" ; olc["ASP"]="D"
        olc["CYS"]="C" ; olc["GLU"]="E" ; olc["GLN"]="Q" ; olc["GLY"]="G"
        olc["HIS"]="H" ; olc["ILE"]="I" ; olc["LEU"]="L" ; olc["LYS"]="K"
        olc["MET"]="M" ; olc["PHE"]="F" ; olc["PRO"]="P" ; olc["SER"]="S"
        olc["THR"]="T" ; olc["TRP"]="W" ; olc["TYR"]="Y" ; olc["VAL"]="V"
    }
    
    FNR == 1 { ch = $5 }
   
    $1~/^(ATOM|HETATM)$/ && $3=="CA" && $NF!="CA" {
        if (ch != $5) { seq = seq ":" ; ch = $5 }
        aa = olc[$4] ? olc[$4] : "X"
        seq = seq aa
    }

    ENDFILE {
        bn = FILENAME
        sub(/^.*[/]/, "", bn)
        sub(/[.]pdb$/, "", bn)
        print ">" bn
        print seq
        bn = seq = ""
    }
    ' $inpdbs
}

pdb_to_fasta ${args[@]:1}
