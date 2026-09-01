#!/usr/bin/env bash
usage () { cat << EOF
Usage: gen_coord_csts PDB
Generate CoordinateConstraint Rosetta contraints

Parameters:
    PDB                     Input PDB

Options:
  --chri [str]              Residues to fix in chain+resi format
                            E.g. --chri A31-33 A71
  --resn [str]              Residues to fix in resn format
                            E.g. --resn ZN BMH
  --name [str]              Atoms to fix with name
                            E.g. --name N CA C
  --crd_stdev [float]       Strength/deviation of constraint
                            0.5 for strong, 2.0 for weak
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

val_opts=(chri resn name crd_stdev)
bool_opts=(hetatm help)
chri=""
resn=""
name=""
crd_stdev="0.5"

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

inpdb=${args[1]}

awk -v crd_stdev="$crd_stdev" '
    BEGIN {
        n = split("'"$chri"'", a, " ")
        for (i=1; i<=n; i++) {
            ch = a[i] ; gsub(/[^A-Z]/, "", ch)

            m = split(a[i], b, "-")
            gsub(/[^0-9]/, "", b[1])
            b[2] = b[2] ? b[2] : b[1]

            for (ri=b[1]+0; ri<=b[2]+0; ri++) {
                chri[ch ri] = 1
            }
        }

        n = split("'"$resn"'", a, " ")
        for (i=1; i<=n; i++) resn[a[i]] = 1

        n = split("'"$name"'", a, " ")
        for (i=1; i<=n; i++) name[a[i]] = 1
    }

    $1!~/^(ATOM|HETATM)$/ || $NF=="H" { next }

    $5 $6 in chri || $4 in resn || $3 in name {
        printf "CoordinateConstraint %s %s CA 1 ", $3, $6 $5
        printf "%s %s %s HARMONIC 0 %s\n", $7, $8, $9, crd_stdev
    }
' $inpdb
