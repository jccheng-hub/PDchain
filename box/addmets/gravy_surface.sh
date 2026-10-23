#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.' >&2
    exit 1
}

usage () { cat << EOF
Usage: gravy_surf PDB
Calculates GRAVY (grand average of hydropathy) from Kyte & Doolittle, 1982 of
the input PDB.

Also calculates GRAVY weighted by relative solvent accessibility (RSA) using
maximum theoretical solvent accessibility values from Tien et. al, 2013.

Also calculates GRAVY for only surface residues (RSA >= 0.2).

This is meant to be a proxy for crystallographic tendency. Higher values of
gravy_rsa and gravy_surface are desirable, but overall gravy ideally remains
negative to prevent aggregation (in theory, anyways).

Output format:
gravy [value]
gravy_rsa [value]
gravy_surface [value]

Parameters:
    PDB                     Input PDB

Options:
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

val_opts=(

)

bool_opts=(
    help
)

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
[[ $help == 1 ]] && { usage ; default_vals ; exit 1 ;}
# }}}

inpdb=${args[1]}

tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX)
trap 'rm -rf $tmpdir' EXIT

python << EOF
from pymol import cmd
cmd.load('$inpdb', 'oriprot')
cmd.remove('not polymer')          # Only keep protein
cmd.remove('elem H')               # Strip hydrogens
cmd.set('dot_solvent', 1)          # SASA
cmd.set('dot_density', 4)          # Highest accuracy
cmd.set('solvent_radius', 1.4)     # Ensure default value of 1.4
cmd.get_area('oriprot', load_b=1)  # Calculate surface area, then load into b_factor col
cmd.save('$tmpdir/tmp.pdb')        # Saved output has SASA in b_factor column
EOF

awk '
$1=="ATOM" { asa[$5][$6]+=$11 ; res[$5][$6] = $4 }
END {
    max_asa_arr(mxa) ; gravy_arr(grv)

    n = n_surface = 0
    for (ch in asa) for (ri in asa[ch]) {
        if (res[ch][ri] == "") continue
        rsa = asa[ch][ri] / mxa[res[ch][ri]]
        s_gravy += grv[res[ch][ri]]
        s_gravy_rsa += grv[res[ch][ri]] * rsa
        n++

        if (rsa < 0.2) continue
        s_gravy_surface += grv[res[ch][ri]]
        n_surface++
    }

    print "gravy " s_gravy / n
    print "gravy_rsa " s_gravy_rsa / n
    print "gravy_surface " s_gravy_surface / n_surface
}

# awk functions {{{
function max_asa_arr(arr) {
    # Initialize array for maximum accessible surface area (MaxASA)
    # From Tien et al., 2013 ; this is the *theoretical* set
    # Format: arr["three-letter-code"] = max-asa
    delete arr
    arr["ALA"] = 129.0 ; arr["LEU"] = 201.0
    arr["ARG"] = 274.0 ; arr["LYS"] = 236.0
    arr["ASN"] = 195.0 ; arr["MET"] = 224.0
    arr["ASP"] = 193.0 ; arr["PHE"] = 240.0
    arr["CYS"] = 167.0 ; arr["PRO"] = 159.0
    arr["GLU"] = 223.0 ; arr["SER"] = 155.0
    arr["GLN"] = 225.0 ; arr["THR"] = 172.0
    arr["GLY"] = 104.0 ; arr["TRP"] = 285.0
    arr["HIS"] = 224.0 ; arr["TYR"] = 263.0
    arr["ILE"] = 197.0 ; arr["VAL"] = 174.0
}

function gravy_arr(arr) {
    # Initialize array for GRAVY
    # From Kyte & Doolittle, 1982
    # Format: arr["three-letter-code"] = hydropathicity-index
    delete arr
    arr["ALA"] = 1.80  ; arr["LEU"] = 3.80 
    arr["ARG"] = -4.50 ; arr["LYS"] = -3.90
    arr["ASN"] = -3.50 ; arr["MET"] = 1.90 
    arr["ASP"] = -3.50 ; arr["PHE"] = 2.80 
    arr["CYS"] = 2.50  ; arr["PRO"] = -1.60
    arr["GLN"] = -3.50 ; arr["SER"] = -0.80
    arr["GLU"] = -3.50 ; arr["THR"] = -0.70
    arr["GLY"] = -0.40 ; arr["TRP"] = -0.90
    arr["HIS"] = -3.20 ; arr["TYR"] = -1.30
    arr["ILE"] = 4.50  ; arr["VAL"] = 4.20 
}
# }}}
' $tmpdir/tmp.pdb
