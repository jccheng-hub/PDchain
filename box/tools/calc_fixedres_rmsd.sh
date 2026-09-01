#!/usr/bin/env bash
usage () { cat << EOF
Usage: calc_fixedres_rmsd REF:INPDB [REF:INPDB...]
Use PyMOL to calculate fixed residue RMSD between reference (REF) and input PDB.
Useful for comparing designs with their predicted structures while taking only
critical residues into account.

Requirements:
  - The "fixedres" line that specifies the list of residues to align with should
    be in the reference PDB. E.g. "fixedres A1 A3 A7 A20"
  - The chain and numbering of both reference and input PDBs must be identical.

Parameters:
    REF                     Path to reference PDB
    INPDB                   Path to input PDB

Options:
  --atom_names [str]        Atom names to use for alignment. Default is all
                            heavy atoms.
                            E.g. --atom_names N CA C O
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

val_opts=(atom_names)
bool_opts=(help)
atom_names=""

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

# Initialize
inpairs="${args[@]:1}"

[[ -n $atom_names ]] && echo "Calculating RMSD with only following atoms: $atom_names"

# Initialize python dictionary
pydict=$(mktemp ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX.txt) ; trap 'rm $pydict' EXIT
cat << EOF > $pydict
pdbpaths = {}
fixedres = {}
EOF

for inpair in $inpairs ; do
    unset p f

    p=($(sed 's/:/ /' <<< "$inpair"))
    if [[ ${#p[@]} -lt 2 ]] ; then
        echo "Could not parse $inpair. Needs format of [REF:INPDB]."
        continue
    elif ! grep '^ATOM' "${p[0]}" &> /dev/null ; then
        echo "Reference PDB ${p[0]} does not exist or does not have ATOM lines."
        continue
    elif ! grep '^ATOM' "${p[1]}" &> /dev/null ; then
        echo "Input PDB ${p[1]} does not exist or does not have ATOM lines."
        continue
    fi

    f=($(sed -n 's/^fixedres //p' ${p[0]}))
    if [[ -z $f ]] ; then
        echo "Did not find fixedres in the reference PDB ${p[0]}."
        continue
    fi

    f_exp=$(
        sed 's/ /\n/g' <<< ${f[*]} |
        sed -E 's|^([A-Z]+)([0-9]+)$|(\1/\2/'"${atom_names// /+}"')|' |
        paste -sd '+'
    )

    grep '^fixedres_rmsd ' ${p[1]} &> /dev/null && {
        sed -i '/^fixedres_rmsd /d' ${p[1]} &&
        echo "Removed fixedres_rmsd value in input PDB ${p[1]}."
    }

    echo "Aligning ${p[1]} to ${p[0]} with ${f[*]}"
    echo "pdbpaths['${p[0]}'] = '${p[1]}'" >> $pydict
    echo "fixedres['${p[0]}'] = '$f_exp'" >> $pydict
done

# Run PyMOL in python
define_pydict=$(cat $pydict)
python << EOF
from pymol import cmd
$define_pydict

for refpdb in pdbpaths:
    cmd.load(refpdb, 'refpdb')
    cmd.load(pdbpaths[refpdb], 'inpdb')
    cmd.remove('elem H')
    aln = cmd.align(
        f'inpdb and ({fixedres[refpdb]})',
        f'refpdb and ({fixedres[refpdb]})',
        cycles=0
    )
    print(f'RMSD between {refpdb} and {pdbpaths[refpdb]}: {aln[0]:.4f}')

    # Append info
    with open(pdbpaths[refpdb], 'a') as file:
        file.write(f'fixedres_rmsd {aln[0]:.4f}\n')
        
    cmd.delete('refpdb')
    cmd.delete('inpdb')
EOF
