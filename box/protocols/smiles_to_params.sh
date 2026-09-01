#!/usr/bin/env bash
usage () { cat << EOF
Usage: smiles_to_params SMILES TLC
Convert SMILES string to Rosetta params using rdkit_to_params

Parameters:
    SMILES                  Input smiles string
    TLC                     Three-letter code for ligand

Options:
  --forcefield [str]        Forcefield to use
                            Options: UFF, MMFF
  --outdir [str]            Output directory
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

val_opts=(forcefield outdir)
bool_opts=(help)

forcefield="UFF"
outdir="."

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

smiles=${args[1]}
tlc=${args[2]}
outdir=${outdir%/}
mkdir -p $outdir

$pixirun rdkit python << EOF
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params

# Make mol object from smiles
mol = Chem.MolFromSmiles('$smiles')

if mol is None:
    raise ValueError('Failed to load molecule. Check input file.')

mol = Chem.AddHs(mol)    # Add hydrogens

# Make 3D structure, then optimize
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.${forcefield}OptimizeMolecule(mol)

# Compute charges, then save
AllChem.ComputeGasteigerCharges(mol)
Chem.MolToMolFile(mol, f'$outdir/${tlc}.mol')

# Run rdkit_to_params
p = Params.load_mol(mol, name='${tlc}')
p.convert_mol()
p.dump('$outdir/${tlc}.params')
p.dump_pdb('$outdir/${tlc}_0001.pdb')
print('Generated $outdir/${tlc}.mol, $outdir/${tlc}.params, and $outdir/${tlc}_0001.pdb')
EOF
