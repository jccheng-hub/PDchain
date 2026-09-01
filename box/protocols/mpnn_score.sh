#!/usr/bin/env bash
usage () { cat << EOF
Usage: mpnn_score PDB...
Run LigandMPNN while forcing 100% sequence recovery to score input sequence.
This "mpnn_score" is the "overall_confidence" metric in the output .fa file.

Because "overall_confidence" can differ depending on seed, this metric doesn't
yield the exact same number for different runs, but they will fall into the same
ballpark.

Parameters:
    PDB                     Input PDB

Options:
  --model_type [str]        The MPNN model used
  --seed [int]              Use a specific seed for running MPNN
                            If not invoked, a random seed will be used instead.
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

val_opts=(model_type seed)
bool_opts=(help)

model_type="protein_mpnn"
seed=""

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

gen_multi_injson () {
    local inpdbs="$@"
    sed 's/ /\n/g' <<< "$inpdbs" |
    jq -R -s 'split("\n") | map(select(length>0)) | map({ (.): "" }) | add'
}

gen_multi_biasjson () {
    local inpdbs="$@"
    awk '
        BEGIN {
            olc["ALA"]="A" ; olc["ARG"]="R" ; olc["ASN"]="N" ; olc["ASP"]="D"
            olc["CYS"]="C" ; olc["GLU"]="E" ; olc["GLN"]="Q" ; olc["GLY"]="G"
            olc["HIS"]="H" ; olc["ILE"]="I" ; olc["LEU"]="L" ; olc["LYS"]="K"
            olc["MET"]="M" ; olc["PHE"]="F" ; olc["PRO"]="P" ; olc["SER"]="S"
            olc["THR"]="T" ; olc["TRP"]="W" ; olc["TYR"]="Y" ; olc["VAL"]="V"
        }

        BEGINFILE { printf FILENAME }
        $1=="ATOM" && $3=="CA" && $4 in olc { printf " " $5 $6 " " olc[$4] " 10.0" }
        ENDFILE { printf "\n" }

        # Output Format:
        # pdb1 chainresi1 aa1 10.0 chainresi2 aa2 10.0 ...
        # pdb2 chainresi1 aa1 10.0 chainresi2 aa2 10.0 ...
        # ...
    ' $inpdbs |
    jq -R -s '
        split("\n") | map(select(length>0)) |
        map(
            split(" ") as $cols | {
                ($cols[0]): reduce range(1 ; $cols|length ; 3) as $i (
                    {} ; . + { ($cols[$i]): { ($cols[$i+1]): $cols[$i+2]|tonumber } }
                )
            }
        ) | add
    '
}

# Assign variables
inpdbs="${args[@]:1}"
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX) ; trap 'rm -r $tmpdir' EXIT

# Make jsons
gen_multi_injson $inpdbs > $tmpdir/inputs.json
gen_multi_biasjson $inpdbs > $tmpdir/bias.json

# Report the model and seed used
echo -n "Running LigandMPNN with --model_type $model_type with "
if [[ -z $seed ]] ; then
    echo "a random seed for scoring."
else
    echo "--seed $seed for scoring."
    mpnn_opts+=("--seed $seed")
fi

# Run MPNN
eval "
    $pixirun ligandmpnn ligandmpnn \
        --pdb_path_multi '$tmpdir/inputs.json' \
        --bias_AA_per_residue_multi '$tmpdir/bias.json' \
        --model_type '$model_type' \
        --out_folder '$tmpdir' ${mpnn_opts[*]}
"

# Get scores
for inpdb in $inpdbs ; do
    bn=$(basename ${inpdb%.*})
    mpnn_score=$(awk '$2 == "id=1," { gsub(",", "", $0)
        for(i=1;i<=NF;i++) {
            if ($i~/^overall_confidence=/) {
                split($i, a, "=") ; print a[2] ; exit
            }
        }
    }' $tmpdir/seqs/${bn}.fa)
    awk -i inplace -v mpnn_score="$mpnn_score" -v model_type="$model_type" '
        $1 != model_type "_score" { print } ENDFILE {
            printf "%s_score %s\n", model_type, mpnn_score
        }
    ' $inpdb
done

# Report
grep -H "^${model_type}_score " $inpdbs
