#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

usage () { cat << EOF
Usage: hqa_convert PDB
Convert OAP to HQA and vice versa

Parameters:
    PDB                     Input PDB

Options:
  --outdir [str]            Output directory
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

val_opts=(outdir)
bool_opts=(help)

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
# }}}

echo "Executing command:" "$0" "$@"
inpdbs=${args[@]:1}

tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX)
trap 'rm -r $tmpdir' EXIT

for inpdb in ${inpdbs[@]} ; do
    if grep "OAP" $inpdb &>/dev/null ; then
        tmppdb="$tmpdir/tmp.pdb"
        awk '
        BEGIN {
            at["C1"]="CZ1" ; at["C5"]="CZ2"
            at["C2"]="CT1" ; at["C6"]="CE1"
            at["C3"]="CI1" ; at["O1"]="OI2"
            at["C4"]="CT2" ; at["N1"]="NE2"
        }
        NR==FNR && $1=="fixedres" { ch=substr($2,1,1); hqa=substr($2,2) } NR==FNR { next }
        $NF=="H" || $4=="OAP" { next }
        $5==ch && $6==hqa { sub("ALA", "HQA", $0) }
        $1=="CoordinateConstraint" && $3=="2E" { $2 = at[$2] ; $3 = hqa ch }
        $1=="AtomPair" && $3=="2E" { $2 = at[$2] ; $3 = hqa ch }
        $1=="AtomPair" && $5=="2E" { $4 = at[$2] ; $5 = hqa ch }
        1
        ' $inpdb $inpdb > $tmppdb
        
        cst="$tmpdir/constraints.cst"
        awk '$0~/###.* constraints .*###/ {x=1-x; next} x' $tmppdb > $cst
        
        mpros $tmppdb --inplace --skip_mpnn --skip_mpnn_score
        mpros $tmppdb --inplace --skip_mpnn --cstfile $cst --ca_stdev 2 --relax --ignore_metals --ligname BCB ZN WTL
        sed -E 's/^ {8}//' << EOF >> $tmppdb
        ### >>> constraints >>> ###
        $(cat $cst)
        ### <<< constraints <<< ###
        $(grep -h "^fixedres" $inpdb)
EOF
    
        # Score
        dss_disicl $tmppdb --append
        calc_scafscore $tmppdb --resn BCB ZN WTL --chri $(sed -n "s/^fixedres //p" $tmppdb)
    
        # Output
        mkdir -p $outdir
        outpdb="$outdir/$(basename ${inpdb%.*})_hqa.pdb"
        cp $tmppdb $outpdb &&
        echo "Saved output to $outpdb"
    
    elif grep "HQA" $inpdb &>/dev/null ; then
        tmppdb="$tmpdir/tmp.pdb"
        awk '
        BEGIN {
            at["CZ1"]="C1" ; at["CZ2"]="C5"
            at["CT1"]="C2" ; at["CE1"]="C6"
            at["CI1"]="C3" ; at["OI2"]="O1"
            at["CT2"]="C4" ; at["NE2"]="N1"
        }
        NR==FNR && $1=="fixedres" { ch=substr($2,1,1); hqa=substr($2,2) } NR==FNR { next }
        $NF=="H" { next }
        $5==ch && $6==hqa {
            if ($3 ~ /^(N|CA|C|O|CB)$/) {
                sub("HETATM", "ATOM  ", $0)
                sub("HQA", "ALA", $0)
            } else if ($3 in at) {
                $0 = substr($0, 1, 13) sprintf("%-4s", at[$3]) substr($0, 18)
                $0 = substr($0, 1, 21) "E" substr($0, 23)
                $0 = substr($0, 1, 22) sprintf("%4s", "2" ) substr($0, 27)
                sub("HQA", "OAP", $0)
                oap[++oap_i] = $0
                next
            } else { next }
        }
        $1=="CoordinateConstraint" && $3 == hqa ch { $2 = at[$2] ; $3 = "2E" }
        $1=="AtomPair" && $3 == hqa ch && $2 !~ /^(N|CA|C|O|CB)$/ {
            $2 = at[$2] ? at[$2] : $2 ; $3 = "2E"
        }
        $1=="AtomPair" && $5 == hqa ch && $4 !~ /^(N|CA|C|O|CB)$/ {
            $4 = at[$2] ? at[$2] : $4 ; $5 = "2E"
        }
        { print }
        END {
            for (i=1; i<=oap_i; i++) { print oap[i] }
        }
        ' $inpdb $inpdb > $tmppdb
        
        cst="$tmpdir/constraints.cst"
        awk '$0~/###.* constraints .*###/ {x=1-x; next} x' $tmppdb > $cst
    
        mpros $tmppdb --inplace --skip_mpnn --cstfile $cst --ca_stdev 2 --relax --ignore_metals --ligname BCB OAP ZN WTL
    
        sed -E 's/^ {4}//' << EOF >> $tmppdb
        ### >>> constraints >>> ###
        $(cat $cst)
        ### <<< constraints <<< ###
        $(grep -h "^fixedres" $inpdb)
EOF

        # Score
        dss_disicl $tmppdb --append
        calc_scafscore $tmppdb --resn BCB OAP ZN WTL --chri $(sed -n "s/^fixedres //p" $tmppdb)
    
        # Output
        mkdir -p $outdir
        outpdb="$outdir/$(basename ${inpdb%.*})_oap.pdb"
        cp $tmppdb $outpdb &&
        echo "Saved output to $outpdb"
    else
        echo "Input PDB $inpdb does not have OAP or HQA"
    fi
done
