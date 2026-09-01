#!/usr/bin/env bash
usage () { cat << EOF
Usage: clean_pdb PDB [PDB...]
Cleans PDB

Parameters:
    PDB                     Input pdb

Options:
  --keep_resn [str]         Residue names to keep outside of 20 AAs
                            E.g. --keep_resn ZN BMH
  --chains [str]            List of chains to preserve. Default is all chains.
                            E.g. --chains A B C
  --outdir [str]            Output directory
  --suffix [str]            Suffix
  --renumber                Renumber each chain starting at 1
  --help                    Display this help and exit
EOF
}

# >>> Defaults >>> {{{
args=($(echo $0 $@ | sed "s/--.*//"))
opts=($(echo $0 $@ | sed "s/--/\n--/" | sed "1d"))
optchk () { [[ " ${opts[*]} " =~ " $1 " ]] && echo 1 || echo 0 ;}
optarg () {
    sed "s/--/\n--/g" <<< "${opts[*]}" |
    sed -n "1,/^$1 /s/^$1 //p" | sed 's/ \+$//'
}

val_opts=(keep_resn chains outdir suffix)
bool_opts=(renumber help)
pdb_list=(${args[@]:1})
keep_resn=""
outdir="."
suffix="_clean"

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

inpdbs="${args[@]:1}"
mkdir -p $outdir

# >>> clean_pdb () >>> {{{
clean_pdb () {
    local pdb=$1
    awk '
    "'"$chains"'" != "" && $5 !~ /^('"${chains// /|}"')$/ {next}

    $1 ~ /^(ATOM|HETATM)$/ {
        alt=substr($0,17,1)
        if (!(alt == "A" || alt == " ")) next
        $0 = substr($0,1,16) " " substr($0,18)  # Gets rid of alt conformers
        $0 = substr($0,1,72) " " substr($0,74)  # Gets rid of segi
        lin[++i] = $0
        if ($3 == "CA") pass[$5$6] = 1
        if ("'"$keep_resn"'" != "" && $4~/^('${keep_resn// /|}')$/) pass[$5$6] = 1
    }

    END {
        for (j=1;j<=i;j++) {
            delete a ; split(lin[j], a)
            if (pass[a[5]a[6]] == 1) {
                if ('$renumber') {
                    if (ch=="" || ch!=a[5]) { d = a[6] - 1 ; x = 0 ; ch = a[5] }

                    # Check for insertion numbers
                    if (a[6]~/[A-Z]/) {
                        if (x_id!=a[6]) { x_id=a[6]; x++ }
                    }
                    sub(/[A-Z]/,"",a[6])
                    resi = sprintf("%4d", a[6] - d + x)
                    lin[j] = substr(lin[j],1,22) resi " " substr(lin[j],28)
                }
                print lin[j]
            }
        }
    }
    ' $pdb
}
# <<< clean_pdbs () <<< }}}

for inpdb in $inpdbs ; do
    cleaned="${outdir%/}/$(basename ${inpdb%.*})${suffix}.pdb"
    clean_pdb $inpdb > $cleaned
    [[ -s $cleaned ]] && echo "Cleaned: $cleaned"
done
