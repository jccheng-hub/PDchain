#!/usr/bin/env bash
usage () { cat << EOF
Usage: gen_natbiasjson PDB
Generate a biasjson for ligandmpnn that applies a global weight for native
residues

Parameters:
    PDB                     Input PDB file

Options:
  --global_weight [float]   Apply specified weight for all residues
  --change_weight [str]     Apply specified weight to select residues
                            Overwrites --global_weight
                            Format: CHRESI:WEIGHT ...
                            E.g. --change_weight A17:10.0 A27:10.0
  --change_resn [str]       Changes native residue to specified residue
                            Format: CHRESI:OLC ...
                            E.g. --change_resn A17:W A30:Y
  --focus_set [str]         Only apply bias towards this focused set of residues
                            Format: CHRESI ...
                            E.g. --focus_set A17 A30 A79
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

val_opts=(global_weight change_weight change_resn focus_set)
global_weight="1.0"
change_weight=""
change_resn=""
focus_set=""

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

# Useful functions
calc    () { awk "BEGIN{print($@)}" ; }
randstr () { tr -cd '0-9a-zA-Z' < /dev/urandom | head -c8 | xargs ;}

inpdb=${args[1]}

# Sed script for converting TLC to OLC
olc_sedf=$(echo "
    s/ALA/:A:/g s/ARG/:R:/g s/ASN/:N:/g s/ASP/:D:/g s/CYS/:C:/g
    s/GLU/:E:/g s/GLN/:Q:/g s/GLY/:G:/g s/HIS/:H:/g s/ILE/:I:/g
    s/LEU/:L:/g s/LYS/:K:/g s/MET/:M:/g s/PHE/:F:/g s/PRO/:P:/g
    s/SER/:S:/g s/THR/:T:/g s/TRP/:W:/g s/TYR/:Y:/g s/VAL/:V:/g
    s/[A-Z][A-Z][A-Z]/:X:/g s/:\([A-Z]\):/\1/g
" | sed -e '1d' -e 's/^ *//g' -e '$d' -e 's/ /\n/g')

# Make initial biasjson
biasjson=$(awk -v wt=$global_weight '
    BEGIN { print "{" }
    $1~/^(ATOM|HETATM)$/ && $3=="CA" {
        ch=$5 ; resi=$6 ; tlc=$4
        print "\"" ch resi "\": {\"" tlc "\": " wt "},"
}' $inpdb | sed '$s/,/\n}/' | sed -f <(echo "$olc_sedf"))

# Trim out residues not within focus set
if [[ ! -z "$focus_set" ]] ; then
    focus_set=$(gawk '
        function expand_range(res,    a, b, i, j, n, m, ch, range) {
            n = split(res, a)
            for (i=1; i<=n; i++) {
                ch = a[i]; gsub(/[^A-Z]/, "", ch)
                gsub(/[A-Z]/, "", a[i])
                m = split(a[i], b, "-")
                for (j=b[1]; j<=b[m]; j++) {
                    range = range == "" ? ch j : range " " ch j
                }
            }
            return range
        } { print expand_range($0) }
    ' <<< "$focus_set")
    biasjson=$( 
        awk -F'"' '$2~/^('${focus_set// /|}')$/' <<< $biasjson |
        sed -e '1s/^/{\n/' -e '$s/,$//' -e '$s/$/\n}/'
    )
fi

# Change weight
if [[ ! -z "$change_weight" ]] ; then
    chwt_sedf=""
    for chwt in $change_weight ; do
        chresi=${chwt%:*} ; wt=${chwt#*:}
        chwt_sedf=$(echo "$chwt_sedf" ; echo "/^\"$chresi\"/s/ [^ ]*}(,{0,1})$/ $wt}\1/")
    done
    biasjson=$(sed -E -f <(echo "$chwt_sedf") <(echo "$biasjson"))
fi

# Change resn
if [[ ! -z "$change_resn" ]] ; then
    chrn_sedf=""
    for chrn in $change_resn ; do
        chresi=${chrn%:*} ; rn=${chrn#*:}
        chrn_sedf=$(echo "$chrn_sedf" ; echo "/^\"$chresi\"/s/ {\"[A-Z]\"/ {\"$rn\"/")
    done
    biasjson=$(sed -f <(echo "$chrn_sedf") <(echo "$biasjson"))
fi

echo "$biasjson"
