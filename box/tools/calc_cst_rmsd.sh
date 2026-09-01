#!/usr/bin/env bash
usage () { cat << EOF
Usage: calc_cst_rmsd PDB...
Calculates RMSD between input PDB and its constraints.
Constraints in input PDB are specified with...

### >>> constraints >>> ###
AtomPair ...
CoordinateConstraint ...
### <<< constraints <<< ###

Currently only supports AtomPair and CoordinateConstraint

Parameters:
    PDB                     Input PDB file

Options:
  --append                  Append info to input PDB
                            E.g. cst_rmsd 0.2793
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
bool_opts=(append help)

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

inpdbs=${args[@]:1}

[[ $append == 1 ]] && inplace_flag="-i inplace"

awk $inplace_flag '
    $1 !~ /^cst_rmsd$/ { lin[++maxlin] = $0 }

    $1 ~ /^(ATOM|HETATM)$/ {
        crd[$6 $5][$3][1] = $7
        crd[$6 $5][$3][2] = $8
        crd[$6 $5][$3][3] = $9
    }
    
    $0 ~ /^### >>> constraints >>> ###/ { cstblock = 1 ; next }
    $0 ~ /^### <<< constraints <<< ###/ { cstblock = 0 ; next }
    
    cstblock && $1~/^AtomPair$/ {
        act = calc_dist(crd[$3][$2], crd[$5][$4])
        ide = $7
        sum += (ide - act) * (ide - act)
        n++
    }
    
    cstblock && $1~/^CoordinateConstraint$/ {
        arr[1] = $6 ; arr[2] = $7 ; arr[3] = $8
        sum += calc_dist_sq(arr, crd[$3][$2])
        n++
    }
    
    ENDFILE {
        if (inplace::enable) {
            for (i=1; i<=maxlin; i++) print lin[i]
        }
        print "cst_rmsd " sqrt(sum / n)

        delete lin ; maxlin = 0
        sum = 0 ; n = 0
    }
    
    function calc_dist(v1, v2,    dx, dy, dz) {
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        dx = v1[1]-v2[1] ; dy = v1[2]-v2[2] ; dz = v1[3]-v2[3]
        return sqrt( dx*dx + dy*dy + dz*dz )
    }

    function calc_dist_sq(v1, v2,    dx, dy, dz) {
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        dx = v1[1]-v2[1] ; dy = v1[2]-v2[2] ; dz = v1[3]-v2[3]
        return dx*dx + dy*dy + dz*dz
    }
' ${inpdbs}
