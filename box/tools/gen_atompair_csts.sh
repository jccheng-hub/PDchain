#!/usr/bin/env bash
usage () { cat << EOF
Usage: gen_atompair_csts PDB...
Generate AtomPair constraints for Rosetta

Need to specify anchors and targets

Every atom of every target residue will be constrained to every atom of every
anchor residue.

Residues can be both anchors and targets, but they will not be constrained
to itself.

Parameters:
    PDB                     Input PDB

Options:
  --anchors [str]           Anchor residues in chain+resi format (e.g. A1 A2 A3)
  --targets [str]           Target residues in chain+resi format (e.g. A4 A5 A6)
  --ap_stdev [float]        Strength/deviation of AtomPair constraints
                            0.5 is strong, 2.0 is weak
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

val_opts=(anchors targets ap_stdev)
bool_opts=(help)
anchors="REQUIRED"
targets="REQUIRED"
ap_stdev="0.5"

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
    [[ ${!opt} == "REQUIRED" ]] && help="1"
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

inpdbs=${args[@]:1}

awk -v anchors="$anchors" -v targets="$targets" -v ap_stdev="$ap_stdev" '
    BEGIN {
        # Get list of anchors
        n = split(anchors, a, " ")
        for (i=1; i<=n; i++) {
            ch = ri = a[i]
            gsub(/[^A-Z]/, "", ch)
            
            # Expand resi
            m = split(a[i], b, "-")
            gsub(/[^0-9]/, "", b[1])
            b[2] = b[2] ? b[2] : b[1]
            for (ri=b[1]+0 ; ri<=b[2]+0 ; ri++) {
                anch_key[ch ri] = 1
            }
        }

        # Get list of targets
        n = split(targets, a, " ")
        for (i=1; i<=n; i++) {
            ch = ri = a[i]
            gsub(/[^A-Z]/, "", ch)

            # Expand resi
            m = split(a[i], b, "-")
            gsub(/[^0-9]/, "", b[1])
            b[2] = b[2] ? b[2] : b[1]
            for (ri=b[1]+0 ; ri<=b[2]+0 ; ri++) {
                targ_key[ch ri] = 1
            }
        }
    }

    $1~/^(ATOM|HETATM)$/ && $3!="OXT" && $NF!="H" {

        if ($5 $6 in anch_key) {
            anch[$5][$6][$3][1] = $7
            anch[$5][$6][$3][2] = $8
            anch[$5][$6][$3][3] = $9
        }

        if ($5 $6 in targ_key) {
            targ[$5][$6][$3][1] = $7
            targ[$5][$6][$3][2] = $8
            targ[$5][$6][$3][3] = $9
        }
    }

    ENDFILE {
        # Loop through all anchors and targets atoms + calculate distance
        for (ch1 in anch) for (ri1 in anch[ch1]) {
            for (ch2 in targ) for (ri2 in targ[ch2]) {
                if (ch1 == ch2 && ri1 == ri2) continue
                if (ch1 ri1 ":" ch2 ri2 in pairs) continue

                for (at1 in anch[ch1][ri1]) for (at2 in targ[ch2][ri2]) {
                    d = calc_dist(anch[ch1][ri1][at1], targ[ch2][ri2][at2])
                    printf "AtomPair %s %s %s %s ", at1, ri1 ch1, at2, ri2 ch2
                    printf "HARMONIC %s %s\n", d, ap_stdev
                }

                # Track pairs
                pairs[ch1 ri1 ":" ch2 ri2] = 1
                pairs[ch2 ri2 ":" ch1 ri1] = 1
            } 
        }

        delete anch ; delete targ ; delete pairs
    }

    function calc_dist (v1, v2,    dx, dy, dz) {
        dx = v1[1] - v2[1]
        dy = v1[2] - v2[2]
        dz = v1[3] - v2[3]
        return sqrt( dx * dx + dy * dy + dz * dz )
    }
' ${inpdbs[@]}
