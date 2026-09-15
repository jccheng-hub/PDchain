#!/usr/bin/env bash
usage () { cat << EOF
Usage: update_contigs --oldcontigs OLDCONTIGS --newcontigs NEWCONTIGS
Update contigs string

Options:
  --oldcontigs [str]        Old contigs string
  --newcontigs [str]        New contigs string
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

val_opts=(oldcontigs newcontigs)
bool_opts=(help)
oldcontigs="REQUIRED"
newcontigs="REQUIRED"

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
[[ $((${#args[@]}-1)) -lt 0 ]] && help="1"
[[ $help == 1 ]] && { usage ; default_vals ; exit ;}
# <<< Defaults <<< }}}

get_newres () {
    # Take input contigs and oldres list and return new resi
    # Output format: oldres:newres

    # If contigs is NONE, then execute this block (requires inpdb defined)
    [[ $1 == NONE ]] && {
        awk -v oldres="${*:2}" '
            BEGIN {
                n = split(oldres, a, " ")
                for (i=1; i<=n; i++) {
                    split(a[i], b, "-")     ; ch_o = substr(b[1], 1, 1)
                    sub(/^[A-Z]/, "", b[1]) ; if (b[2]=="") b[2]=b[1]
                    for (j=b[1]+0; j<=b[2]+0; j++) res_o[ch_o j] = 1
                }
            }

            $1~/^(ATOM|HETATM)$/ && $3=="CA" && $NF!~/CA/ {
                if (oldres && !($5 $6 in res_o)) next
                print $5 $6 ":" $5 $6
            }
        ' $inpdb
        return
    }

    # For all other contigs, run the following
    local contigs=$(
        sed 's|/0 | |' <<< $1 | sed 's/ /\n/g' |
        sort | awk '$0!~/^[ ]*$/ {a[++n]=$0} END {
            for (i=1; i<n; i++) {
                print a[i]"/0 "
            }
            print a[n]
        }' | paste -sd ' '
    )

    sed -e 's|/|\n|g' -e 's/ /\n/g' <<< "$contigs" |
    awk -v oldres="${*:2}" -F '-' '
        BEGIN {
            n = split(oldres, a, " ")
            for (i=1; i<=n; i++) {
                split(a[i], b, "-")     ; ch_o = substr(b[1], 1, 1)
                sub(/^[A-Z]/, "", b[1]) ; if (b[2]=="") b[2]=b[1]
                for (j=b[1]+0; j<=b[2]+0; j++) res_o[ch_o j] = 1
            }
            s = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            split(s, ch, "") ; ch_i = 1 ; resi = 1
        }
        $1 == 0       { ch_i++ ; next }
        $1 ~ /^[0-9]/ { resi += $1 ; next }
        $1 ~ /^[A-Z]/ {
            ch_o = substr($1, 1, 1)  ;  sub(/^[A-Z]/, "", $1)
            for (i=$1+0; i<=$2+0; i++) {
                if ( oldres && !(ch_o i in res_o) ) { resi++ ; continue }
                print ch_o i ":" ch[ch_i] resi++
            }
        }
    '
}

update_contigs () {
    local oldcontigs="$1" newcontigs="$2" oldtab newtab fixbbres num
    oldtab=$(get_newres "$oldcontigs" | sed 's/:/ /')
    newtab=$(get_newres "$newcontigs" | sed 's/:/ /')
    fixbbres=$(
        echo $oldcontigs | sed -e 's|/|\n|g' -e 's| |\n|g' | sed -n '/^[A-Z]/p' | paste -sd ' '
    )

    # Count total number of residues in newcontigs
    num=$(sed 's|/| |g' <<< "$newcontigs" | awk '{
        for (i=1; i<=NF; i++) {
            if ($i ~ /^[A-Z]/) {
                ori=$i; sub(/-.*$/,"",ori); sub(/^[A-Z]/,"",ori)
                fin=$i; sub(/^.*-/,"",fin)
                num += fin - ori + 1
            } else {
                ori=$i; sub(/-.*$/,"",ori)
                num += ori
            }
        }
    } END { print num }')

    awk -v fix="$fixbbres" -v num="$num" '
        NR == FNR { a[$1]=$2; next }
        NR > FNR { for (i in a) if (a[i]==$1) b[i]=$2 }
        END {
            n = split(fix, afix)
            idx = 0
            for (i=1; i<=n; i++) {
                ch = substr(afix[i], 1, 1)
                pair = afix[i] ; sub(/-/, " " ch, pair) ; split(pair, apair)
                ori = substr(b[apair[1]], 2)
                fin = substr(b[apair[2]], 2)

                ch2_f = substr(b[apair[1]], 1, 1)
                ch2_i = (ch2_i == "") ? ch2_f : ch2_i

                printf "/%s-%s", ori-idx-1, ori-idx-1
                if (ch2_f == ch2_i) {
                    printf "/%s", afix[i]
                } else {
                    printf "/0 %s", afix[i]
                }

                idx=fin
            }
            printf "/%s-%s\n", num-idx, num-idx
        }
    ' <(echo "$oldtab") <(echo "$newtab") |
    sed -e 's|^/||' -e 's|^0-0/||' -e 's|/0-0$||'
}

if [[ $newcontigs == NONE ]] ; then
    echo $oldcontigs
else
    update_contigs "$oldcontigs" "$newcontigs"
fi
