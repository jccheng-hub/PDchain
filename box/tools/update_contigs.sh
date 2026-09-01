#!/usr/bin/env bash
usage () { cat << EOF
Usage: update_contigs OLDCONTIGS NEWCONTIGS
Update contigs string

Parameters:
    OLDCONTIGS              Original contigs
    NEWCONTIGS              New contigs

Options:
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

val_opts=()
bool_opts=(help)

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
[[ $((${#args[@]}-1)) -lt 2 ]] && help="1"
[[ $help == 1 ]] && { usage ; default_vals ; exit ;}
# <<< Defaults <<< }}}

get_newres () {
    local contigs=$1
    local oldres=$(echo $@ | awk '{for(i=2;i<=NF;i++)print$i}' | paste -sd ' ')
    [[ -z $oldres ]] && oldres=$(echo $contigs | sed 's|/|\n|g' | sed -n '/^[A-Z]/p' | paste -sd ' ')
    oldres=$(echo $oldres | awk '{
        for (i=1; i<=NF; i++) {
            ch = substr($i, 1, 1); sub(/^[A-Z]/, "", $i)
            delete a; split($i, a, "-") ; if (!(a[2])) a[2]=a[1]
            for (j=a[1]; j<=a[2]; j++) print ch j
        }
    }' | paste -sd ' ')
    echo $contigs | sed 's|/|\n|g' | awk '
    !/^[A-Z]/ {
        delete a ; split($0, a, "-")
        for (i=1; i<=a[1]; i++) { tot++ }
    }
    /^[A-Z]/ {
        ch = substr($0, 1, 1) ; sub(/^[A-Z]/, "", $0)
        delete a ; split($0, a, "-")
        for (i=a[1]; i<=a[2]; i++) {
            tot++
            if (match(" '"$oldres"' ", " "ch i" ") != 0) {
                print ch i ":A" tot
            }
        }
    }'
}

update_contigs () {
    local oldcontigs=$1 newcontigs=$2 oldtab newtab fixbbres num
    oldtab=$(get_newres $oldcontigs | sed 's/:/ /')
    newtab=$(get_newres $newcontigs | sed 's/:/ /')
    fixbbres=$(echo $oldcontigs | sed 's|/|\n|g' | sed -n '/^[A-Z]/p' | paste -sd ' ')
    num=$(awk -F '/' '{
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
    } END { print num }' <<< "$newcontigs")
    awk -v fix="$fixbbres" -v num="$num" '
        NR == FNR { a[$1]=$2; next }
        NR > FNR { for (i in a) if (a[i]==$1) b[i]=$2 }
        END {
            n = split(fix, afix)
            idx = 0
            for (i=1; i<=n; i++) {
                ch=substr(afix[i],1,1); pair=afix[i]; sub(/-/," "ch,pair); split(pair, apair)
                ori=b[apair[1]] ; sub(/[A-Z]/,"",ori)
                fin=b[apair[2]] ; sub(/[A-Z]/,"",fin)
                printf "/%s-%s", ori-idx-1, ori-idx-1
                printf "/%s", afix[i]
                idx=fin
            }
            printf "/%s-%s\n", num-idx, num-idx
        }
    ' <(echo "$oldtab") <(echo "$newtab") | sed 's|^/||'
}

oldcontigs="${args[1]}"
newcontigs="${args[2]}"
if [[ $newcontigs == NONE ]] ; then
    echo $oldcontigs
else
    update_contigs $oldcontigs $newcontigs
fi
