#!/usr/bin/env bash
usage () { cat << EOF
Usage: calc_od600 ODIN
Calculate OD600 for the specified number of hours.
Uses the logistic growth equation.

To get the initial density from colony count, use the equation:
ODIN = 1 / (8 * VOL)
where ODIN is the initial density, and VOL is the volume of media

The equation above assumes the following values:
1 colony = 10^8 cells
1.0 OD600 = 8 x 10^8 cells

Parameters:
    ODIN                    Initial OD600

Options:
  --capacity [float]        Maximum achievable density
  --doubling_time [float]   Time in hours for E. coli to double
  --rate [float]            Specific growth rate
  --lag_time [int]          Number of hours in lag phase
  --tmax [float]            Maximum time to report (in hours)
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

val_opts=(capacity doubling_time lag_time tmax)
bool_opts=(help)
capacity="3.0"
doubling_time="0.5"
lag_time="1"
tmax="4"

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

odin="${args[1]}"

awk -v odin="$odin" -v capacity="$capacity" \
    -v doubling_time="$doubling_time" -v lag_time="$lag_time" -v tmax="$tmax" \
'
    function get_od(t) {
        if (t <= lag_time) return odin
        od_numer = capacity
        od_denom = 1 + ( (capacity - odin)/odin ) * exp( (-log(2)/doubling_time) * (t-lag_time))
        return od_numer / od_denom
    }

    BEGIN {
        printf "%-8s%s\n", "Time", "OD600"
        for (i=0; i<=tmax; i+=0.25) {
            if ( (capacity - get_od(i)) < 0.001 ) continue
            printf "%-8s%s\n", i, get_od(i)
        }
    }
'
