#!/usr/bin/env bash
usage () { cat << EOF
Usage: filter_pool INDIR METRIC...
Filter designs in pool using specified metrics.

Approach is basically a venn diagram thing.
1.  Make an ordered list of each metric
2.  Designs that appear in every list "passes"
3.  If the number of passed designs is over the specified limit, shorten each
    list by 1 design (i.e. take the worst one out from each list), then
    repeat.

Parameters:
    INDIR                   Input directory
    METRIC                  Metric string
                            Needs to have the format of min:metric, max:metric,
                            or val:metric=ideal_val
                            E.g. min:fin_total_energy_per_res
                            E.g. max:protein_mpnn_confidence
                            E.g. val:dsasa=0.75

Options:
  --perc [float]            Target top percent
  --num [int]               Target top number (overrides --perc)
  --outdir [str]            Path to output directory to store passed designs
  --outfile [str]           Path to output text file listing passed designs
                            If invoked, will not copy passed designs into output
                            directory (i.e. --outdir is ignored)
  --inplace                 Will archive away failed designs in input directory
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

val_opts=(num perc outfile outdir)
bool_opts=(inplace help)
num=""
perc="10"
outfile=""
outdir="filtered"

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

# Input arguments
indir=${args[1]%/}
mets=(${args[@]:2})

# Output list of available metrics
[[ -z ${mets[*]} ]] && {
    gawk '
        $1~/^[a-z]/ {
            met[$1]++
            max = met[$1]>max ? max=met[$1] : max
            sum[$1]+=$2
            val[$1][met[$1]] = $2
        }

        END {
            print "List of metrics + averages + sd:"
            for (i in met) if (met[i]==max) {
                printf "%-28s%-12s%s\n", i, sum[i] / max, get_sd(val[i])
            }
        }

        function get_sd(v,    key, mean, n, std) {
            for (key in v) { mean += v[key] ; n++ }
            mean = mean / n
            for (key in v) {
                std += (v[key] - mean) * (v[key] - mean)
            }
            return sqrt(std / n)
        }
    ' $indir/*.pdb ; exit
}

# Check for metrics formatting
fmt_check=$(sed 's/ /\n/g' <<< "${mets[*]}" | gawk -F ':' '
    $1~/^(min|max|val)$/{i++} END {if (i==NR) {print 1} else {print 0}}'
)

# Report inputs
if [[ $fmt_check == 1 ]] ; then
    echo "Using the following metrics for $indir/*.pdb"
    sed 's/ /\n/g' <<< "${mets[*]}" | gawk -F ':' '{printf "%-32s%s\n", $2, $1}'
else
    echo "All metrics provided must follow the format of max:metric, min:metric, or val:metric=ideal_val"
    exit
fi

# Initialize
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX) ; trap 'rm -r $tmpdir' EXIT

# Get metrics
for met in ${mets[@]} ; do
    met=($(echo $met | sed 's/[:|=]/ /g'))
    if [[ ${met[0]} =~ ^min$ ]] ; then
        echo "Generating ascending sorted list for ${met[1]}..."
        grep -H "^${met[1]} " $indir/*.pdb 2>/dev/null | sort -k2,2n -k1R |
        sed 's/:/ /' > $tmpdir/${met[1]}.txt
    elif [[ ${met[0]} =~ ^max$ ]] ; then
        echo "Generating descending sorted list for ${met[1]}..."
        grep -H "^${met[1]} " $indir/*.pdb 2>/dev/null | sort -k2,2rn -k1R |
        sed 's/:/ /' > $tmpdir/${met[1]}.txt
    elif [[ ${met[0]} =~ ^val$ ]] ; then
        [[ -z ${met[2]} ]] && {
            echo "Ideal value for ${met[1]} not provided. Omitting this metric!" ; continue
        }
        echo "Generating ascending sorted list for abs(${met[1]} - ${met[2]}) ..."
        grep -H "^${met[1]} " $indir/*.pdb 2>/dev/null | awk -v idval=${met[2]} '
            function abs(x) { return x < 0 ? -x : x }
            { print $1 "_dev_from_" idval, abs(idval-$2) }
        ' | sort -k2,2n -k1R | sed 's/:/ /' > ${tmpdir}/${met[1]}.txt
    else
        echo "Error when parsing val: ${met[1]}."
    fi

    # Trim out files
    [[ -s $tmpdir/${met[1]}.txt ]] || {
        echo "Metric ${met[1]} not found in input PDBs... Omitting this metric!"
        rm $tmpdir/${met[1]}.txt
    }

    # Ignore metric if best and worst values are identical
    awk 'NR==1 {best=$NF} {worst=$NF} END { if (NR>10 && best==worst) exit 1 }' $tmpdir/${met[1]}.txt || {
        echo "The best and worst values for ${met[1]} are identical... Omitting this metric!"
        rm $tmpdir/${met[1]}.txt
    }

    unset met
done

# Get threshold
totnum=$(echo $indir/*.pdb | wc -w)
topnum=$(awk 'BEGIN{printf "%.0f", ('"$perc * $totnum * 0.01"')}')
[[ -n $num ]] && topnum=$num
echo "Threshold: $topnum / $totnum"

gawk -v topnum=$topnum -v ranked="$tmpdir/ranked.txt" '
{ dt[FILENAME][FNR]["pdb"] = $1 ; dt[FILENAME][FNR]["val"] = $3 } END {
    lim = FNR ; metcount = length(dt)
    while (lim > 0) {
        for (txt in dt) {
            pdb = dt[txt][lim]["pdb"]
            if ( !(pdb in elimkey) ) {
                elimkey[pdb] = 1 ; elim[++f] = pdb
            }
        }

        # Record max line
        if (!maxlim && topnum > FNR - f) maxlim = lim

        # Break loop if all designs were eliminated
        if (f == FNR) { break } else { lim-- }
    }
    print "Max percentile: " 100 * maxlim / FNR

    # Reverse elimination list to give ranked list
    for (i=f; i>=1; i--) print elim[i] >> ranked
}
' $tmpdir/*.txt

# Output
if [[ -n $outfile ]] ; then
    mkdir -p $(dirname $outfile)
    cp $tmpdir/ranked.txt $outfile
    echo "Ranked list (best to worst) saved to $outfile"
elif [[ $inplace == 1 ]] ; then
    tormv=($(gawk -v topnum=$topnum 'FNR > topnum' $tmpdir/ranked.txt))
    if [[ -n $tormv ]] ; then
        mkdir -p $indir/archive
        if [[ $(echo $indir/*.pdb | wc -w) -gt $topnum ]] ; then
            mv ${tormv[@]} $indir/archive
            echo "Archived away ${#tormv[@]} failed designs"
        else
            echo "Number of designs in $indir is below threshold of $topnum"
            echo "No designs were archived"
        fi
    else
        echo "No designs were archived"
    fi
    cp $tmpdir/ranked.txt $indir
    echo "Ranked list (best to worst) saved to $indir/ranked.txt"
    
    echo "${mets[@]}" | sed 's/ /\n/g' > $indir/ranked_metrics.txt
else
    mkdir -p $outdir
    cp $tmpdir/ranked.txt $outdir
    cp $(head -n $topnum $tmpdir/ranked.txt) $outdir
    echo "Transferred top $topnum designs to $outdir"
    echo "Ranked list (best to worst) saved to $outdir/ranked.txt"

    echo "${mets[@]}" | sed 's/ /\n/g' > $outdir/ranked_metrics.txt
fi
