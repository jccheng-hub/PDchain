#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

usage () { cat << EOF
Usage: evo_rfd_chain INDIR METRIC...
Evolve rfd_chain designs found in INDIR

Parameters:
    INDIR                   Input directory with rfd_chain designs
    METRIC                  Metric string
                            Needs to have the format of min:metric, max:metric,
                            or val:metric=ideal_val
                            E.g. min:fin_reu_per_res
                            E.g. max:protein_mpnn_confidence
                            E.g. val:dsasa=0.7

Options:
  --outdir [str]            Output directory for design pool
  --poolsize [int/str]      Determines the design pool size. The design pool is
                            where designs will be aggregated and continuously
                            evolved. The design pool also determines how designs
                            are ranked, as performance on each metric is
                            relative. Use keyword "ALL" to set the poolsize to
                            the number of designs found in the input directory.
  --poolperc [float]        Percent of poolsize to use as the subset pool. The
                            subset pool is randomly chosen from the design pool
                            in each design run, and the top ranking design from
                            the subset pool is chosen for redesign.
                            For example, --poolsize 200 --poolperc 5 would
                            result in a subset poolsize of 10, meaning that 10
                            designs are randomly chosen from the pool of 200,
                            and the top design of that 10-member pool is chosen
                            for redesign.
                            Increase this value to increase the odds of higher
                            ranked designs being chosen (at the cost of
                            diversity).
                            Accepted range: 0-100
                            Note: a subset poolsize of less than 1 will be
                            treated as 1, and a subset poolsize of 1 would mean
                            that the design chosen for redesign is truly random.
  --famiperc [float]        Percent of poolsize a single family's maximum
                            allowed occupancy. For example, a famiperc of 20 for
                            a poolsize of 200 would mean that a maximum of 40
                            family members are allowed in the pool. If the
                            number of family members exceed 40, then the worst
                            performing one is automatically archived.
  --gen0perc [float]        If the percent of gen0 designs in pool is below this
                            threshold, then quit. If unset, then the continuous
                            evolution will persist until the time limit.
  --tlim [str]              Time limit. Use format of [float][d/h/m/s].
                            E.g. --tlim 8.5h
  --auto_update             Automatically updates output directory with input
                            directory every round. Useful when the input
                            directory is dynamically being updated with designs.
  --skip_diffusion          Skip RFdiffusion step during redesign.
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

val_opts=(
    outdir              poolsize            clash_cut
    tlim                ss_trim             vary_linkers        natbias
    backrub             timesteps           threads             temperature
    relax_repeats       redes_dist          model_type          design_cycles
    poolperc            gen0perc            nterm_trim
    cterm_trim          ca_stdev            nterm_add           cterm_add
    max_totres          helix_cap           reset_perc          loop_cut
    rog_cut             famiperc            select_met
)
bool_opts=(
    active_site         disallow_cys        dec_only            fix_bb
    fix_chi             idealize            inc_only            monomer_ROG
    partial             ppi_mode            relax               sc_context
    inpaint_seq         help                skip_diffusion      ignore_metals
    auto_update
)

outdir="./evolved"
poolsize="ALL"
poolperc="0"
famiperc="100"
clash_cut=""
loop_cut=""
rog_cut=""
tlim="2h"
design_cycles="1"
ss_trim="4"
nterm_trim="1"
cterm_trim="1"
nterm_add="0"
cterm_add="0"
vary_linkers="1"
natbias=""
backrub=""
timesteps="15"
threads="1"
temperature="0.1"
relax_repeats="1"
redes_dist=""
gen0perc=""
model_type="protein_mpnn"
ca_stdev=""
max_totres=""
helix_cap=""
reset_perc="0"
select_met=""

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

# Initialize
shopt -s nullglob
echo "Executing command:" "$0" "$@"
indir="${args[1]%/}"
mets="${args[@]:2}"
outdir="${outdir%/}" ; mkdir -p $outdir
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX) ; trap 'rm -r $tmpdir' EXIT

# Add max_totres
[[ -n $max_totres ]] && {
    all_opts+=("--mintot $max_totres" "--addtot 0")
}

# Check time limit
tlim_s=$(gawk -v tlim="$tlim" 'BEGIN {
    mult["s"] = 1 ; mult["m"] = 60 ; mult["h"] = 3600 ; mult["d"] = 86400
    if ( match(tlim, /^([0-9.]+)([a-z])$/, a) && a[2] in mult ) {
        print int(a[1] * mult[a[2]])
    } else {
        print "ERROR"
    }
}')

# Handle time limit error
if [[ $tlim_s != "ERROR" ]] ; then
    echo "Time limit set to $tlim ($tlim_s seconds)"
else
    echo "Error with parsing time. Make sure --tlim follows the format of [float][s/m/h/d]."
    exit
fi

# Get final time
time_f=$(echo "$SECONDS + $tlim_s" | bc)

# Set contigs to NONE with --skip_diffusion
[[ $skip_diffusion == 1 ]] && contigs="NONE"

# Determine pool size and famisize
[[ $poolsize == "ALL" ]] && poolsize=$(echo $indir/*pdb | wc -w)
famisize=$(echo "$poolsize * $famiperc * 0.01" | bc | sed 's/^0$/1/')
echo "Pool size: $poolsize"
echo "Family size: $famisize"

# Define function for removing worst family members
trim_family () {
    local rankfile="$1" famisize="$2"
    local tormv=$(awk -v famisize=$famisize '
        {
            f++ ; file[f] = fami[f] = $1
            sub(/_gen[0-9]+_[a-zA-Z]{6}[.]pdb$/, "", fami[f])
            c[fami[f]]++
        }

        c[fami[f]] > famisize { print file[f] }
    ' $rankfile | paste -sd ' ')
    [[ -n $tormv ]] && {
        mv $tormv $outdir/archive &&
        echo "Archived the following designs for exceeding the maximum family size of $famisize: $tormv"
    }
}

# Check metrics
[[ -z $mets ]] && {
    echo "Need to enter metrics for evolution"
    exit
}

# Make initial design pool
echo "Generating initial design pool..."
for inpdb in $indir/*.pdb ; do
    bnpdb=$(basename ${inpdb%.*})
    [[ $(echo $outdir/$bnpdb*.pdb $outdir/archive/$bnpdb*.pdb | wc -w) == 0 ]] && {
        cp -v $inpdb $outdir/${bnpdb}_gen0_orides.pdb
    }
done
filter_pool $outdir $mets --num $poolsize --inplace

# Start evolving
while [[ $SECONDS -lt $time_f ]] ; do

    # Get ranked file
    rankfile="$tmpdir/ranked.txt"
    filter_pool $outdir $mets --outfile $rankfile

    # Pull random subset of designs from pool, then choose best one from subset
    cursize=$(echo $outdir/*.pdb | wc -w)
    subsize=$(echo "$cursize * $poolperc / 100" | bc | sed 's/^0$/1/')
    echo "Subpool size with current poolsize of $cursize: $subsize"

    # Take top 2 from random subset, then take the one with smaller family size
    inpdb=($(
        awk '
            {
                file[++f]=$1
                fami[f]=$1
                sub(/_gen[0-9]+_[a-zA-Z]{6}[.]pdb$/, "", fami[f])
                c[fami[f]]++
            }

            END {
                for (i=1; i<=f; i++) {
                    print file[i], i, c[fami[i]]
                }
            }
        ' $rankfile | shuf -n$subsize | sort -k2n | head -n2 | awk '
            { lin[NR]=$0 ; rank[NR]=$2 ; fami[NR]=$3 }
            END {    
                if (!fami[2] || fami[1] <= fami[2]) {
                    print lin[1]
                } else {
                    print lin[2]
                }
            }
        '
    ))
    echo "Picked $inpdb (rank ${inpdb[1]}) for evolution"

    # Prepare values for rog_cut and loop_cut (tolerate small increase in base value)
    unset extra_opts
    rog_cut_tol=$(awk '
        $1=="rog_ala" && $2 { sum += $2 ; idx++ ; val[idx] = $2 }
        END {
            if (idx + 0 > 20) {
                avg = sum / idx
                for (i=1; i<=idx; i++) { ssq += ( val[i] - avg )^2 }
                print sqrt(ssq / idx)
            } else if (idx + 0 > 0) {
                print 9999
            } else {
                print 0
            }
        }
    ' $outdir/*.pdb)
    loop_cut_tol=$(awk '
        $1=="perc_loop" && $2 { sum += $2 ; idx++ ; val[idx] = $2 }
        END {
            if (idx + 0 > 20) {
                avg = sum / idx
                for (i=1; i<=idx; i++) { ssq += ( val[i] - avg )^2 }
                print sqrt(ssq / idx)
            } else if (idx + 0 > 0) {
                print 9999
            } else {
                print 0
            }
        }
    ' $outdir/*.pdb)
    rog_cut_true=$(awk -v tol="$rog_cut_tol" '$1=="rog_ala" && $2 {print $2 + tol}' $inpdb)
    loop_cut_true=$(awk -v tol="$loop_cut_tol" '$1=="perc_loop" && $2 {print $2 + tol}' $inpdb)
    [[ -n $rog_cut_true ]] && extra_opts+=("--rog_cut $rog_cut_true")
    [[ -n $loop_cut_true ]] && extra_opts+=("--loop_cut $loop_cut_true")
    [[ -n $extra_opts ]] && echo "Appending the following options: ${extra_opts[@]}"

    # Run rfd_chain
    outbn=($(basename ${inpdb%.*} | sed -E 's|_([^_]+)_([^_]+)$| \1 \2|'))
    oldgen=$(tr -cd '[0-9]' <<< "${outbn[1]}")
    newgen=$(echo "$oldgen + 1" | bc)
    rfd_chain $contigs \
        --inpdb $inpdb \
        --outprefix $outdir/${outbn}_gen${newgen} \
        --randsuffix \
        ${all_opts[@]} ${extra_opts[@]}

    # Auto update
    [[ $auto_update == 1 ]] && {
        echo "Updating pool because of --auto_update..."
        for inpdb in $indir/*.pdb ; do
            bnpdb=$(basename ${inpdb%.*})
            [[ $(echo $outdir/$bnpdb*.pdb $outdir/archive/$bnpdb*.pdb | wc -w) == 0 ]] && {
                cp -v $inpdb $outdir/${bnpdb}_gen0_orides.pdb
            }
        done
    }

    # Filter
    filter_pool $outdir $mets --num $poolsize --inplace
    trim_family $outdir/ranked.txt $famisize

    # If gen0perc is set and the gen0 pool is below threshold, then exit
    [[ -n $gen0perc ]] && {
        all_gen0=($(
            echo $outdir/*pdb | sed 's/ /\n/g' | grep -E 'gen0_[^_]{6}[.]pdb$'
        ))
        cur_gen0perc=$(echo "100 * ${#all_gen0[@]} / $poolsize" | bc -l)
        if [[ $(echo "$cur_gen0perc > $gen0perc" | bc) == 1 ]] ; then
            echo "Current gen0perc is $cur_gen0perc%, above the cutoff of $gen0perc."
        else
            echo "Current gen0perc is $cur_gen0perc%, meeting the cutoff of $gen0perc."
            echo "Exiting..."
            exit
        fi
    }
done
