#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
paramsdir="$pixiroot/box/params"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

usage () { cat << EOF
Usage: rfd_chain [CONTIGS]
Run RFdiffusion and LigandMPNN with Rosetta refinement

Parameters:
    CONTIGS                 RFdiffusion contigs.
                            E.g. 180-180
                            E.g. 20-30/A5-76/20-30
                            If CONTIGS is not provided, then one will be
                            generated using the contigs options.
                            To skip RFdiffusion, input NONE as your contigs

General options:
  --inpdb [str]             Input PDB path
  --indir [str]             Input directory with PDBs. A random PDB will be
                            pulled from this directory if specified.
  --ligname [str]           Ligand name (3-letter code)
                            Can specify multiple ligands, but only the first is
                            considered for RFdiffusion's auxiliary potential.
                            All ligands are considered for MPNN and Rosetta.
  --outprefix [str]         Output prefix
                            E.g. --outprefix ./outputs/rfd_des
  --randsuffix              Generate a random string for the suffix
  --help                    Display this help and exit

Contigs options:
  --fixbbres [str]          Residues that will have their backbones fixed during
                            diffusion (chain+residue format).
                            E.g. --fixbbres A7-10 A39 A72
                            If --fixedres is used but --fixbbres is not, then
                            --fixbbres will take --fixedres as input.
  --random_order            If --fixbbres (or --fixedres) is provided, then the
                            order of the fixed residues are shuffled during
                            contigs generation.
  --mingap [int]            Minimum number of diffusable residues in each gap
  --addgap [int]            Added to --mingap determines maximum residue gap
  --mintot [int]            Minimum total residues of the entire protein
  --addtot [int]            Added to --mintot determines maximum total residues
  --vary_linkers [int]      Allow diffusable regions to vary by this number
  --inc_only                If --vary_linkers is provided, then the diffusable
                            regions cannot decrease in length.
  --dec_only                If --vary_linkers is provided, then the diffusable
                            regions cannot increase in length.
  --ss_to_contigs           Converts secondary structure in input PDB to contigs
                            used for RFdiffusion. Secondary structure residues
                            are fixed, and loops are rediffused.
  --ss_trim [int]           If --ss_to_contigs is used, allows the specified
                            number of ss residues neighboring loops to be
                            rediffused by RFdiffusion.
                            E.g. If DSSP string is LLHHHHHLL with --ss_trim 2,
                            then it would be treated as LLLLHLLLL for
                            RFdiffusion.
  --nterm_add [int]         Number of residues to add to N-terminus
  --cterm_add [int]         Number of residues to add to C-terminus
  --nterm_trim [int]        Trim specified number of residues from N terminus
  --cterm_trim [int]        Trim specified number of residues from C terminus
  --regap                   If --ss_to_contigs is used, then the fixed regions
                            will be treated as --fixbbres during contigs
                            generation.

RFdiffusion options:
  --timesteps [int]         Number of timesteps for RFdiffusion
  --rog_cut [float]         Cutoff for poly-ala ROG metric. Lower values are
                            more compact and globular. For a 200-residue
                            protein, --rog_cut 16 is a good start. Adjust value
                            for larger/smaller proteins.
  --loop_cut [float]        Cutoff for percent loop. Lower values have more
                            alpha helices and beta sheets. Range from 0-100.
                            Recommend --loop_cut 30 for most applications.
  --clash_cut [int]         Cutoff for clashes with ligand. If not provided,
                            will use number of non-H ligand atoms divided by 2.
  --backrub [int]           After running RFdiffusion, rerun RFdiffusion again
                            with partial diffusion. This specifies the number
                            of timesteps for this partial diffusion rerun.
  --partial                 Use partial diffusion
  --monomer_ROG             Apply monomer radius of gyration potential
  --ppi_mode                Apply settings recommended for PPI
                            I.e. noise_scale_ca=0.5 and noise_scale_frame=0.5
  --inpaint_seq             If a residue has its backbone fixed for RFdiffusion
                            but is not fixed for sequence design, then hide the
                            residue identity during RFdiffusion.

LigandMPNN + PyRosetta options:
  --design_cycles [int]     Number of MPNN-FastRelax cycles.
                            Every cycle, the new design is accepted or rejected
                            depending on whether they improve the Rosetta score
                            and/or MPNN score.
                            Recommendations:
                            --design_cycles 2 for ok/cheap run
                            --design_cycles 5 for rigorous/expensive run
  --select_met [str]        Selection metric. The specified metric is used to
                            accept/reject new designs in between design cycles.
                            If not specified, then will default to Rosetta and
                            MPNN scores.
                            Specify with the following formats:
                            min:metric, max:metric, or val:metric=ideal_val
                            E.g. --select_met min:ddg
                            E.g. --select_met val:dsasa=0.75
  --model_type [str]        MPNN model
                            E.g. ligand_mpnn, protein_mpnn, soluble_mpnn, etc.
  --temperature [float]     MPNN temperature
  --fixedres [str]          Residues to fix (exclusive with --redesres)
  --redesres [str]          Residues to redesign (exclusive with --fixedres)
  --natbias [float]         Bias design towards each residue at each position
                            found in the input PDB. A value of 1.0 is a good
                            start. Values 3.0 and above tend to just fix the
                            residues entirely.
  --relax_repeats [int]     Number of relax repeats to use
  --cstfile [str]           Path to cst file
  --fix_stdev [float]       Strength of constraints for fixed residues
  --bb_stdev [float]        Strength of constraints for backbone atoms
  --ca_stdev [float]        Strength of constraints for CA atoms
  --lig_stdev [float]       Strength of constraints for ligand atoms
                            Need to invoke --ligname
  --ap_stdev [float]        Strength of constraints between fixed residue atoms
                            and ligand atoms (using AtomPair)
                            Need both --ligname and --fixedres
  --threads [int]           Number of threads to use for PyRosetta
  --nbr_dist [float]        Neighbor distance for LigandMPNN redesign around
                            ligands
  --disallow_cys            Disallow cysteines during design
  --sc_context              Use side chain atoms as ligands during LigandMPNN
  --idealize                Idealize bond lengths and bond angles
  --relax                   Use relax as a refinement
  --fix_bb                  Fix backbones
  --fix_chi                 Fix chi angles
  --skip_mpnn               Skip MPNN
  --skip_refine             Skip refinement with PyRosetta
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
    addgap              addtot              bb_stdev            ca_stdev            
    clash_cut           cstfile             design_cycles       fix_stdev
    fixbbres            fixedres            indir               inpdb
    lig_stdev           ligname             loop_cut            mingap
    mintot              model_type          natbias             numdes
    outprefix           ppi_hotspots        ppi_target
    redesres            relax_repeats       rog_cut             ss_trim
    temperature         threads             timesteps           vary_linkers 
    ap_stdev            backrub             nterm_trim          cterm_trim
    nterm_add           cterm_add           helix_cap           reset_perc
    nbr_dist            select_met
)

bool_opts=(
    disallow_cys        dec_only            ss_to_contigs     
    fix_bb              fix_chi             idealize            inc_only            
    monomer_ROG         partial             ppi_mode            random_order        
    randsuffix          regap               relax               inpaint_seq
    sc_context          skip_mpnn           skip_refine         ignore_metals
    help
)

addgap="100"
addtot="80"
ap_stdev=""
backrub=""
bb_stdev=""
ca_stdev=""
clash_cut=""
cstfile=""
design_cycles="2"
fix_stdev=""
fixbbres=""
fixedres=""
indir=""
inpdb=""
lig_stdev=""
ligname=""
loop_cut="9999"
mingap="20"
mintot="120"
model_type="protein_mpnn"
natbias=""
numdes="1"
outprefix="outputs/rfd_des"
ppi_hotspots=""
ppi_target=""
redesres=""
relax_repeats="1"
rog_cut="9999"
ss_trim="1"
nterm_trim="1"
cterm_trim="1"
nterm_add="0"
cterm_add="0"
temperature="0.1"
threads="1"
timesteps="15"
vary_linkers="0"
helix_cap=""
reset_perc="0"
nbr_dist=""
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
[[ $help == 1 ]] && { usage ; default_vals ; exit ;}
# }}}

shopt -s nullglob
echo "Executing command:" "$0" "$@"

# >>> General functions >>> {{{
expand_res_range () {
    [[ $# == 0 ]] && return
    echo $@ | sed 's/ /\n/g' | awk -F '-' '{
        ch = substr($0,1,1) ; sub(/^[A-Z]/,"",$0) ; $2 = (!$2)? $1 : $2
        for (i=$1+0; i<=$2+0; i++) print ch i
    }' | paste -sd ' '
}

# Convert chainresi format to pymol selection (e.g. A73 to (///A/73))
res_to_sel () {
    sed 's|[A-Z]|/&/|g' <<< $@ |
    awk '{for(i=1;i<=NF-1;i++)printf"(//%s) or ",$i;printf"(//%s)\n",$NF}'
}

# Take input contigs and oldres list and gets new resi | format oldres:newres
get_newres () {

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

calc_perc_loop () {
    local inpdb=($@)
    dss_disicl ${inpdb[0]} --append --fixedres ${inpdb[@]:1}
    echo -n "$(basename ${inpdb%.*}) "
    grep -h "^perc_loop " $inpdb | head -n1
}

collapse_contigs () {
    local contigs="$@"
    echo $contigs | sed 's/ /\n/g' | awk '
        BEGIN { srand() } {
            contigs = ""
            n = split($0, a, "/") 
            for (i=1; i<=n; i++) {
                m = split(a[i], b, "-")
                if (b[1] == 0) {
                    contigs = contigs "/0" ; continue
                }
                if (match(a[i], /^[0-9]/)) {
                    r = get_rand_int(b[1], b[2])
                    contigs = contigs == "" ? r "-" r : contigs "/" r "-" r
                } else {
                    contigs = contigs == "" ? a[i] : contigs "/" a[i]
                }
            }
            print contigs
        }
    
        function get_rand_int(min, max) {
            return int(rand() * (max - min + 1)) + min
        }
    ' | paste -sd ' '
}
# <<< General functions <<< }}}

# Initialize
mkdir -p $(dirname $outprefix)
contigs=${args[@]:1}
outdir=$(realpath $outprefix | xargs dirname)
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX) ; trap 'rm -r $tmpdir' EXIT

# >>> rfd_chain () >>> {{{
rfd_chain () {
    local contigs=$1 outprefix=$2
    local tmpdir=$tmpdir
    local inpdb=$inpdb
    local ppi_target=$ppi_target
    local redesres=$redesres
    local fixedres=$fixedres
    local fixbbres=$fixbbres

    # Generate fixbbres if not provided
    [[ -z $fixbbres && -z $contigs ]] && {
        local fixbbres=$(
            sed 's|/|\n|g' <<< $contigs | grep '^[A-Z]' | paste -sd ' '
        )
    }
    [[ -z $fixbbres && -n $fixedres ]] && local fixbbres="$fixedres"

    # If ppi_target was provided with contigs, then generate contigs
    [[ -z $contigs && -n $ppi_target ]] && {
        echo "PPI target provided without contigs. Generating contigs..."
        local fixedres=$(awk ' 
            $1=="ATOM" && $3=="CA" { chain[++n]=$5; resi[n]=$6 } END {
                ch = chain[1] ; r1 = resi[1]
                for (i=1; i<=n; i++) if (chain[i]==ch) r2 = resi[i]
                print ch r1 "-" r2
            }
        ' $ppi_target)
        local contigs=$(gen_contigs ${all_opts[@]})
        local contigs="$contigs/0 $fixedres"
        local inpdb="$ppi_target"
    }

    # Collapse input contigs if it is provided
    [[ -n $contigs ]] && {
        local contigs=$(collapse_contigs "$contigs")
    }
    
    # Convert redesres to fixedres
    [[ -n $redesres && -n $inpdb ]] && {
        local redesres=$(expand_res_range $redesres)
        local fixedres=$(awk '
            $1=="ATOM" && $3=="CA" && $5$6!~/^('${redesres// /|}')$/ {print$5$6}
        ' $inpdb | paste -sd ' ')
    }
    
    # Get data from previous runs
    [[ -n $inpdb ]] && {
        local oldcontigs=$(sed -n "s/^contigs //p" $inpdb)
        local oldfixed=$(sed -n "s/^fixedres //p" $inpdb)
        local oldlig=$(sed -n "s/^ligname //p" $inpdb)
        local oldcsts=$(awk '/>>> constraints >>>/ {f=1;next} /<<< constraints <<</ {f=0} f' $inpdb)
        [[ $ss_to_contigs == 1 && -z $oldcontigs ]] && {
            local oldcontigs=$(gen_contigs $fixbbres --inpdb $inpdb --ss_to_contigs --ss_trim 100)
        }
        [[ -z $fixedres && -n $oldfixed ]] && {
            echo "Using old fixedres residues found in input PDB: $oldfixed"
            local fixedres="$oldfixed"
            all_opts+=("--fixedres $fixedres")
        }
        [[ -z $ligname && -n $oldlig ]] && {
            echo "Using old ligands found in input PDB: $oldlig"
            local ligname=($oldlig)
            all_opts+=("--ligname ${ligname[*]}")
        }
    }
    
    # Generate contigs if one isn't provided
    [[ -z $contigs ]] && {
        echo "No contigs were provided. Generating contigs..."
        local contigs_opts=()
        [[ -n $oldcontigs && $oldcontigs != NONE ]] && {
            echo "Using old contigs found in input PDB for fixing backbone: $oldcontigs"
            local oldfixbb=$(echo $oldcontigs | sed 's|/|\n|g' | sed -n '/^[A-Z]/p' | paste -sd ' ')
            local fixbbres=$(get_newres "$oldcontigs" $oldfixbb | cut -d: -f2 | paste -sd ' ')
            local contigs_opts+=("--ss_to_contigs")
        }
        [[ -n $inpdb && $partial == 1 ]] && {
            local contigs_opts+=("--ss_to_contigs")
            [[ -z $ss_trim ]] && local contigs_opts+=("--ss_trim 100")
        }
        local contigs=$(gen_contigs $fixbbres ${contigs_opts[@]} ${all_opts[@]})
        [[ $contigs == FAILED ]] && {
            sleep 1 ; echo "Failed to find contigs" ; return
        }
    }
    echo "Using the following contigs for RFdiffusion: $contigs"
    
    # Generate fused ligands input PDB if ligands are provided
    [[ -n $inpdb ]] && {
        local inpdb_tmp="$tmpdir/$(basename $inpdb)"
        local ligstr=${ligname[*]}
        awk '
            $NF=="H" { next }
            $1=="ATOM" { print $0 ; next }
            $1=="HETATM" && $4~/^'"${ligstr// /|}"'$/ {
                print substr($0, 1, 17) "LIG X   1" substr($0, 27)
                next
            }
        ' $inpdb > $inpdb_tmp
    }
    
    # Prepare options for RFdiffusion
    local rfd_opts=()
    [[ -n $inpdb ]] && local rfd_opts+=("inference.input_pdb=$inpdb_tmp")
    
    if [[ $partial == 1 ]] ; then
        local rfd_opts+=("diffuser.partial_T=$timesteps")
    else
        local rfd_opts+=("diffuser.T=$timesteps")
    fi
    
    local ligname=(${ligname[@]})
    local contigs_fixbb=$(sed 's|/|\n|g' <<< "$contigs" | grep '^[A-Z]' | paste -sd ' ')

    [[ -n $ligname && $partial == 0 ]] && {
        local rfd_opts+=("potentials.guide_scale=1 potentials.substrate=LIG")
        local rfd_opts+=("'potentials.guiding_potentials=[\"type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1\"]'")
    }
    
    [[ $monomer_ROG == 1 && $partial == 0 ]] && {
        local rfd_opts+=("'potentials.guiding_potentials=[\"type:monomer_ROG,weight:1,min_dist:10\"]'")
        local rfd_opts+=("potentials.guide_scale=2 potentials.guide_decay=\"quadratic\"")
    }
    
    [[ -n $ppi_hotspots ]] && {
        local ppi_hotspots=$(expand_res_range $ppi_hotspots)
        local rfd_opts+=("'ppi.hotspot_res=[${ppi_hotspots// /,}]'")
    }
    
    [[ $ppi_mode == 1 ]] && {
        local rfd_opts+=("denoiser.noise_scale_ca=0.5" "denoiser.noise_scale_frame=0.5")
    }
    
    [[ -n $contigs_fixbb ]] && {
        local rfd_opts+=("inference.ckpt_override_path=$pixiroot/models/ActiveSite_ckpt.pt")
    }
    
    [[ $inpaint_seq == 1 && -n $contigs_fixbb ]] && {
        echo "The flag --inpaint_seq was invoked."
        local exp_fixed=$(expand_res_range $fixedres)
        local exp_fixbb=$(expand_res_range $contigs_fixbb)
        local inpainted=$(sed 's| |\n|g' <<< "$exp_fixbb" | awk -v sc_str="$exp_fixed" '
            BEGIN {
                n = split(sc_str, a, " ")
                for (i=1; i<=n; i++) sc[a[i]] = 1
            }
            $1 in sc { next } 1
        ' | paste -sd '/')
        [[ -n $inpainted ]] && {
            echo -n "The following residues will have their sequence hidden during diffusion: "
            sed 's|/| |g' <<< "$inpainted"
            local rfd_opts+=("'contigmap.inpaint_seq=[$inpainted]'")
        }
    }

    [[ -n $inpdb && ! -e $inpdb ]] && {
        echo "Input PDB does not exist..."
        return
    }

    # Run RFdiffusion
    local step1="$tmpdir/step1" ; mkdir -p $step1
    local rfd_suppress_logs=(
        "hydra.output_subdir=null"
        "'hydra.job_logging.root.handlers=[console]'"
        "~hydra.job_logging.handlers.file"
        "hydra/hydra_logging=disabled"
        "hydra.run.dir=."
    )
    if [[ -n $inpdb && $contigs == NONE ]] ; then
        echo "Input pdb detected and contigs set to NONE. Skipping diffusion."
        grep "^ATOM" $inpdb > $step1/Diffused_0.pdb
    else
        eval "
            $pixirun rfdiffusion rfdiffusion \
                'contigmap.contigs=[$contigs]' \
                inference.output_prefix=$step1/Diffused \
                inference.num_designs=1 \
                ${rfd_opts[*]} ${rfd_suppress_logs[*]}
        "
    fi
    
    # Rerun RFdiffusion if --backrub was enabled
    [[ -n $inpdb && -n $backrub ]] && {
        echo -n "The option --backrub was enabled. "
        echo "Running partial diffusion with $backrub timesteps..."
        if [[ -n $fixedres ]] ; then
            local newfixbb=$(get_newres "$contigs" $fixedres | cut -d: -f2 | paste -sd ' ')
            local re_contigs=$( 
                gen_contigs $newfixbb --ss_to_contigs --ss_trim 100 --inpdb $step1/Diffused_0.pdb
            )
        else
            # Without fixed residues, get contigs string that matches the exact length of the first chain
            re_contigs=$(
                awk '$1=="ATOM" && $3=="CA"' $step1/Diffused_0.pdb |
                awk '
                    ch == "" { ch = $5 }
                    ch { ri[NR] = $6 }
                    END { len = ri[NR] - ri[1] + 1 ; print len "-" len }
                '
            )
        fi
        eval "
            $pixirun rfdiffusion rfdiffusion \
                'contigmap.contigs=[$re_contigs]' \
                inference.input_pdb=$step1/Diffused_0.pdb \
                inference.output_prefix=$step1/ReDiffused \
                inference.num_designs=1 \
                diffuser.partial_T=${backrub} \
                inference.ckpt_override_path=$pixiroot/models/ActiveSite_ckpt.pt \
                ${rfd_suppress_logs[*]}
        "
        mv $step1/ReDiffused_0.pdb $step1/Diffused_0.pdb
    }
    
    local diffused=($step1/Diffused_*.pdb)
    [[ -s $diffused ]] || return

    # Prepare MPNN options
    local mpnn_opts=(${all_opts[*]} "--repack")
    
    # If --idealize is used, then idealize diffused backbone, then take out idealize flag
    [[ $idealize == 1 ]] && {
        mpros $diffused --idealize --inplace --skip_mpnn --skip_mpnn_score
        mpnn_opts=($(
            sed 's/--/\n--/g' <<< "${mpnn_opts[@]}" | awk 'NF && $1!="--idealize"'
        ))
    }

    # Align diffused output to original input PDB
    [[ -n $inpdb ]] && {
        local restab=$(get_newres "$contigs" $fixedres)
        local orires_list=$(echo "$restab" | cut -d: -f1 | paste -sd ' ')
    
        # Write pair fit string
        local fit_str=$(echo "$restab" | awk -F ':' '
            {
                ch[NR][1] = ri[NR][1] = $1
                ch[NR][2] = ri[NR][2] = $2
                sub(/[^A-Z]+$/, "", ch[NR][1])
                sub(/[^A-Z]+$/, "", ch[NR][2])
                sub(/^[^0-9]+/, "", ri[NR][1])
                sub(/^[^0-9]+/, "", ri[NR][2])
            }
            END {
                for (i=1; i<=NR; i++) {
                    for (j=i; j<=NR; j++) {
                        if (ch[j+1][1] != ch[j][1] || ch[j+1][2] != ch[j][2]) break
                        if (ri[j+1][1]-ri[j][1] != 1 || ri[j+1][2]-ri[j][2] != 1) break
                    }
                    
                    printf "!/diffused//%s/%s-%s/N+CA+C!,", ch[i][2], ri[i][2], ri[j][2]
                    printf "!/inpdb//%s/%s-%s/N+CA+C!\n",  ch[i][1], ri[i][1], ri[j][1]
                    i=j
                }
            }
        ' | sed "s/!/'/g" | paste -sd ',')
        
        # Add ligand
        local ligstr=${ligname[*]}
        [[ -n $ligname ]] && local add_ligand="or (inpdb and resn ${ligstr// /+})"
    
        # Pair fit
        sed -E 's/^ {12}//' << EOF | python
            import glob
            from pathlib import Path
            from pymol import cmd
            diffused_list = glob.glob('$step1/Diffused*.pdb')
            
            cmd.load('$inpdb', 'inpdb')
            diffused = diffused_list[0]
            cmd.load(diffused, 'diffused')
            pf = cmd.pair_fit($fit_str)
            cmd.create('obj01', 'diffused $add_ligand')
            cmd.save(diffused, 'obj01')
            cmd.delete('diffused')
            cmd.delete('obj01')
            with open(diffused, 'a') as f:
                f.write(f'\npf_rmsd {pf:.4f}')
            
            bn = Path(diffused).stem
            print(f'{bn} pf_rmsd {pf:.4f}')
EOF
    }
    
    # Get fixed residues for scafscore
    local scafscore_opts=()
    if [[ -n $ligname ]] ; then
        local scafscore_opts+=("--resn ${ligname[@]}")
    elif [[ -n $ppi_hotspots ]] ; then
        local hspots=$(get_newres "$contigs" | cut -d: -f2 | paste -sd ' ')
        local scafscore_opts+=("--chri $hspots")
    elif [[ -n $fixedres ]] ; then
        local fixednew=$(get_newres "$contigs" $fixedres | cut -d: -f2 | paste -sd ' ')
        local scafscore_opts+=("--chri $fixednew")
    fi
    
    # Calculate scafscore
    local diffused=($step1/Diffused_*.pdb)
    calc_perc_loop $diffused
    calc_scafscore $diffused ${scafscore_opts[@]}

    # Automatically calculate clash_cut
    [[ -n $ligname ]] && {
        if [[ -z $clash_cut ]] ; then
            local clash_cut_true=$(awk '
                $1=="HETATM" && $NF!="H" { n++ } END { print int(n/2) }
            ' $diffused)
            echo "Using --clash_cut $clash_cut_true (half of all non-hydrogen ligand atoms)"
        else
            local clash_cut_true="$clash_cut"
            echo "Using --clash_cut $clash_cut_true"
        fi
    }

    # Filter out
    local perc_loop=$(awk '$1=="perc_loop"{print$2}' $diffused)
    local rog_ala=$(awk '$1=="rog_ala"{print$2}' $diffused)
    [[ -n $ligname ]] && local clash=$(awk '$1=="clash"{print$2}' $diffused)
    if [[ $(echo "$rog_ala > $rog_cut" | bc) == 1 ]] ; then
        echo "Failed rog_ala cutoff of $rog_cut"
        rm $diffused && return
    elif [[ $(echo "$perc_loop > $loop_cut" | bc) == 1 ]] ; then
        echo "Failed perc_loop cutoff of $loop_cut"
        rm $diffused && return
    elif [[ -n $ligname && $(echo "$clash > $clash_cut_true" | bc) == 1 ]] ; then
        echo "Failed clash cutoff of $clash_cut_true"
        rm $diffused && return
    fi
    
    local diffused=($step1/Diffused_*.pdb)
    [[ -s "$diffused" ]] || {
        echo "No scaffolds passed filters. Exiting." && return
    }
    
    # Generate cstfile
    local tmpcst="$tmpdir/constraints.cst"
    [[ -n $cstfile ]] && cp $cstfile $tmpcst
    [[ -z $cstfile && -n $oldcsts ]] && {
        echo "Using old constraints found in input PDB: $(wc -l <<< "$oldcsts") constraints"
        echo "$oldcsts" > $tmpcst
    }
    [[ -n $inpdb && -n $fixedres && -n $fix_stdev ]] && {
        gen_coord_csts $inpdb --chri $fixedres --crd_stdev $fix_stdev >> $tmpcst
        echo "Generated coordinate constraints for $fixedres"
    }
    [[ -n $inpdb && -n $ligname && -n $lig_stdev ]] && {
        gen_coord_csts $inpdb --resn ${ligname[*]} --crd_stdev $lig_stdev >> $tmpcst
        echo "Generated coordinate constraints for ${ligname[*]}"
    }
    [[ -n $inpdb && -n $fixedres && -n $ligname && -n $ap_stdev ]] && {
        local ligchri=$(
            awk -v ligname="${ligname[*]}" '
                BEGIN {
                    n = split(ligname, a, " ")
                    for (i=1; i<=n; i++) lig_key[a[i]] = 1
                }
                $4 in lig_key { print $5 $6 }
            ' $inpdb | sort -u | paste -sd ' '
        )
        gen_atompair_csts $inpdb \
            --anchors $ligchri   \
            --targets $fixedres  \
            --ap_stdev $ap_stdev >> $tmpcst
        echo "Generated atompair constraints between ${ligname[*]} and $fixedres"
    }
    
    # Update cstfile with new residue numbering
    [[ -s $tmpcst ]] && {
        local garr=$(get_newres "$contigs" $fixedres | awk -F ':' -v dq='"' '{
            ch1 = substr($1,1,1) ; ch2 = substr($2,1,1)
            sub(/^[A-Z]/,"",$1)  ; sub(/^[A-Z]/,"",$2)
            printf "a[%s]=%s\n", dq $1 ch1 dq, dq $2 ch2 dq
        }' | paste -sd ';')
    
        awk -i inplace '
            BEGIN{'"$garr"'}
            $1=="CoordinateConstraint" { $3 = a[$3]!="" ? a[$3] : $3 }
            $1=="AtomPair" { $3 = a[$3]!="" ? a[$3] : $3 ; $5 = a[$5]!="" ? a[$5] : $5 }
            { print }
        ' $tmpcst
    
        echo "Adding a total of $(wc -l < $tmpcst) constraints as a cstfile..."
        local mpnn_opts=("--cstfile $tmpcst" ${mpnn_opts[@]})
    }
    
    # Update fixedres
    [[ -n $fixedres ]] && {
        local fixednew=$(get_newres "$contigs" $fixedres | cut -d: -f2 | paste -sd ' ')
        local mpnn_opts=("--fixedres $fixednew" ${mpnn_opts[@]})
    }
    
    # Add native biasjson
    [[ -n $natbias ]] && {
        local biasjson="$tmpdir/natbias.json"
        local newres=$(get_newres "$contigs" | cut -d: -f2 | paste -sd ' ')
        gen_natbiasjson ${diffused[0]} \
            --focus_set $newres \
            --global_weight $natbias > $biasjson
        local mpnn_opts=("--biasjson $biasjson" ${mpnn_opts[@]})
    }
    
    # Run MPNN-FastRelax
    local step2="$tmpdir/step2" ; mkdir -p $step2
    local topreu="1000"
    local topconf="0"
    local i
    for i in $(seq 1 $design_cycles) ; do

        # Run MPNN-Rosetta
        mpros $diffused --batch_size 1 --outdir $step2 ${mpnn_opts[@]}
        local design="$step2/$(basename ${diffused%.*})_0001.pdb"
        [[ -s $design ]] || continue
    
        # Check for selection metric
        local metpass="0"
        unset metarr
        [[ -n $select_met ]] && {
            echo "Processing the option: --select_met $select_met"
            local metarr=($(sed 's/[:|=]/ /g' <<< $select_met))  # (min/max/val metname ideal_val)
            if [[ ! ${metarr[0]} =~ ^(min|max|val)$ ]] ; then
                echo "WARNING: Could not parse selection logic! Refer to format in --help."
            elif [[ ${metarr[0]} == val && -z ${metarr[2]} ]] ; then
                echo "WARNING: Selection logic 'val' was invoked without an ideal value. Refer to format in --help."
            elif [[ -z ${metarr[1]} ]] ; then
                echo "WARNING: No metric detected! Refer to format in --help."
            elif ! grep "^${metarr[1]} " $design &>/dev/null ; then
                echo "WARNING: Could not find the metric ${metarr[1]} in $design"
                echo "List of available metrics:"
                awk '$1~/^[a-z]/ && $2==$2+0 {print$1}' $design
            else
                echo "Found the metric ${metarr[1]}" in $design.
                local metpass="1"
            fi
        }

        # Store best design
        local reu=$(awk '$1 == "fin_reu" {print$2}' $design)
        local conf=$(awk -v model="$model_type" '$1 == model "_score" {print$2}' $design)
        echo "fin_reu $topreu => $reu"
        echo "${model_type}_score $topconf => $conf"
        if [[ $i == 1 ]] ; then
            local topdes=$(cat $design) topreu="$reu" topconf="$conf"
            mv $design $diffused
        elif [[ $metpass == 1 ]] ; then
            echo "Accepting/rejecting designs using --select_met ${select_met}"
            local orimet=$(echo "$topdes" | awk -v met=${metarr[1]} '$1==met{print$2}')
            local finmet=$(awk -v met=${metarr[1]} '$1==met{print$2}' $design)
            echo "${metarr[1]} $orimet => $finmet"
            if [[ ${metarr[0]} == "min" ]] ; then
                echo "Selecting for lower values of ${metarr[1]}"
                if [[ $(echo "$finmet < $orimet" | bc) == 1 ]] ; then
                    echo "Score improved. Accepting new design."
                    local topdes=$(cat $design) topreu="$reu" topconf="$conf"
                    mv $design $diffused
                else
                    echo "Score worsened. Rejecting new design."
                fi
            elif [[ ${metarr[0]} == "max" ]] ; then
                echo "Selecting for higher values of ${metarr[1]}"
                if [[ $(echo "$finmet > $orimet" | bc) == 1 ]] ; then
                    echo "Score improved. Accepting new design."
                    local topdes=$(cat $design) topreu="$reu" topconf="$conf"
                    mv $design $diffused
                else
                    echo "Score worsened. Rejecting new design."
                fi
            elif [[ ${metarr[0]} == "val" ]] ; then
                echo "Selecting for values of ${metarr[1]} closest to ${metarr[2]}"
                local oridiff=$(awk "BEGIN{print ($orimet) - (${metarr[2]}) }" | sed 's/^-//')
                local findiff=$(awk "BEGIN{print ($finmet) - (${metarr[2]}) }" | sed 's/^-//')
                echo "diff ${metarr[1]} $oridiff => $findiff"
                if [[ $(echo "$findiff < $oridiff" | bc) == 1 ]] ; then
                    echo "Score improved. Accepting new design."
                    local topdes=$(cat $design) topreu="$reu" topconf="$conf"
                    mv $design $diffused
                else
                    echo "Score worsened. Rejecting new design."
                fi
            else
                echo "Something went wrong... Because John sucks at coding. Rejecting new design."
            fi
        elif [[ $(echo "$reu < $topreu" | bc) == 1 && $(echo "$conf > $topconf" | bc) == 1 ]] ; then
            echo "Both scores improved. Accepting new design."
            local topdes=$(cat $design) topreu="$reu" topconf="$conf"
            mv $design $diffused
        elif [[ $(echo "$reu < $topreu" | bc) == 1 || $(echo "$conf > $topconf" | bc) == 1 ]] ; then
            echo "Only one of the scores improved. Flipping coin to accept/reject."
            if [[ $(shuf -e -n1 1 2) == 1 ]] ; then
                echo "Accepting new design."
                local topdes=$(cat $design) topreu="$reu" topconf="$conf"
                mv $design $diffused
            else
                echo "Rejecting new design."
                rm $design
            fi
        else
            echo "Both scores worsened. Rejecting new design."
            rm $design
        fi
    done
    echo "$topdes" > $design
    
    local design=($step2/Diffused*.pdb)
    [[ -s $design ]] || {
        echo "Something went wrong with MPNN-FastRelax. Exiting." && return
    }
    
    # Record constraints
    [[ -s $tmpcst ]] && {
        echo "### >>> constraints >>> ###"
        cat $tmpcst
        echo "### <<< constraints <<< ###"
    } >> $design
    [[ -s $tmpcst ]] && calc_cst_rmsd $design --append
   
    # Output
    [[ -n $oldcontigs ]] && local contigs=$(update_contigs $oldcontigs $contigs)
    local bn=$(basename ${design%.*})
    calc_perc_loop $design $fixednew
    calc_scafscore $design ${scafscore_opts[@]}
    echo "contigs $contigs" | tee -a $design | sed "s|^|$bn |"
    [[ -n $inpdb ]] && {
        echo "inpdb $(basename ${inpdb%.*})" | tee -a $design | sed "s|^|$bn |"
    }
    [[ -n $fixednew ]] && {
        echo "fixedres $fixednew" | tee -a $design | sed "s|^|$bn |"
    }
    [[ -n $ligname ]] && {
        echo "ligname ${ligname[*]}" | tee -a $design | sed "s|^|$bn |"
    }

    grep "^ATOM" $design &>/dev/null && cp $design ${outprefix}.pdb
}
# <<< rfd_chain () <<< }}}

# Generate ids
if [[ $randsuffix == 1 ]] ; then
    i=$(awk "BEGIN{print(($numdes + 1) * 6)}")
    id=($(tr -cd 'A-Z' </dev/urandom | head -c ${i} | fold -w 6))
else
    id=($(seq -w 0 9999))
fi

# Loop through main function
i=1
while [[ $i -le $numdes ]] ; do

    # If output file already exists, skip to next round
    [[ -s ${outprefix}_${id[i]}.pdb ]] && {
        echo "${outprefix}_${id[i]}.pdb already exists, skipping..."
        ((i++)) ; continue
    }

    # Check for input directory
    [[ -n $indir ]] && {
        echo "Input directory (--indir) was specified."
        inpdb=$(shuf -e -n1 ${indir%/}/*.pdb)
        if [[ -z $inpdb ]] ; then
            echo "No PDBs detected in input directory." && exit
        else
            echo "Random pull as input PDB: $inpdb"
        fi
    }

    # Run rfd_chain function
    rfd_chain "$contigs" "${outprefix}_${id[i]}" ; rm -r $tmpdir/*

    # Check outputs
    if [[ -s ${outprefix}_${id[i]}.pdb ]] ; then
        echo "${outprefix}_${id[i]}.pdb was generated successfully."
        ((i++))
    else
        echo "${outprefix}_${id[i]}.pdb failed to generate."
    fi
done

runtime=($(awk "BEGIN{print($SECONDS/60,$SECONDS/60/$numdes)}"))
echo "Total runtime for $numdes designs: $runtime minutes"
echo "Total runtime per design: ${runtime[1]} minutes"
