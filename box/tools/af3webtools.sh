#!/usr/bin/env bash
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

usage () { cat << EOF
Usage: af3webtools SUBCMD INPUTS...
Command-line toolkit for preparing inputs and processing outputs for
alphafoldserver.com

Parameters:
    SUBCMD                  Subcommand for af3toolkit
                            Accepts 'prep' or 'unzip'
    INPUTS                  If SUBCMD is 'prep', then inputs can be .pdb and/or
                            .fasta files.
                            If SUBCMD is 'unzip', then inputs are the .zip files
                            from alphafoldserver.com

Options:
  --outjson [str]           Path to output json, which can be uploaded to the
                            AF3 server for job submission.
                            Only relevant for the 'prep' subcommand
  --refdir [str]            Path to reference directory, which should contain
                            the .pdb files for aligning AF3 outputs.
                            Only relevant for the 'unzip' subcommand
  --outdir [str]            Path to output directory, which will contain the
                            original .cif files from AF3, aligned .pdb and .pse
                            files.
                            Only relevant for the 'unzip' subcommand
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

val_opts=(outjson refdir outdir)
bool_opts=(help)

outjson="af3job_$(date +%y%m%d-%H%M%S).json"
refdir="."
outdir="."

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

# Return command line input
echo "Executing the follwing command:" "$0" "$@"

# Assign variables
subcmd=${args[1]}
infiles=${args[@]:2}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX")
trap 'rm -r $tmpdir' EXIT

# Check input command
[[ $subcmd =~ ^prep|unzip$ ]] || {
    echo "Input subcommand must be 'prep' or 'unzip'"
    exit
}

# --- pdb_to_csv () --- # {{{
pdb_to_csv () {
    local inpdb="$1"
    awk -v OFS=',' '
    BEGIN {
        prot["ALA"]="A" ; prot["ARG"]="R" ; prot["ASN"]="N" ; prot["ASP"]="D"
        prot["CYS"]="C" ; prot["GLU"]="E" ; prot["GLN"]="Q" ; prot["GLY"]="G"
        prot["HIS"]="H" ; prot["ILE"]="I" ; prot["LEU"]="L" ; prot["LYS"]="K"
        prot["MET"]="M" ; prot["PHE"]="F" ; prot["PRO"]="P" ; prot["SER"]="S"
        prot["THR"]="T" ; prot["TRP"]="W" ; prot["TYR"]="Y" ; prot["VAL"]="V"

        dna["DA"]="A" ; dna["DT"]="T" ; dna["DC"]="C" ; dna["DG"]="G"
        rna["A"]="A"  ; rna["U"]="U"  ; rna["C"]="C"  ; rna["G"]="G"

        lig["ADP"]="CCD_ADP" ; lig["ATP"]="CCD_ATP" ; lig["AMP"]="CCD_AMP"
        lig["GTP"]="CCD_GTP" ; lig["GDP"]="CCD_GDP" ; lig["FAD"]="CCD_FAD"
        lig["NAD"]="CCD_NAD" ; lig["NAP"]="CCD_NAP" ; lig["NDP"]="CCD_NDP"
        lig["HEM"]="CCD_HEM" ; lig["HEC"]="CCD_HEC" ; lig["PLM"]="CCD_PLM"
        lig["OLA"]="CCD_OLA" ; lig["MYR"]="CCD_MYR" ; lig["CIT"]="CCD_CIT"
        lig["CLA"]="CCD_CLA" ; lig["CHL"]="CCD_CHL" ; lig["BCL"]="CCD_BCL"
        lig["BCB"]="CCD_BCB"

        ion["MG"]="MG" ; ion["ZN"]="ZN" ; ion["CL"]="CL" ; ion["CA"]="CA"
        ion["NA"]="NA" ; ion["MN"]="MN" ; ion["K"]="K"   ; ion["FE"]="FE"
        ion["CU"]="CU" ; ion["CO"]="CO"
    }
    
    $1!~/^(ATOM|HETATM)$/ { next }  # Skip non ATOM or HETATM lines
    chri == $5$6 { next } { chri = $5$6 }  # Skip repeated chain+resi ID

    # Track chain
    ! ($5 in ch_key) { ch_key[$5] = 1 ; ch[++ch_i] = $5 }

    # Store sequence
    $4 in prot {
        type[ch_i] = "proteinChain"
        seq[ch_i] = seq[ch_i] prot[$4]
        next
    } $4 in dna {
        type[ch_i] = "dnaSequence"
        seq[ch_i] = seq[ch_i] dna[$4]
        next
    } $4 in rna {
        type[ch_i] = "rnaSequence"
        seq[ch_i] = seq[ch_i] rna[$4]
        next
    } $4 in lig {
        type[ch_i] = "ligand"
        seq[ch_i] = ( seq[ch_i] ? seq[ch_i] "," : "" ) lig[$4]
        ch_i++ ; next
    } $4 in ion {
        type[ch_i] = "ion"
        seq[ch_i] = ( seq[ch_i] ? seq[ch_i] "," : "" ) ion[$4]
        ch_i++ ; next
    }

    END {
        bn = FILENAME
        sub(/^.*[/]/, "", bn)
        sub(/[.]pdb$/, "", bn)
        print "name", bn
        for (i=1; i<=ch_i; i++) {
            if (seq[i]) { print type[i], seq[i] }
        }
    }
    ' $inpdb
}
# }}}

# --- run_prep () --- # {{{
run_prep () {
    local infiles="$1" outjson="$2"

    # Check inputs
    local infile
    for infile in $infiles ; do
        if [[ ${infile##*.} =~ ^fa|fasta$ ]] ; then
            local csvdat=$(awk -v OFS=',' '
                $0 !~ /^>/ { seq = seq $0 }
                END {
                    bn = FILENAME ; sub(/^.*[/]/, "", bn) ; sub(/[.](fa|fasta)$/, "", bn)
                    gsub(/ +/, "", seq)
                    print "name", bn
                    print "proteinChain", seq
                }
            ' $infile)
        elif [[ ${infile##*.} == pdb ]] ; then
            local csvdat=$(pdb_to_csv $infile)
        else
            continue
        fi
        echo "$csvdat" | jq -R -s '
            split("\n") | map(select(length>0) | split(",")) |
            ( map(select(.[0] == "name")) | .[0][1] ) as $jobName |
            map(select(.[0] != "name")) |
            map(
                if (.[0] | test("^proteinChain|dnaSequence|rnaSequence$")) then
                    { (.[0]): { "sequence": .[1], "count": 1 } }
                elif .[0] == "ligand" then
                    { (.[0]): { "ligand": .[1], "count": 1 } }
                elif .[0] == "ion" then
                    { (.[0]): { "ion": .[1], "count": 1 } }
                else
                    empty
                end
            ) as $sequences |
            {
                "name": $jobName,
                "modelSeeds": [1],
                "sequences": $sequences,
                "dialect": "alphafoldserver",
                "version": 1
            }
        '
    done | jq -s '.' > $outjson &&
    echo "Made $outjson"

    total_jobs=$(jq 'length' $outjson)
    batch_size=30
    for ((i=0; i<$total_jobs; i+=batch_size)); do
        batch_num=$(( i / batch_size + 1 ))
        jq ".[$i : $i+$batch_size]" ${outjson} > "${outjson%.*}_${batch_num}.json"
        echo "Created ${outjson%.*}_${batch_num}.json"
    done
    rm $outjson
}
# }}}

# --- run_unzip () --- # {{{
run_unzip () {
    local infiles="$1" refdir="$2" outdir="$3"
    local rawdir="$outdir/af3_outputs"
    local psedir="$outdir/aln_pses"
    local pdbdir="$outdir/aln_pdbs"
    mkdir -p $rawdir $psedir $pdbdir

    echo "Refdir set to $refdir. If PDBs were used as original input, make sure you set --refdir to the directory with PDBs."

    # Extract files
    local infile
    for infile in $infiles ; do
        [[ ${infile##*.} == zip ]] || continue
        unzip -o $infile -d $rawdir
    done

    # Initialize python script and log file
    local pyscript="$tmpdir/align.py"
    local logfile="$psedir/log.txt"
    echo -n > $logfile
    echo "from pymol import cmd" > $pyscript

    # Get prefixes, then align if they exist in refdir
    local prefixes=$(find $rawdir -maxdepth 2 -name "*_model_[0-4].cif" | sed 's/_model_[0-4].cif$//' | sort -u)
    local prefix
    for prefix in $prefixes ; do
        local bn=$(basename $prefix | sed 's/^fold_//')
        local oripdb=$(find $refdir -maxdepth 1 -iname "${bn}.pdb")
        ls "$oripdb" &>/dev/null || continue
        sed -E 's/^ {8}//' << EOF >> $pyscript
        cmd.load('$oripdb')
        for i in range(5):
            cmd.load(f'${prefix}_model_{i}.cif', f'af3model_{i}')
            cmd.alter(f'af3model_{i}', 'segi=""')
            aln = cmd.align(f'af3model_{i} and name CA', '$bn and name CA', cycles=0)
            cmd.align(f'af3model_{i} and name CA', '$bn and name CA', cycles=5)
            cmd.save(f'$pdbdir/${bn}_af3model_{i}.pdb', f'af3model_{i}')
            with open(f'$pdbdir/${bn}_af3model_{i}.pdb', 'a') as pdb:
                pdb.write(f'\nrmsd_ca {aln[0]}\n')

            print(f'${bn}_af3model_{i}.pdb rmsd_ca: {aln[0]}')
            with open(f'$logfile', 'a') as log:
                log.write(f'${bn}_af3model_{i}.pdb rmsd_ca: {aln[0]}\n')
        cmd.center()
        cmd.show('lines', 'not name N+C+O')
        cmd.set('sphere_scale', 0.5)
        cmd.hide('ev', 'elem H')
        cmd.save('$psedir/${bn}_af3.pse')
        cmd.delete('all')
EOF
    done

    # Run python script
    python $pyscript

    # Make ranked list
    sed 's/_af3model_[0-4].pdb rmsd_ca: / /' $logfile | awk '
    { sum[$1]+=$2 ; n[$1]++ } END {
        for (bn in sum) {
            print bn "_af3.pse", "rmsd_ca_avg:" , sum[bn]/n[bn]
        }
    }
    ' | sort -k3n | tee $psedir/ranked.txt
}
# }}}

# Main
if [[ $subcmd == prep ]] ; then
    run_prep "$infiles" "$outjson"
elif [[ $subcmd == unzip ]] ; then
    run_unzip "$infiles" "$refdir" "$outdir"
fi
