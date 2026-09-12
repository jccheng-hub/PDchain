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
Usage: mpros PDB...
Run LigandMPNN, then run Rosetta

Parameters:
    PDB                     Input PDB

General options:
  --fixedres [str]          List of fixed residues in chain+resi format
                            E.g. A7-10 A12 A19
  --redesres [str]          List of residues to design in chain+resi format
                            Will fix all other residues
                            E.g. A11 A13-18 A27
  --outdir [str]            Path to output directory
  --inplace                 Have output PDB replace input PDB
                            Forces --batch_size 1
  --help                    Display this help and exit

LigandMPNN options:
  --model_type [str]        Model to use. Choose from the following:
                            ligand_mpnn
                            protein_mpnn
                            soluble_mpnn
  --batch_size [int]        Number of designs to generate
  --temperature [float]     Temperature for design
  --biasjson [str]          Json file for biasing design
  --nbr_dist [float]        Distance in angstroms for neighbor detection.
                            When specifying a ligand with --ligname but not
                            using --model_type ligand_mpnn, residues within this
                            distance of the ligand will be redesigned with
                            --model_type ligand_mpnn.
  --repack                  Have MPNN repack every rotamer, including fixed ones
  --disallow_cys            Disallow cysteines during design
  --sc_context              Use side chain atoms as context
  --skip_mpnn               Skip MPNN design
  --skip_mpnn_score         Skip MPNN scoring step

PyRosetta options:
  --paramsdir [str]         Path to extra params directory
  --resfile [str]           Path to resfile
  --cstfile [str]           Path to cstfile
  --ca_stdev [float]        Strength/deviation of CA constraints
  --fix_stdev [float]       Strength/deviation of constraints on fixed residues
  --relax_repeats [int]     Number of relax repeats
  --ligname [str]           Three-letter ligand codes
                            Specify to add constraints during refinement
  --threads [int]           Threads to use during PyRosetta run.
  --idealize                Idealize bond lengths and bond angles
  --relax                   Refine with relax
  --design                  Refine with relax while allowing design
  --minimize                Minimize after refinement
  --ncaa_pal                Use ncaa packer palette
  --fix_bb                  Fix backbone during refinement
  --fix_chi                 Fix chi angles during refinement
  --ignore_metals           Do not set up metal-binding constraints
  --renumber_chain          Renumber each chain for the output pdb
  --keep_script             Keep the PyRosetta script used
  --skip_refine             Skip Rosetta refinement
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
    batch_size          biasjson            ca_stdev            cstfile             
    fix_stdev           fixedres            lig_stdev           ligname 
    model_type          outdir              redesres            nbr_dist
    relax_repeats       resfile             suffix              temperature         
    threads
)

bool_opts=(
    disallow_cys        design              fix_bb              fix_chi             
    idealize            ignore_metals       keep_script         ncaa_pal
    relax               renumber_chain      repack              sc_context
    skip_mpnn           skip_refine         minimize            inplace
    skip_mpnn_score     help
)

batch_size="1"
biasjson=""
ca_stdev=""
cstfile=""
fix_stdev=""
fixedres=""
lig_stdev=""
ligname=""
nbr_dist="8"
model_type="protein_mpnn"
outdir="."
redesres=""
relax_repeats="1"
resfile=""
suffix=""
temperature="0.1"
threads="1"

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

# Shell settings
shopt -s nullglob

# Print input command
echo "Executing command:" "$0" "$@"

# Input arguments
inpdb=(${args[@]:1})
outdir=$(realpath $outdir)

mkdir -p $outdir
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX) ; trap 'rm -r $tmpdir' EXIT

# --- Useful functions --- # {{{
expand_res_range () {
    echo $@ | sed 's/ /\n/g' | gawk -F '-' '{
        # Expected input per line: [chain][i]-[j] or [chain][i]
        ch = substr($0,1,1)   # Get chain 
        sub(/^[A-Z]/,"",$0)   # Get rid of chain letter
        $2 = $2 ? $2 : $1
        for (i=$1; i<=$2; i++) print ch i
    }' | paste -sd ' '
}

find_neighbors () {
    local inpdb="$1" fixedres="$2" lignames="$3"
    [[ -n $fixedres ]] && fixedres=$(expand_res_range $fixedres)
    awk -v fixedres="$fixedres" -v lignames="$lignames" -v dist="$nbr_dist" '
    BEGIN {
        if (dist == "") dist = 4
        if (fixedres) { split(fixedres, a, " ") ; for (i in a) fix[a[i]] }
        if (lignames) { split(lignames, a, " ") ; for (i in a) lig[a[i]] }
    }
    
    $1~/^(ATOM|HETATM)$/ && $NF!="H" {
        if ( (fixedres && $5$6 in fix) || (lignames && $4 in lig) ) {
            foc[$5][$6][$3][1] = $7
            foc[$5][$6][$3][2] = $8
            foc[$5][$6][$3][3] = $9
        } else {
            nfc[$5][$6][$3][1] = $7
            nfc[$5][$6][$3][2] = $8
            nfc[$5][$6][$3][3] = $9
        }
    }
    
    ENDFILE {
        find_nbr(foc, nfc, dist, nbr)
        for (ch in nbr) for (ri in nbr[ch]) print ch ri
        delete foc ; delete nfc ; delete nbr
    }
    
    function find_nbr(obj1, obj2, dcut, nbr,
                      dcut_sq, ch1, ch2, ri1, ri2, at1, at2, dx, dy, dz, d_sq)
    {
        # Find all atoms in obj2 within dcut of obj1, then store in nbr
        # obj1: object 1 | obj[chain][resi][atom][1/2/3] = x/y/z
        # obj2: object 2 | obj[chain][resi][atom][1/2/3] = x/y/z
        # dcut: distance cutoff | float
        # nbr: nbr array | nbr[chain][resi][atom][1/2/3] = x/y/z
        delete nbr ; dcut_sq = dcut * dcut
        for (ch1 in obj1) for (ri1 in obj1[ch1]) for (at1 in obj1[ch1][ri1]) {
            for (ch2 in obj2) for (ri2 in obj2[ch2]) for (at2 in obj2[ch2][ri2]) {
                dx = obj1[ch1][ri1][at1][1] - obj2[ch2][ri2][at2][1]
                dy = obj1[ch1][ri1][at1][2] - obj2[ch2][ri2][at2][2]
                dz = obj1[ch1][ri1][at1][3] - obj2[ch2][ri2][at2][3]
                d_sq = dx*dx + dy*dy + dz*dz
                if (d_sq < dcut_sq) {
                    nbr[ch2][ri2][at2][1] = obj2[ch2][ri2][at2][1]
                    nbr[ch2][ri2][at2][2] = obj2[ch2][ri2][at2][2]
                    nbr[ch2][ri2][at2][3] = obj2[ch2][ri2][at2][3]
                }
            }
        }
    }
    ' $inpdb | sort -k1,1.1 -k1.2n | uniq | paste -sd ' '
}
# }}}

# LigandMPNN options 
mpnn_opts=(
    "--pack_side_chains 1"
    "--number_of_packs_per_design 1"
    "--pack_with_ligand_context 1"
    "--temperature $temperature"
    "--repack_everything $repack"
    "--ligand_mpnn_use_side_chain_context $sc_context"
)
[[ -n $redesres ]] && {
    redesres=$(expand_res_range $redesres)
    fixedres=$(gawk '
        $1=="ATOM" && $3=="CA" && $5$6!~/^('${redesres// /|}')$/ {print$5$6}
    ' $inpdb | paste -sd ' ')
}
[[ -n $fixedres ]] && {
    fixedres=$(expand_res_range $fixedres)
    mpnn_opts+=("--fixed_residues \"$fixedres\"")
}
[[ -n $biasjson ]] && {
    biasjson=$(realpath $biasjson)
    echo "Using json for biasing residues $biasjson"
    mpnn_opts+=("--bias_AA_per_residue \"$biasjson\"")
}
[[ $disallow_cys == 1 ]] && {
    echo "Excluding cysteines from sequence design"
    mpnn_opts+=("--omit_AA \"C\"")
}

# --- run_mpnn () --- # {{{
run_mpnn () {
    local inpdbs="$1" outdir=$2
    local injson="$tmpdir/inpdbs.json"
    sed 's/ /\n/g' <<< "$inpdbs" |
    jq -R -s 'split("\n") | map(select(length>0)) | map({ (.): "" }) | add' > $injson

    eval "
        $pixirun ligandmpnn ligandmpnn \
            --pdb_path_multi $injson \
            --out_folder $outdir \
            --model_type $model_type \
            --batch_size $batch_size \
            ${mpnn_opts[@]}
    "

    # If ligname is provided, then use LigandMPNN on residues within contact distance of ligand
    [[ -n $ligname && $model_type != "ligand_mpnn" ]] && {
        echo "Ligands are specified, but the --model_type is NOT ligand_mpnn."
        echo "Will now run LigandMPNN on residues within ${nbr_dist} angstroms of ligands."

        [[ $(echo $inpdbs | wc -w) -gt 1 ]] && {
            echo "WARNING: Multiple input PDBs were specified. Unless all input PDBs are homologous and have the exact numbering, this feature will not work as intended!"
        }

        # Move initial outputs into temporary dir
        local redesdir="$tmpdir/redesdir" ; mkdir -p $redesdir
        mv $outdir/packed $outdir/backbones $outdir/seqs $redesdir

        # Get additional fixedres
        local nbrres=$(find_neighbors "$(echo $redesdir/packed/*.pdb)" "" "$ligname")
        local finres=$(
            echo $nbrres        | sed 's/ /\n/g'                       |
            sort -k1,1.1 -k1.2n | awk '$1!~/^('"${fixedres// /|}"')$/' |
            paste -sd ' '
        )
        echo "The following residues were within ${nbr_dist} angstroms of $ligname: $nbrres"
        echo "List of fixed residues: ${fixedres}"
        echo "List of residues to redesign: ${finres}"

        # Update injson and fixedres
        mpnn_opts=($(echo ${mpnn_opts[@]} | sed 's/--/\n--/g' | grep -v '^--fixed_residues' | paste -sd ' '))
        mpnn_opts+=("--redesigned_residues \"$finres\"")
        echo $redesdir/packed/*.pdb | sed 's/ /\n/g' |
        jq -R -s 'split("\n") | map(select(length>0)) | map({ (.): "" }) | add' > $injson

        # Run LigandMPNN if neighboring residues exist. If not, copy back the original output.
        if [[ -n $finres ]] ; then
            eval "
                $pixirun ligandmpnn ligandmpnn \
                    --pdb_path_multi $injson \
                    --out_folder $outdir \
                    --model_type ligand_mpnn \
                    --batch_size 1 \
                    ${mpnn_opts[@]}
            "

            # Rename outputs
            local outpdb
            for outpdb in $outdir/packed/*.pdb ; do
                mv $outpdb ${outpdb%_packed_1_1.pdb}.pdb
            done
        else
            echo "No residues were eligible for LigandMPNN redesign. Retaining original designs."
            mv $redesdir/packed $redesdir/backbones $redesdir/seqs $outdir
        fi
    }
}
# }}}

# --- run_pyrosetta () --- # {{{
run_pyrosetta () {
    local indir=$1 outdir=$2
    
    mkdir -p $outdir
    local tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX)
    
    # Make initial template
    local pyrosetta_script="$tmpdir/run_pyrosetta.py"
    sed -E 's/^ {4}//' << EOF > $pyrosetta_script
    #!/usr/bin/env python
    import shutil, pyrosetta
    from pathlib import Path
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
    
    pyrosetta.init('ROSETTA_CMD_OPTIONS')
    inpdbs  = [INPDB_PYLIST]
    outpdbs = [OUTPDB_PYLIST]
    fixedres = 'COMMA_LIST_FIXEDRES'
    ligands = 'COMMA_LIST_LIGANDS'
    ncaas   = 'COMMA_LIST_NCAAS'
    resfcmd = 'RESFILE_CMD'
    resfile = 'PATH_TO_RESFILE'
    cstfile = 'PATH_TO_CSTFILE'
    
    xml_string = f'''
    <ROSETTASCRIPTS>
        <SCOREFXNS>
            <ScoreFunction name="r15" weights="ref2015.wts"/>
            <ScoreFunction name="r15_cst" weights="ref2015_cst.wts"/>
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
            <True name="FullPose"/>
             Index name="FixedRes" resnums='{fixedres}'/>
             ResidueName name="Ligs" residue_name3="{ligands}"/>
        </RESIDUE_SELECTORS>
        <PACKER_PALETTES>
            <DefaultPackerPalette name="pakpal"/>
             CustomBaseTypePackerPalette name="pakpal" additional_residue_types="{ncaas}"/>
        </PACKER_PALETTES>
        <TASKOPERATIONS>
            <IncludeCurrent name="icr"/>
             ResfileCommandOperation name="maintask" residue_selector="FullPose" command="{resfcmd}"/>
             ReadResfile name="maintask" filename="{resfile}"/>
        </TASKOPERATIONS>
        <MOVE_MAP_FACTORIES>
            <MoveMapFactory name="mmf" bb="$((1-fix_bb))" chi="$((1-fix_chi))"/>
        </MOVE_MAP_FACTORIES>
        <SIMPLE_METRICS>
            <TotalEnergyMetric name="total_reu" scoretype="total_score" residue_selector="FullPose" scorefxn="r15"/>
             TotalEnergyMetric name="fixed_reu" scoretype="total_score" residue_selector="FixedRes" scorefxn="r15"/>
             SasaMetric name="sasa_metric" residue_selector="FixedRes" sasa_metric_mode="all_sasa"/>
             HbondMetric name="hbond_metric" output_as_pdb_nums="true" residue_selector="Ligs"/>
        </SIMPLE_METRICS>
        <FILTERS>
             DSasa name="dsasa" lower_threshold="0" upper_threshold="9999"/>
             Ddg name="ddg" threshold="9999" repeats="5"/>
        </FILTERS>
        <MOVERS>
            <VirtualRoot name="add_vroot"/>
            <SetupMetalsMover name="setup_metals" metals_detection_LJ_multiplier="1" metals_distance_constraint_multiplier="$((1-ignore_metals))" metals_angle_constraint_multiplier="$((1-ignore_metals))"/>
             AddConstraints name="add_csts">
                 FileConstraintGenerator name="filecst" filename="{cstfile}"/>
             /AddConstraints>
             Idealize name="idealize"/>
            <MinMover name="minimize" movemap_factory="mmf" scorefxn="r15_cst"/>
             FastDesign name="fastdes" task_operations="icr,maintask" packer_palette="pakpal" movemap_factory="mmf" scorefxn="r15_cst"/>
        </MOVERS>
        <PROTOCOLS>
            <Add metrics="total_reu" labels="ori_reu"/>
            <Add mover_name="add_vroot"/>
            <Add mover_name="setup_metals"/>
             Add mover_name="add_csts"/>
             Add mover_name="idealize"/>
             Add mover_name="fastdes"/>
             Add mover_name="minimize"/>
            <Add metrics="total_reu" labels="fin_reu"/>
             Add metrics="fixed_reu" labels="fixedres_reu"/>
             Add metrics="sasa_metric" labels="fixedres_sasa"/>
             Add metrics="hbond_metric" labels="hbonds_to_focus"/>
             Add filter_name="dsasa"/>
             Add filter_name="ddg"/>
        </PROTOCOLS>
        <OUTPUT scorefxn="r15"/>
    </ROSETTASCRIPTS>
    '''
    
    # Run PyRosetta
    def pyrosetta_xml(inpdb, outpdb):
        pose = pyrosetta.pose_from_pdb(inpdb)
        xml = XmlObjects.create_from_string(xml_string).get_mover('ParsedProtocol')
        xml.apply(pose)
        pose.dump_pdb(outpdb)
    
    cstdir = Path(cstfile).parent
    for i in range(len(inpdbs)):
        inpdb = inpdbs[i]
        outpdb = outpdbs[i]
    
        # Update cst based on input PDB
        bn = inpdb.split('/')[-1].split('.')[0]
        shutil.move(f'{cstdir}/{bn}.cst', cstfile)
    
        pyrosetta_xml(inpdb, outpdb)
EOF

    local sed_script="$tmpdir/sub_pyrosetta.sed"
    
    # Prepare Rosetta command line options
    local rosetta_opts=(
        "-in:ignore_waters false"
        "-in:ignore_unrecognized_res false"
        "-run:preserve_header true"
        "-multithreading:total_threads $threads"
        "-relax:default_repeats $relax_repeats"
    )
    [[ $extra_chis == 1 ]] && rosetta_opts+=("-ex1" "-ex2")
    [[ $renumber_chain == 1 ]] && rosetta_opts+=("-renumber_pdb" "-per_chain_renumbering")
    
    # Add params from paramsdir
    local params=(${paramsdir%/}/*.params)
    [[ -n ${params[@]} ]] && {
        ncaas=$( 
            grep -H '^TYPE POLYMER' ${paramsdir%/}/*.params | cut -d: -f1 |
            xargs grep -h '^NAME' | gawk '{print$2}' | paste -sd ' '
        )
        [[ -n $ncaas ]] && {
            for ncaa in $ncaas ; do
                rotlib="$paramsdir/$ncaa.rotlib"
                [[ ! -f $rotlib ]] && echo "Could not find $rotlib" && exit
                rotlibpath=$(grep '^NCAA_ROTLIB_PATH ' $paramsdir/$ncaa.params | gawk '{print$2}')
                if [[ ! -f $rotlibpath ]] ; then
                    echo "$ncaa.rotlib does not exist at $rotlibpath"
                    echo "Editing params file so that it points to $rotlib"
                    contents=$(sed -E "s;^(NCAA_ROTLIB_PATH ).*$;\1$rotlib;" $paramsdir/$ncaa.params)
                    [[ -n $contents ]] && echo "$contents" > $paramsdir/$ncaa.params
                fi
            done
        }
        rosetta_opts+=("-extra_res_fa ${params[@]}")
    }
    echo "s|ROSETTA_CMD_OPTIONS|${rosetta_opts[@]}|" >> $sed_script
    
    # Format python lists
    local inpdb_list=(${indir%/}/*.pdb)
    local inpdb_pylist=$( 
        sed 's/ /\n/g' <<< ${inpdb_list[@]} |
        sed -E "s/^|$/'/g" | paste -sd ','
    )
    local outpdb_pylist=$( 
        basename -a ${inpdb_list[@]} | sed "s|^|${outdir%/}/|" |
        sed -E "s/^|$/'/g" | paste -sd ','
    )
    echo "s|INPDB_PYLIST|$inpdb_pylist|"   >> $sed_script
    echo "s|OUTPDB_PYLIST|$outpdb_pylist|" >> $sed_script
    
    # Packer palette
    [[ $ncaa_pal == 1 && -n $ncaas ]] && {
        echo "Including ncaa packer palette."
        echo 's|<(DefaultPackerPalette name="pakpal")| \1|'         >> $sed_script
        echo 's| (CustomBaseTypePackerPalette name="pakpal")|<\1|'  >> $sed_script
        echo "s|COMMA_LIST_NCAAS|${ncaas// /,}|"                    >> $sed_script
    }
    
    # Add constraints
    local tmpcst=$(mktemp $tmpdir/tmp_XXXXXX.cst)
    [[ -n $cstfile ]] && cp $cstfile $tmpcst
    echo "s|PATH_TO_CSTFILE|$tmpcst|" >> $sed_script
    
    local pdb
    for pdb in ${inpdb[@]} ; do
        for i in $(seq 1 $batch_size) ; do
            local pdbcst="$tmpdir/$(basename ${pdb%.*})_packed_${i}_1.cst"
            cp $tmpcst $pdbcst
            [[ -n $ca_stdev ]] && {
                gen_coord_csts $pdb --name CA --crd_stdev $ca_stdev >> $pdbcst
            }
            [[ -n $fix_stdev && -n $fixedres ]] && {
                gen_coord_csts $pdb --chri $fixedres --crd_stdev $fix_stdev >> $pdbcst
            }
            [[ -n $lig_stdev && -n $ligname ]] && {
                gen_coord_csts $pdb --resn $ligname --crd_stdev $lig_stdev >> $pdbcst
            }
        done
    done
    
    [[ -n $cstfile || -n $ca_stdev || ( -n $ligname && -n $lig_stdev ) || -n $fix_stdev ]] && {
        echo "Including constraints in pyrosetta script."
        echo 's| (AddConstraints name="add_csts")|<\1|'         >> $sed_script
        echo 's| (FileConstraintGenerator name="filecst")|<\1|' >> $sed_script
        echo 's| (/AddConstraints)|<\1|'                        >> $sed_script
        echo 's| (Add mover_name="add_csts")|<\1|'              >> $sed_script
    }

    # Additional metrics for ligands
    [[ -n $ligname ]] && {
        echo "Including additional metrics for ligands."
        echo 's| (ResidueName name="Ligs")|<\1|'            >> $sed_script
        echo 's| (HbondMetric name="hbond_metric")|<\1|'    >> $sed_script
        echo 's| (Add metrics="hbond_metric")|<\1|'         >> $sed_script
        echo 's| (DSasa name="dsasa")|<\1|'                 >> $sed_script
        echo 's| (Add filter_name="dsasa")|<\1|'            >> $sed_script
        echo "s|COMMA_LIST_LIGANDS|${ligname// /,}|"        >> $sed_script
    }

    # Addition metric for ppi
    [[ $(awk '$1~/^ATOM|HETATM$/ && !ch[$5] {ch[$5]++; n++} END {print n}' $inpdb) -gt 1 ]] && {
        echo "Including metrics for multiple chains."
        echo 's| (Ddg name="ddg")|<\1|'                     >> $sed_script
        echo 's| (Add filter_name="ddg")|<\1|'              >> $sed_script
    }    

    # Additional metrics for fixed residues
    [[ -n $fixedres ]] && {
        echo "Including metrics for fixed residues."
        local invfixed=$(
            echo $fixedres | sed 's/ /\n/g' | awk '{
                ch=$0 ; sub(/[^A-Z]+$/, "", ch) ; sub(/^[A-Z]+/, "", $0) ; print $0 ch
            }' | paste -sd ','
        )
        echo 's| (Index name="FixedRes")|<\1|'              >> $sed_script
        echo 's| (SasaMetric name="sasa_metric")|<\1|'      >> $sed_script
        echo 's| (TotalEnergyMetric name="fixed_reu")|<\1|' >> $sed_script
        echo 's| (Add metrics="sasa_metric")|<\1|'          >> $sed_script
        echo 's| (Add metrics="fixed_reu")|<\1|'            >> $sed_script
        echo "s|COMMA_LIST_FIXEDRES|$invfixed|"             >> $sed_script
    }

    # Idealize
    [[ $idealize == 1 ]] && {
        echo "Including idealize mover."
        echo 's| (Idealize name="idealize")|<\1|'  >> $sed_script
        echo 's| (Add mover_name="idealize")|<\1|' >> $sed_script
    }
    
    # Refinement
    if [[ -n $resfile ]] ; then 
        echo "Including design mover with resfile."
        echo 's| (ReadResfile name="maintask")|<\1|' >> $sed_script
        echo 's| (FastDesign name="fastdes")|<\1|'   >> $sed_script
        echo 's| (Add mover_name="fastdes")|<\1|'    >> $sed_script
        echo "s|PATH_TO_RESFILE|$resfile|"           >> $sed_script
    elif [[ $design == 1 ]] ; then
        echo "Including design mover."
        echo 's| (ResfileCommandOperation name="maintask")|<\1|' >> $sed_script
        echo 's| (FastDesign name="fastdes")|<\1|'               >> $sed_script
        echo 's| (Add mover_name="fastdes")|<\1|'                >> $sed_script
        echo "s|RESFILE_CMD|ALLAA|"                              >> $sed_script
    elif [[ $relax == 1 ]] ; then
        echo "Including relax mover."
        echo 's| (ResfileCommandOperation name="maintask")|<\1|' >> $sed_script
        echo 's| (FastDesign name="fastdes")|<\1|'               >> $sed_script
        echo 's| (Add mover_name="fastdes")|<\1|'                >> $sed_script
        echo "s|RESFILE_CMD|NATAA|"                              >> $sed_script
    fi
    
    # Minimize
    [[ $minimize == 1 ]] && {
        echo "Including minimize mover."
        echo 's| (Add mover_name="minimize")|<\1|' >> $sed_script
    }

    # Make substitutions, then run pyrosetta script
    sed -i -E -f $sed_script $pyrosetta_script
    $pixirun pyrosetta python $pyrosetta_script
    [[ $keep_script == 1 ]] && cp $pyrosetta_script $outdir
    
    # Calculate energy per res
    local outpdb
    for outpdb in $outdir/*.pdb ; do
        local nres=$(gawk '$1~/^(ATOM|HETATM)$/{print $5$6}' $outpdb | sort -u | wc -l)
        local met
        for met in ori_reu fin_reu ; do
            val=$(sed -n "s/^${met} //p" $outpdb)
            val=$(gawk "BEGIN{print($val/$nres)}")
            echo ${met}_per_res $val >> $outpdb
            echo $(basename $outpdb) ${met}_per_res $val
        done
    done
    
    rm -r $tmpdir
}
# }}}

# Check for --inplace flag
[[ $inplace == 1 ]] && batch_size=1

# Run LigandMPNN
step1="$tmpdir/step1" ; mkdir -p $step1
if [[ $skip_mpnn == 0 ]] ; then
    run_mpnn "${inpdb[*]}" "$step1"
else
    mkdir -p $step1/packed
    for pdb_i in ${inpdb[@]} ; do
        bnpdb=$(basename ${pdb_i%.*})
        for i in $(seq 1 $batch_size) ; do
            awk '$1~/^(ATOM|HETATM)$/' $pdb_i > $step1/packed/${bnpdb}_packed_${i}_1.pdb
        done
    done
fi

# Run PyRosetta for refinement
step2="$tmpdir/step2" ; mkdir -p $step2
if [[ $skip_refine == 0 ]] ; then
    run_pyrosetta $step1/packed $step2
else
    cp $step1/packed/*.pdb $step2
fi

# Score outputs with LigandMPNN
[[ $skip_mpnn_score == 0 ]] &&
mpnn_score ${step2}/*.pdb --model_type ${model_type}

# Output
for i_pdb in ${inpdb[@]} ; do
    i=1; j=$(printf "%04d" $i)
    bn=$(basename ${i_pdb%.*})
    f_pdb=($step2/${bn}*.pdb)
    for pdb in ${f_pdb[@]} ; do
        [[ $inplace == 1 ]] && { mv $pdb $i_pdb ; continue ; }
        mv $pdb $outdir/${bn}${suffix}_${j}.pdb
        i=$((i+1)) ; j=$(printf "%04d" $i)
    done
done
[[ $keep_script == 1 ]] && cp $step2/run_pyrosetta.py $outdir
