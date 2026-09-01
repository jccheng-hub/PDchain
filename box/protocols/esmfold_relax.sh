#!/usr/bin/env bash
usage () { cat << EOF
Usage: esmfold_relax PDB/FASTA...
Uses ESMFold webserver for structure prediction, then relax with PyRosetta.

Once a sequence is predicted, the raw prediction will be stored locally.
If the same sequence is requested again, then the local PDB will be retrieved
instead of requesting the ESMFold webserver.

Parameters:
    PDB                     Input pdb file
    FASTA                   Input fasta file

Options:
  --outdir [str]            Output directory
  --foldrepo [str]          Directory for storing newly predicted structures and
                            fetching previously predicted structures
  --ca_stdev [float]        Strength of CA constraints during relax
  --skip_refine             Skip PyRosetta refinement
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

val_opts=(ca_stdev foldrepo outdir)
bool_opts=(skip_refine help)

ca_stdev="1"
outdir="."
foldrepo=""

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

# Check Pixi
pixiroot="$PIXI_PROJECT_ROOT"
pixitoml="$pixiroot/pixi.toml"
pixirun="pixi run -m $pixitoml -e"
[[ -n $pixiroot ]] || {
    echo 'Need $PIXI_PROJECT_ROOT to be defined.'
    exit
}

# Functions
pdb_to_fasta () {
    local pdb=$1 ch=$2
    local bn="$(basename ${pdb%.*})_$ch"
    awk '
    BEGIN {
        olc["ALA"]="A" ; olc["ARG"]="R" ; olc["ASN"]="N" ; olc["ASP"]="D"
        olc["CYS"]="C" ; olc["GLU"]="E" ; olc["GLN"]="Q" ; olc["GLY"]="G"
        olc["HIS"]="H" ; olc["ILE"]="I" ; olc["LEU"]="L" ; olc["LYS"]="K"
        olc["MET"]="M" ; olc["PHE"]="F" ; olc["PRO"]="P" ; olc["SER"]="S"
        olc["THR"]="T" ; olc["TRP"]="W" ; olc["TYR"]="Y" ; olc["VAL"]="V"
        print ">'$bn'"
    }

    $1~/^(ATOM|HETATM)$/ && $3=="CA" && $5=="'$ch'" {
        o = (olc[$4]) ? olc[$4] : "X" ; printf o
    } END {print ""}' $pdb
}

# Initialize
infiles=${args[@]:1}
outdir=${outdir%/} ; mkdir -p $outdir
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX) ; trap 'rm -r $tmpdir' EXIT

# Keep only input files with recognized extension
infiles=($(sed 's/ /\n/g' <<< "$infiles" | grep -E '[.](fasta|fa|pdb)$'))

# Report inputs
if [[ ${#infiles[@]} -gt 0 ]] ; then
    echo "Folding the following with ESMFold webserver:"
    sed 's/ /\n/' <<< "${infiles[*]}"
else
    echo "No input .pdb, .fasta, or .fa files detected."
    exit
fi

# If foldrepo is not given, then make it in $pixiroot
[[ -z $foldrepo ]] && {
    foldrepo="$pixiroot/foldrepo/esmfold"
    echo "The option --foldrepo was not specified."
}
mkdir -p $foldrepo
echo "Setting foldrepo to $foldrepo."

# Run ESMFold for all inputs
esmdir="$tmpdir/esmfold" ; mkdir -p $esmdir
for infile in ${infiles[@]} ; do

    # Get sequence
    if [[ $infile =~ [.]pdb$ ]] ; then
        ch=$(awk '$1=="ATOM" {print$5; exit}' $infile)
        seq=$(pdb_to_fasta $infile $ch | tail -n1)
    elif [[ $infile =~ [.]fasta$ || $infile =~ [.]fa$ ]] ; then
        seq=$(awk '/^>/ {x=1-x; next} x' $infile | paste -sd '')
    fi
    seqlen=$(echo -n $seq | wc -c)
    echo "Submitting sequence ($seqlen residues) from $infile."
    echo "$seq"

    # Check foldrepo
    prevfold=""
    prevfold=($(grep -H "^seq $seq$" $foldrepo/*.pdb 2>/dev/null | cut -d: -f1))

    # Fold
    bn=$(basename ${infile%.*})
    esmpdb="$esmdir/${bn}_esmfold.pdb"
    if [[ -s $prevfold ]] ; then
        echo -e "Found $prevfold with identical sequence\nSkipping ESMFold webserver request"
        grep "^ATOM" $prevfold > $esmpdb
    else
        curl -k -X POST --data "$seq" https://api.esmatlas.com/foldSequence/v1/pdb/ > $esmpdb
    fi

    # Calculate pLDDT
    atomlin=$(grep "^ATOM" $esmpdb 2>/dev/null | wc -l)
    if [[ $atomlin -gt 0 ]] ; then
        awk -v esmpdb="$esmpdb" '
        $1=="ATOM" && $3=="CA" { s+=$(NF-1); n++ }
        END {
            plddt_ca = 100 * s/n
            printf "\nplddt_ca %s\n", plddt_ca >> esmpdb
        }
        ' $esmpdb
        echo "seq $seq" >> $esmpdb
    else
        echo "ESMFold failed for $infile..."
        rm $esmpdb
        continue
    fi

    # Store into foldrepo
    [[ $atomlin -gt 0 && -z $prevfold ]] && {
        cp $esmpdb $foldrepo && echo "Stored predicted structure at $foldrepo/$(basename $esmpdb)"
    }

    # Generate constraints for later Rosetta refinement
    awk '
    $1=="ATOM" && $3=="CA" {
        printf "CoordinateConstraint %s %s CA 1 ", $3, $6 $5
        printf "%s %s %s HARMONIC 0 '$ca_stdev'\n", $7, $8, $9
    }
    ' $esmpdb > ${esmpdb%.*}.cst
done

# >>> pyrosetta_relax () >>> {{{
pyrosetta_relax () {
    local inpdbs=($1) outdir=$2
    local pyscript="$tmpdir/relax.py"
    mkdir -p $outdir
    sed -E 's/^ {4}//' << EOF | $pixirun pyrosetta python
    import pyrosetta
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

    def get_xml_string(cstfile):
        return f'''
        <ROSETTASCRIPTS>
            <SCOREFXNS>
                <ScoreFunction name="r15" weights="ref2015.wts"/>
                <ScoreFunction name="r15_cst" weights="ref2015_cst.wts"/>
            </SCOREFXNS>
            <RESIDUE_SELECTORS>
                <True name="FullPose"/>
            </RESIDUE_SELECTORS>
            <PACKER_PALETTES>
            </PACKER_PALETTES>
            <TASKOPERATIONS>
                <IncludeCurrent name="icr"/>
                <ResfileCommandOperation name="maintask" residue_selector="FullPose" command="NATAA"/>
            </TASKOPERATIONS>
            <MOVE_MAP_FACTORIES>
                <MoveMapFactory name="mmf" bb="1" chi="1"/>
            </MOVE_MAP_FACTORIES>
            <SIMPLE_METRICS>
                <TotalEnergyMetric name="total_reu" scoretype="total_score" residue_selector="FullPose" scorefxn="r15"/>
                <SecondaryStructureMetric name="ss_metric" residue_selector="FullPose"/>
            </SIMPLE_METRICS>
            <FILTERS>
            </FILTERS>
            <MOVERS>
                <VirtualRoot name="add_vroot"/>
                <AddConstraints name="add_csts">
                    <FileConstraintGenerator name="filecst" filename="{cstfile}"/>
                </AddConstraints>
                <MinMover name="minimize" movemap_factory="mmf" scorefxn="r15_cst"/>
                <FastDesign name="fastdes" task_operations="icr,maintask" movemap_factory="mmf" scorefxn="r15_cst"/>
            </MOVERS>
            <PROTOCOLS>
                <Add metrics="total_reu" labels="ori_reu"/>
                <Add mover_name="add_vroot"/>
                <Add mover_name="add_csts"/>
                <Add mover_name="fastdes"/>
                <Add metrics="total_reu" labels="fin_reu"/>
            </PROTOCOLS>
            <OUTPUT scorefxn="r15"/>
        </ROSETTASCRIPTS>
        '''

    # Define input variables
    input_str = '${inpdbs[*]}'
    inpdbs = input_str.split()
    outdir = '$outdir'

    # Initialize PyRosetta
    pyrosetta.init('-run:preserve_header -relax:default_repeats 1')

    # Loop through all input PDBs
    for inpdb in inpdbs:

        # Update xml_string with cstfile based on input pdb
        cstfile = inpdb.rsplit('.', 1)[0] + '.cst'
        xml_string = get_xml_string(cstfile)

        # Run PyRosetta
        pose = pyrosetta.pose_from_pdb(inpdb)
        xml = XmlObjects.create_from_string(xml_string).get_mover('ParsedProtocol')
        xml.apply(pose)

        # Generate output pdb name, then save pose as output
        bn = inpdb.rsplit('/', 1)[-1].rsplit('.', 1)[0]
        outpdb = outdir + '/' + bn + '_rlx.pdb'
        pose.dump_pdb(outpdb)
EOF
}
# <<< pyrosetta_relax () <<< }}}

# Relax with PyRosetta
rlxdir="$tmpdir/relaxed" ; mkdir -p $rlxdir
esmpdbs=($esmdir/*.pdb)
if [[ $skip_refine == 0 ]] ; then
    pyrosetta_relax "${esmpdbs[*]}" $rlxdir
else
    echo "The option --skip_refine was invoked. Skipping Rosetta refinement."
    for esmpdb in ${esmpdbs[@]} ; do
        cp $esmpdb $rlxdir/$(basename ${esmpdb%.*})_rlx.pdb
    done
fi

# Align with PyMOL
alndir="$tmpdir/aligned" ; mkdir -p $alndir
inpdbs=($(sed 's/ /\n/g' <<< "${infiles[*]}" | grep '[.]pdb$' | paste -sd ' '))
python << EOF
from pathlib import Path
from pymol import cmd

inpdb_str = '${inpdbs[*]}'
mobdir = '$rlxdir'
suffix = '_esmfold_rlx'
alndir = '$alndir'

inpdbs = inpdb_str.split()
for target in inpdbs:
    bn = target.split('/')[-1].rsplit('.', -1)[0]
    mobile = mobdir + '/' + bn + suffix + '.pdb'

    # Check if mobile exists
    if not Path(f'{mobile}').exists():
        continue

    # Load structures
    cmd.load(target, 'target')
    cmd.load(mobile, 'mobile')

    # Align without outlier rejection (to get true value), then align normally
    aln = cmd.align('mobile and name CA', 'target and name CA', cycles=0)
    cmd.align('mobile and name CA', 'target and name CA', cycles=5)

    # Output
    outpdb = alndir + '/' + bn + suffix + '_aln.pdb'
    cmd.save(outpdb, 'mobile')

    # Append CA-RMSD to aligned structure
    with open(outpdb, 'a') as file:
        file.write(f'rmsd_ca {aln[0]:.4f}\n')

    print(f'Aligned {mobile} to {target}')
    cmd.delete('all')
EOF

# Output
for infile in ${infiles[@]} ; do
    bn=$(basename ${infile%.*})
    esmpdb="$esmdir/${bn}_esmfold.pdb"
    rlxpdb="$rlxdir/${bn}_esmfold_rlx.pdb"
    alnpdb="$alndir/${bn}_esmfold_rlx_aln.pdb"

    # Skip if ESMFold output does not exist
    [[ $(grep '^ATOM' $esmpdb 2>/dev/null | wc -l) == 0 ]] && {
        echo "ESMFold failed for $infile..."
        continue
    }

    # Grab footer
    pdblist=()
    for pdb in $rlxpdb $esmpdb $alnpdb ; do
        [[ $(grep "^ATOM" $pdb 2>/dev/null | wc -l) -gt 0 ]] && pdblist+=($pdb)
    done
    footer=$(awk '
        BEGINFILE { x = 0 }
        $1 ~ /^(ATOM|HETATM|TER|END|[ ]*)$/ { x = 1 ; next } x
    ' ${pdblist[@]})

    # If extension is .pdb, then use alnpdb over rlxpdb as output
    tmppdb=$(mktemp $tmpdir/tmp_XXXXXX.pdb)
    if [[ $infile =~ [.]pdb$ ]] ; then
        grep '^ATOM' $alnpdb 2>/dev/null > $tmppdb
    else
        grep '^ATOM' $rlxpdb 2>/dev/null > $tmppdb
    fi
    echo "$footer" >> $tmppdb

    # Output
    outpdb="$outdir/${bn}_esmfold.pdb"
    [[ $(grep "^ATOM" $tmppdb | wc -l) -gt 0 ]] && mv $tmppdb $outpdb
done
