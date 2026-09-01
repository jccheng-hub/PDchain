#!/usr/bin/env bash
usage () { cat << EOF
Usage: gen_invrots PDB CHRI:SS...
Generate inverse rotamers. Can scaffold AAs onto secondary structure.
Allows conformer generation for ligands as well.
Torsion angles are randomized, and each output will have a random suffix.

Parameters:
    PDB                     Input PDB file
    CHRI:SS                 Side chains to preserve in chain+resi format
                            E.g. A71 A73
                            Can specify number of secondary structure residues
                            on each side
                            E.g. A71:E1 A73:H2
                            ...will put A71 on a sheet (sheet-A71-sheet) and A73
                            on a helix (helix-helix-A73-helix-helix)
                            If you use "N" to specify secondary structure, H or
                            E will be randomly picked
                            E.g. A71:N2 A73:N3

Options:
  --keepres [str]           Residues to keep (chain+resi)
                            E.g. --keepres A71-73
  --numrots [int]           Number of inverse rotamers to generate
  --batchsize [int]         Number of inverse rotamers to attempt to generate at
                            a time (before pruning)
  --clashcut [int]          Distance cutoff for clash detection
  --parallel [int]          Number of jobs to spawn concurrently
  --inter_ctol [int]        Number of intermolecular clashes allowed
  --intra_ctol [int]        Number of intramolecular clashes allowed
  --torsfile [str]          Path to extra torsion text file
                            Format: tors["ABC"]["C1-C2-C3-C4"]="60,180,300"
                            ...where ABC is your three-letter-code, C1..C4 are
                            the four atoms that make up the dihedral, and a list
                            of comma-delimited angles.
                            Keywords can be used to replace the list of angles
                            E.g. tors["ABC"]["C1-C2-C3-C4"]=sp33
                            Angles for sp33: 60,180,300
                            Angles for sp23: 0,60,90,120,180,240,270,300
                            Angles for free: 0,5,10,15,...,355
  --keep_hydrogens          Keep hydrogens
  --outdir [str]            Path to output directory
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

val_opts=(numrots batchsize parallel keepres torsfile clashcut inter_ctol intra_ctol outdir)
bool_opts=(keep_hydrogens help)

numrots="1"
batchsize="1000"
parallel="1"
keepres=""
torsfile=""
clashcut="2.5"
inter_ctol="0"
intra_ctol="0"
outdir="./invrots"

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

shopt -s nullglob

# Get inputs
inpdb=${args[1]} ; bn=$(basename ${inpdb%.*})
res_str=${args[@]:2}
outdir=${outdir%/} ; mkdir -p $outdir
tmpdir=$(mktemp -d ${TMPDIR:-/tmp}/tmp_${USER}_XXXXXX)
trap 'kill $(jobs -p) 2>/dev/null; rm -rf "$tmpdir"' EXIT

# >>> gawk_torslib >>> {{{
[[ -n $torsfile ]] && {
    extra_tors=$(cat $torsfile)
    echo "Adding the following torsions:"
    echo "$extra_tors"
}
freerot=$(seq 0 5 355 | paste -sd ',')
gawk_torslib="$tmpdir/torslib.awk"
cat << EOF > $gawk_torslib
BEGIN {
    # Predfined angle sets
    sp33 = "60,180,300"
    sp23 = "0,60,90,120,180,240,270,300"
    free = "$freerot"

    # Allowed torsions for canonical amino acids (inverted rotamers)
    tors["ALA"]["1HB-CB-CA-N"]  = sp33  ; tors["ASN"]["OD1-CG-CB-CA"]  = sp23       
    tors["CYS"]["HG-SG-CB-CA"]  = free  ; tors["ASN"]["CG-CB-CA-N"]    = sp33       
    tors["CYS"]["SG-CB-CA-N"]   = sp33  ; tors["PRO"]                               
    tors["ASP"]["OD1-CG-CB-CA"] = sp23  ; tors["GLN"]["OE1-CD-CG-CB"]  = sp23       
    tors["ASP"]["CG-CB-CA-N"]   = sp33  ; tors["GLN"]["CD-CG-CB-CA"]   = sp33       
    tors["GLU"]["OE1-CD-CG-CB"] = sp23  ; tors["GLN"]["CG-CB-CA-N"]    = sp33       
    tors["GLU"]["CD-CG-CB-CA"]  = sp33  ; tors["ARG"]["NH1-CZ-NE-CD"]  = "0,180"    
    tors["GLU"]["CG-CB-CA-N"]   = sp33  ; tors["ARG"]["CZ-NE-CD-CG"]   = "180"      
    tors["PHE"]["CD1-CG-CB-CA"] = sp23  ; tors["ARG"]["NE-CD-CG-CB"]   = "180"      
    tors["PHE"]["CG-CB-CA-N"]   = sp33  ; tors["ARG"]["CD-CG-CB-N"]    = sp33       
    tors["GLY"]                         ; tors["SER"]["HG-OG-CB-CA"]   = free       
    tors["HIS"]["ND1-CG-CB-CA"] = sp23  ; tors["SER"]["OG-CB-CA-N"]    = sp33       
    tors["HIS"]["CG-CB-CA-N"]   = sp33  ; tors["THR"]["1HG-OG1-CB-CA"] = free       
    tors["ILE"]["CG1-CB-CA-N"]  = sp33  ; tors["THR"]["OG1-CB-CA-N"]   = sp33       
    tors["LYS"]["1HZ-NZ-CE-CD"] = free  ; tors["VAL"]["CG1-CB-CA-N"]   = sp33       
    tors["LYS"]["NZ-CE-CD-CG"]  = "180" ; tors["TRP"]["CD1-CG-CB-CA"]  = sp23       
    tors["LYS"]["CE-CD-CG-CB"]  = "180" ; tors["TRP"]["CD1-CG-CB-CA"]  = sp23       
    tors["LYS"]["CD-CG-CB-N"]   = "180" ; tors["TYR"]["HH-OH-CZ-CE1"]  = free       
    tors["LEU"]["CD1-CG-CB-CA"] = sp33  ; tors["TYR"]["CD1-CG-CB-CA"]  = sp23       
    tors["LEU"]["CG-CB-CA-N"]   = sp33  ; tors["TYR"]["CG-CB-CA-N"]    = sp33       
    tors["MET"]["1HE-CE-SD-CG"] = free
    tors["MET"]["CE-SD-CG-CB"]  = "180"

    # Additional torsions
    $extra_tors

    # one-letter-code
    olc["ALA"] = "A" ; olc["GLN"] = "Q" ; olc["LEU"] = "L" ; olc["SER"] = "S"
    olc["ARG"] = "R" ; olc["GLU"] = "E" ; olc["LYS"] = "K" ; olc["THR"] = "T"
    olc["ASN"] = "N" ; olc["GLY"] = "G" ; olc["MET"] = "M" ; olc["TRP"] = "W"
    olc["ASP"] = "D" ; olc["HIS"] = "H" ; olc["PHE"] = "F" ; olc["TYR"] = "Y"
    olc["CYS"] = "C" ; olc["ILE"] = "I" ; olc["PRO"] = "P" ; olc["VAL"] = "V"
}
EOF
# <<< gawk_torslib <<< }}}

# >>> write_py_ss_invrot () >>> {{{
write_py_ss_invrot () {
    gawk -f $gawk_torslib \
    -v keepres=$keepres -v res_str="$res_str" \
    -v batchsize="$batchsize" -v keepH="$keep_hydrogens" -v outdir="$1" \
    -v randseed=$(tr -cd '1-9' </dev/urandom | head -c8) -f - << 'EOF' $inpdb
    BEGIN {
        srand(randseed)
    
        # Atoms used for pair-fitting
        pf_atoms["ALA"]="N C CA CB"             ; pf_atoms["MET"]="CG SD CE"
        pf_atoms["CYS"]="CA CB SG"              ; pf_atoms["ASN"]="CB CG OD1 ND2"
        pf_atoms["ASP"]="CB CG OD1 OD2"         ; pf_atoms["PRO"]="N CA CB CD"
        pf_atoms["GLU"]="CG CD OE1 OE2"         ; pf_atoms["GLN"]="CG CD OE1 NE2"
        pf_atoms["PHE"]="CG CD1 CD2 CE1 CE2 CZ" ; pf_atoms["ARG"]="NE CZ NH1 NH2"
        pf_atoms["GLY"]="N CA C"                ; pf_atoms["SER"]="CA CB OG"
        pf_atoms["HIS"]="CG ND1 CD2 CE1 NE2"    ; pf_atoms["THR"]="CB OG1 CG2"
        pf_atoms["ILE"]="CB CG1 CD1"            ; pf_atoms["VAL"]="CB CG1 CG2"
        pf_atoms["LYS"]="CD CE NZ"              ; pf_atoms["TRP"]="CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"
        pf_atoms["LEU"]="CB CG CD1 CD2"         ; pf_atoms["TYR"]="CG CD1 CD2 CE1 CE2 CZ OH"
    
        # Split inputs (e.g. A71:E2 A73:H1)
        tot_res = split(res_str, a)
        for (i=1; i<=tot_res; i++) {
            split(a[i], b, ":")
            ch[i]     = substr(b[1], 1, 1)    # chain
            resi[i]   = substr(b[1], 2)       # residue number
            ss_id[i]  = substr(b[2], 1, 1)    # secondary structure (E, H, or N)
            ss_rep[i] = substr(b[2], 2, 1)    # number of ss stubs on each side
        }
    
        # Track chains
        i = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ; split(i, ch_fin, "")
    }
    
    $1~/^(ATOM|HETATM)$/ && res[$5][$6]=="" && ( $NF!="H" || keepH ) { resn[$5][$6] = $4 }
    
    ENDFILE {
        # Get basename
        n = split(FILENAME, a, "/") ; bn = a[n] ; sub(/\.pdb$/,"",bn)
    
        # Get random IDs
        i = rand_chars(batchsize) ; split(i, rand_id) ; id_idx=1
    
        # Get random SS
        i = rand_chars(batchsize * tot_res, 1, "12") ; split(i, rand_ss) ; ss_idx=1
    
        # Each loop iteration generates one inverse rotamer model
        for (r=1; r<=batchsize; r++) {
            ch_idx = 1 ; create = ""
            printf "cmd.create('fin_obj', 'none')\n"
    
            # Loop through each input residue
            for (t=1; t<=tot_res; t++) {
                i_resn = resn[ch[t]][resi[t]]
    
                # If residue does not have defined torsions (e.g. ligands), immediately transfer
                if ( !(i_resn in tors) ) {
                    printf "cmd.create('tmp_obj', '/oripdb//%s/%s')\n", ch[t], resi[t]
                    printf "cmd.alter('tmp_obj', 'chain=\"%s\"')\n", ch_fin[ch_idx++]
                    printf "cmd.create('fin_obj', 'fin_obj or tmp_obj')\n"
                    printf "cmd.delete('tmp_obj')\n"
                    continue
                }
    
                if (i_resn in olc) {
                    # Determine SS for input residue
                    if (ss_id[t] == "H")        { ss_i = 1
                    } else if (ss_id[t] == "E") { ss_i = 2
                    } else if (ss_id[t] == "N") { ss_i = rand_ss[ss_idx++]
                    }
    
                    # Make peptide string
                    pep = ""
                    for (i=1; i<=ss_rep[t]*2; i++) {
                        pep = pep "A" ; if (i==ss_rep[t]) pep = pep olc[i_resn]
                    }
                    pep = pep ? pep : olc[i_resn]
                    printf "cmd.fab('%s', 'tmp_obj',"         , pep
                    printf "chain='%s', resi='%s', ss='%s'"  , ch_fin[ch_idx++], resi[t] - ss_rep[t], ss_i
                    printf ")\n"
    
                    # Align
                    n = split(pf_atoms[i_resn], a)
                    nameplus = ""
                    for (i=1; i<=n; i++) {
                        nameplus = nameplus ? nameplus "+" a[i] : a[i]
                    }
                    printf "cmd.align("
                    printf "'/tmp_obj///%s and name %s'," , resi[t], nameplus
                    printf "'/oripdb//%s/%s and name %s'" , ch[t], resi[t], nameplus
                    printf ")\n"
                } else {
                    printf "cmd.create('tmp_obj', '/oripdb//%s/%s')\n", ch[t], resi[t]
                    printf "cmd.alter('tmp_obj', 'chain=\"%s\"')\n", ch_fin[ch_idx++]
                }
    
                # Generate rotamer
                for (key in tors[i_resn]) {
                    split(key, atoms, "-")
                    tmp = tors[i_resn][key]   ; n = split(tmp, angle, ",")
                    tmp = int(rand()*n) + 1   ; ang = angle[tmp]
                    printf "cmd.set_dihedral("
                    printf "'/tmp_obj///%s/%s'," , resi[t], atoms[1] 
                    printf "'/tmp_obj///%s/%s'," , resi[t], atoms[2]
                    printf "'/tmp_obj///%s/%s'," , resi[t], atoms[3]
                    printf "'/tmp_obj///%s/%s'," , resi[t], atoms[4]
                    printf "'%s')\n" , ang
                }
    
                printf "cmd.create('fin_obj', 'fin_obj or tmp_obj')\n"
                printf "cmd.delete('tmp_obj')\n"
            }
    
            # Transfer keepres
            if (keepres!="") {
                printf "cmd.create('tmp_obj', f'{keepres}')         \n"
                printf "cmd.alter('tmp_obj', 'chain=\"%s\"')        \n", ch_fin[ch_idx++]
                printf "cmd.create('fin_obj', 'fin_obj or tmp_obj') \n"
                printf "cmd.delete('tmp_obj')                       \n"
            }
    
            # Save
            outpdb = outdir "/" bn "_" rand_id[id_idx++] ".pdb"
            if (!keepH) {
                printf "cmd.remove('elem H')                                         \n", outpdb
            }
            printf "cmd.save('%s', 'fin_obj')                                    \n", outpdb
            printf "#print(f'Output saved at %s (%s/%s)            \\r', end='')  \n", outpdb, r, batchsize
            printf "cmd.delete('fin_obj')                                        \n"
        }
    }
    
    function rand_chars(num,    len, chars, clen, i, j, output) {
        len = len ? len : 6
        chars = chars ? chars : "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        clen = length(chars)
        for (i=1; i<=num; i++) {
            for (j=1; j<=len; j++) {
                output = output substr(chars, int(rand() * clen)+1, 1)
            }
            if (i < num) output = output " "
        }
        return output
    }
EOF
}
# <<< write_py_ss_invrot () <<< }}}

# >>> detect_chain_clashes () >>> {{{
detect_chain_clashes () {
    local args=(${@:2})
    gawk -v clashcut="$1" -v tot=${#args[@]} -f - << 'EOF' ${@:2}
    BEGIN { clashcut_sq = clashcut * clashcut }
    BEGINFILE { inter_clashes = 0 ; intra_clashes = 0 }
    
    $1~/^(ATOM|HETATM)$/ && $NF!="H" {
        x[$5][$6][$3] = $7 ; y[$5][$6][$3] = $8 ; z[$5][$6][$3] = $9
    }
    
    ENDFILE {
        for (ch1 in x) for (ch2 in x) {
            if (ch1 > ch2) continue    # Avoid double-counting chains
            for (ri1 in x[ch1]) for (ri2 in x[ch2]) {
                for (nm1 in x[ch1][ri1]) for (nm2 in x[ch2][ri2]) {
    
                    # Calculate distance
                    dx = x[ch1][ri1][nm1] - x[ch2][ri2][nm2]
                    dy = y[ch1][ri1][nm1] - y[ch2][ri2][nm2]
                    dz = z[ch1][ri1][nm1] - z[ch2][ri2][nm2]
                    d_sq = dx*dx + dy*dy + dz*dz
    
                    if (d_sq > clashcut_sq) continue
                    if (ch1 != ch2) { inter_clashes++ } else { intra_clashes++ }
                }
            }
        }
        id++
        printf "inter_clashes %s %s\n", inter_clashes, clashcut >> FILENAME
        printf "intra_clashes %s %s\n", intra_clashes, clashcut >> FILENAME
        # printf "%s intermolecular clashes found in %s (%s/%s)         \r", inter_clashes, FILENAME, id, tot
        # printf "%s intramolecular clashes found in %s (%s/%s)         \r", intra_clashes, FILENAME, id, tot
    
        # Clean up
        delete x ; delete y ; delete z
    }
EOF
}
# <<< detect_chain_clashes () <<< }}}

# Check input
[[ -s $inpdb ]] || {
    echo "Input PDB $inpdb does not exist" && exit
}

# >>> keepres to keepstr >>> {{{
[[ -n $keepres ]] && keepstr=$(gawk -v keepres="$keepres" '
    BEGIN {
        keepexp = expand_range(keepres)
        n = split(keepexp, a)
        for (i=1; i<=n; i++) {
            b = substr(a[i], 1, 1) "/" substr(a[i], 2) "/"
            keepstr = keepstr ? keepstr " + " b : b
        }
        print keepstr
    }
    
    function expand_range(res,    a, b, i, j, n, m, ch, range) {
        n = split(res, a)
        for (i=1; i<=n; i++) {
            ch = a[i]; gsub(/[^A-Z]/, "", ch)
            gsub(/[A-Z]/, "", a[i])
            m = split(a[i], b, "-")
            for (j=b[1]; j<=b[m]; j++) {
                range = range == "" ? ch j : range " " ch j
            }
        }
        return range
    }
')
# <<< keepres to keepstr <<< }}}

# >>> gen_invrots () >>> {{{
gen_invrots () {
    local batchdir="$1" 
    local groupdir="$2"
    mkdir -p $batchdir $groupdir
    local bnsub=$(basename $batchdir)

    # Generate inverse rotamers
    echo "[$bnsub] Generating $batchsize inverse rotamers for $inpdb..."
    local pyscript=$(mktemp $batchdir/tmp_XXXXXX.py)
    sed -E 's/^ {4}//' << EOF > $pyscript
    from pymol import cmd
    keepres = 'oripdb and (${keepstr})'
    cmd.load('$inpdb', 'oripdb')
EOF
    write_py_ss_invrot $batchdir >> $pyscript
    sed -E 's/^ {4}//' << EOF >> $pyscript
    print(f'[$bnsub] Finished generating $batchsize inverse rotamers.')
EOF
    local outpdbs=($(gawk '/Output saved at/ {print$4}' $pyscript))
    python $pyscript && rm $pyscript

    # Run clash checker
    echo "[$bnsub] Checking for clashes: < $clashcut angstroms (--clashcut $clashcut)..."
    detect_chain_clashes "$clashcut" ${outpdbs[@]}
    
    # Get min intra clash value
    local min_intra=$(grep -h "^intra_clashes " ${outpdbs[@]} | sort -k2n | awk '{print$2;exit}')
    echo "[$bnsub] Reference intramolecular clash value from this batch: $min_intra"
    
    # Remove clashes
    echo "[$bnsub] Number of intermolecular clashes allowed: $inter_ctol (--inter_ctol $inter_ctol)"
    echo "[$bnsub] Number of intramolecular clashes allowed: $min_intra + $intra_ctol (--inter_ctol $intra_ctol)"
    echo "[$bnsub] Pruning inverse rotamers..."
    local inter=($(grep -H "^inter_clashes " ${outpdbs[@]} | sed 's/ /:/g'))
    local intra=($(grep -H "^intra_clashes " ${outpdbs[@]} | sed 's/ /:/g'))
    local tormv=($(gawk -v inter_ctol="$inter_ctol" -v intra_ctol="$intra_ctol" -v min_intra="$min_intra" '
        $1 == "inter_clashes" && $2 > inter_ctol { tormv[FILENAME] = 1 }
        $1 == "intra_clashes" && $2 > intra_ctol + min_intra { tormv[FILENAME] = 1 }
        END { for (i in tormv) print i }
    ' ${outpdbs[@]}))
    rm ${tormv[@]} &>/dev/null && echo "[$bnsub] Removed ${#tormv[@]} inverse rotamers from $batchdir"

    local outpdbs=($batchdir/*.pdb)
    [[ ${#outpdbs[@]} -gt 0 ]] &&
    mv ${outpdbs[@]} $groupdir && echo "[$bnsub] ${#outpdbs[@]} rotamers saved to $groupdir"
}
# <<< gen_invrots () <<< }}}

# Cycle
cycle_gen_invrots () {
    local batchdir=$1
    local groupdir=$2
    local bnsub=$(basename $1)
    local outrots=($groupdir/$(basename ${inpdb%.*})_??????.pdb)
    while [[ ${#outrots[@]} -lt $numrots ]] ; do
        gen_invrots $batchdir $groupdir
        local outrots=($groupdir/$(basename ${inpdb%.*})_??????.pdb)

        echo "[$bnsub] Pruning redundant rotamers..."
        local tormv=($(sha256sum ${outrots[@]} 2>/dev/null | awk 'h[$1]++ {print$2}'))
        rm ${tormv[@]} &>/dev/null && echo "[$bnsub] Removed ${#tormv[@]} inverse rotamers from $groupdir"

        echo "[$bnsub] ${#outrots[@]} total rotamers in $groupdir"
    done
}

# Parallel
parallel=$(awk -v n=$(nproc) -v p="$parallel" 'BEGIN { p = p > n ? n : int(p) ; print p}')
echo "Spawning $parallel jobs to generate $numrots rotamers"
echo "Batchsize before filtering: $batchsize"
for i in $(seq 1 $parallel) ; do
    cycle_gen_invrots $tmpdir/batch$i $tmpdir/group &
done
wait

# Move to outputs
mkdir -p $outdir
outpdbs=($tmpdir/group/*.pdb)
for pdb in ${outpdbs[@]} ; do
    mv $pdb $outdir 
done && echo "Saved ${#outpdbs[@]} rotamers to $outdir"

runtime=$(awk "BEGIN{print($SECONDS/60)}")
echo "gen_invrots runtime: $runtime minutes"
