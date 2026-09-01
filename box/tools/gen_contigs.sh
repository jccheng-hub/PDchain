#!/usr/bin/env bash
usage () { cat << EOF
Usage: gen_contigs [RES...]
Generate specific contigs for RFdiffusion

E.g. gen_contigs A3-3 A7-10 A72-72
...could result in the following contigs
23-23/A3-3/11-11/A7-10/39-39/A72-72/57-57

Depending on the input, the default script might fail to find a contigs string
that matches the specifications. Tune the parameters --mingap, --mintot,
--addgap, and --addtot to increase likelihood of success.

Parameters:
    RES                     Residues to keep in chain+resi format
                            E.g. A9-9
                            E.g. A5-10

Options for default usage:
  --mingap [int]            Mininum residues in each gap
  --addgap [int]            Maximum possible number added to --mingap
                            E.g. --mingap 20 --addgap 60 allows each gap to be
                            between 20-80 residues long.
  --mintot [int]            Minimum residues in diffused protein
  --addtot [int]            Maximum possible rumber added to --mintot
                            E.g. --mintot 160 --addtot 40 allows the total
                            residue count to be 160-200
  --random_order            Randomize the input RES order

Options for ss usage:
  --ss_to_contigs           Uses input pdb's secondary structure to construct
                            contigs. Secondary structure is determined using
                            the DISICL algorithm (bb dihedral-only approach).
                            Residues on secondary structures are preserved in
                            addition to the RES strings provided.
  --inpdb [str]             Input pdb for determining contigs
  --ss_trim [int]           If using ss_to_contigs, then trims secondary
                            structures based on their sequence distance to
                            loops.
                            E.g. --ss_trim 1 with an ss string of LLHHHHLL
                            allows the two H adjacent L to be rediffused. This
                            effectively treats the ss string as LLLHHLLL. If
                            --ss_trim 2 was used instead, then the ss string is
                            effectively LLLLLLLL.
                            Can specify with range (e.g. --ss_term 2-4)
  --nterm_trim [int]        Trim specified number of residues from N terminus.
                            Can specify with range (e.g. --nterm_trim 2-4)
  --cterm_trim [int]        Trim specified number of residues from C terminus.
                            Can specify with range (e.g. --cterm_trim 2-4)
  --helix_cap [int]         Maximum length of a helix before making it entirely
                            rediffusable
  --reset_perc [float]      Percent chance for setting ss_trim to 100, which
                            will make all non-fixed residues rediffusable
  --regap                   Redetermine gaps based on default usage.
                            Preserves all fixed backbones after --ss_trim [int]

General options:
  --vary_linkers [int]      Allow linkers to vary by specified INT
  --nterm_add [int]         Number of residues to add to N-terminus
                            Can specify with range (e.g. --nterm_add 4-8)
  --cterm_add [int]         Number of residues to add to C-terminus
                            Can specify with range (e.g. --cterm_add 4-8)
  --inc_only                When used with --vary_linkers, only allow the linker
                            length to increase
  --dec_only                When used with --vary_linkers, only allow the linker
                            length to decrease
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

val_opts=(
    mingap              addgap              mintot              addtot
    inpdb               ss_trim             nterm_trim          cterm_trim
    vary_linkers        nterm_add           cterm_add           helix_cap
    reset_perc 
)

bool_opts=(
    random_order        ss_to_contigs       dec_only            inc_only
    regap               help
)

mintot="180"
addtot="40"
mingap="20"
addgap="100"
inpdb=""
ss_trim="1"
nterm_trim="1"
cterm_trim="1"
nterm_add="0"
cterm_add="0"
vary_linkers="0"
helix_cap=""
reset_perc="0"

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

# Collapse inputs with ranges
collapse_range () {
    awk -v range="$1" 'BEGIN {
        n = split(range, a, "-")
        a[2] = a[2] ? a[2] : a[1]
        for (i=a[1]+0 ; i<=a[2]+0; i++) print i
    }' | shuf -n1
}

ss_trim=$(collapse_range $ss_trim)
nterm_add=$(collapse_range $nterm_add)
cterm_add=$(collapse_range $cterm_add)
nterm_trim=$(collapse_range $nterm_trim)
cterm_trim=$(collapse_range $cterm_trim)

# If reset_perc is enabled, max out ss_trim at a specified chance
[[ $(echo "$reset_perc > 0" | bc) == 1 ]] && {
    rand=$(( (RANDOM % 100) + 1 ))
    [[ $(echo "$rand <= $reset_perc" | bc) == 1 ]] && ss_trim=100
}

# >>> dss_disicl_tab () >>> {{{
dss_disicl_tab () {
    gawk -f - << 'EOF' $@
    BEGIN {
        # DISICL was made by Nagy and Oostenbrink, J. Chem. Inf. Model. 2013
        # This is my attempt to implement an GNU awk version of DISICL
    
        # Define min and max angles for phi and psi
        minphi["a1"][1] = -95    ; minphi["a2"][1] = -107   ; minphi["b1"][1] = -135 
        maxphi["a1"][1] = -40    ; maxphi["a2"][1] = -40    ; maxphi["b1"][1] = -87  
        minpsi["a1"][1] = -70    ; minpsi["a2"][1] = -32    ; minpsi["b1"][1] = 95   
        maxpsi["a1"][1] = -32    ; maxpsi["a2"][1] = -12    ; maxpsi["b1"][1] = 150  
                                         
        minphi["b2"][1] = -175   ; minphi["b2"][2] = -180   ; minphi["b2"][3] = -135
        maxphi["b2"][1] = -135   ; maxphi["b2"][2] = -135   ; maxphi["b2"][3] = -105
        minpsi["b2"][1] = 95     ; minpsi["b2"][2] = 136    ; minpsi["b2"][3] = 150 
        maxpsi["b2"][1] = 136    ; maxpsi["b2"][2] = 180    ; maxpsi["b2"][3] = 180 
    
        minphi["b2"][4] = -180   ; minphi["d"][1] = -150    ; minphi["d"][2]  = -150
        maxphi["b2"][4] = -135   ; maxphi["d"][1] = -67     ; maxphi["d"][2]  = -107
        minpsi["b2"][4] = -180.1 ; minpsi["d"][1] = 8       ; minpsi["d"][2]  = -32 
        maxpsi["b2"][4] = -160   ; maxpsi["d"][1] = 40      ; maxpsi["d"][2]  = 8   
                                         
        minphi["d2"][1] = -165   ; minphi["d1"][1] = -107   ; minphi["dx"][1] = 38 
        maxphi["d2"][1] = -95    ; maxphi["d1"][1] = -40    ; maxphi["dx"][1] = 140
        minpsi["d2"][1] = -70    ; minpsi["d1"][1] = -12    ; minpsi["dx"][1] = -25
        maxpsi["d2"][1] = -32    ; maxpsi["d1"][1] = 8      ; maxpsi["dx"][1] = 75 
                                         
        minphi["g"][1]  = 55     ; minphi["gx"][1] = -113   ; minphi["p"][1]  = -87
        maxphi["g"][1]  = 95     ; maxphi["gx"][1] = -60    ; maxphi["p"][1]  = -30
        minpsi["g"][1]  = -100   ; minpsi["gx"][1] = 50     ; minpsi["p"][1]  = 95 
        maxpsi["g"][1]  = -40    ; maxpsi["gx"][1] = 95     ; maxpsi["p"][1]  = 180
                                         
        minphi["p"][2] = -100    ; minphi["px"][1] = 35     ; minphi["px"][2] = 60  
        maxphi["p"][2] = -60     ; maxphi["px"][1] = 100    ; maxphi["px"][2] = 120 
        minpsi["p"][2] = -180.1  ; minpsi["px"][1] = -180.1 ; minpsi["px"][2] = 150 
        maxpsi["p"][2] = -150    ; maxpsi["px"][1] = -110   ; maxpsi["px"][2] = 180 
                                         
        minphi["z"][1]  = -170
        maxphi["z"][1]  = -113
        minpsi["z"][1]  = 50
        maxpsi["z"][1]  = 95
    
        # Define secondary structure classification
        ssc["H"]["a1-a1"] = "alpha-helix"
        ssc["H"]["a1-a2"] = "alpha-helix"
        ssc["H"]["a2-a1"] = "alpha-helix"
        ssc["H"]["a1-d2"] = "pi-helix"
        ssc["H"]["d2-d2"] = "pi-helix"
        ssc["H"]["d2-a1"] = "pi-helix"
        ssc["H"]["a2-d2"] = "pi-helix"
        ssc["H"]["a2-d2"] = "helix-cap"
        ssc["H"]["d-d2"]  = "helix-cap"
        ssc["H"]["d1-d2"] = "helix-cap"
        ssc["H"]["d2-d"]  = "helix-cap"
        ssc["H"]["b1-a1"] = "helix-cap"
        ssc["H"]["b2-a1"] = "helix-cap"
        ssc["H"]["b2-a2"] = "helix-cap"
        ssc["H"]["p-a1"]  = "helix-cap"
        ssc["H"]["p-a2"]  = "helix-cap"
        ssc["E"]["b2-b2"] = "extended-beta-strand"
        ssc["E"]["b1-b1"] = "normal-beta-strand"
        ssc["E"]["b1-b2"] = "normal-beta-strand"
        ssc["E"]["b2-b1"] = "normal-beta-strand"
        ssc["E"]["b1-p"]  = "beta-cap"
        ssc["E"]["b2-p"]  = "beta-cap"
        ssc["E"]["p-b1"]  = "beta-cap"
        ssc["E"]["p-b2"]  = "beta-cap"
    
        # Iterate through arrays in ascending numerical order
        PROCINFO["sorted_in"] = "@ind_num_asc"
    
        # Define SS priority
        ssp["H"] = 3 ; ssp["E"] = 2 ; ssp["L"] = 1
    
        # Define pi
        PI = atan2(0, -1)
    }
    
    $1 ~ /^(ATOM|HETATM)$/ {
        tlc[$5][$6] = $4                     # residue name (three-letter code)
        elm[$5][$6][$3] = substr($0, 78, 1)  # element
        bfa[$5][$6][$3] = $11                # b-factor / pLDDT
        crd[$5][$6][$3][1] = $7              # x coordinate
        crd[$5][$6][$3][2] = $8              # y coordinate
        crd[$5][$6][$3][3] = $9              # z coordinate
    }
    
    ENDFILE {
        # Get phi and psi angles, then identify DISCL region
        for (chain in crd) for (resi in crd[chain]) {
    
            # Cannot be N or C terminal residues
            if (!(resi-1 in crd[chain] && resi+1 in crd[chain])) continue
    
            # Must have C at N terminus and N at C terminus
            if (!("N"  in crd[chain][resi] &&
                  "CA" in crd[chain][resi] &&
                  "C"  in crd[chain][resi] &&
                  "C"  in crd[chain][resi-1] &&
                  "N"  in crd[chain][resi+1] )) { continue }
    
            # Calculate phi and psi
            phi = calc_dihedral(crd[chain][resi-1]["C"], crd[chain][resi]["N"], crd[chain][resi]["CA"], crd[chain][resi]["C"])
            psi = calc_dihedral(crd[chain][resi]["N"], crd[chain][resi]["CA"], crd[chain][resi]["C"], crd[chain][resi+1]["N"])
    
            # Identify region using phi and psi
            breaksig = 0
            for (reg_i in minphi) {
                for (i in minphi[reg_i]) {
                    if ( phi >= minphi[reg_i][i] &&
                         phi <= maxphi[reg_i][i] &&
                         psi >= minpsi[reg_i][i] &&
                         psi <= maxpsi[reg_i][i] )
                    {
                        reg[chain][resi] = reg_i ; breaksig = 1 ; break
                    }
                }
                if (breaksig) break
            }
        }
    
        # Assign secondary structure
        for (chain in crd) for (resi in crd[chain]) {
    
            # Assign first and last residue of each chain as loops
            if (!(resi-1 in crd[chain] && resi+1 in crd[chain])) {
                if (ss[chain][resi] == "") ss[chain][resi] = "L"
                continue
            }
    
            # Skip assignment of resi + 2 doesn't exist (should already be assigned)
            if (!(resi+2 in crd[chain])) continue
    
            # Get segment type
            reg_i = reg[chain][resi]
            reg_j = reg[chain][resi+1]
            seg_type = "L"
            for (N in ssc) {
                if (reg_i "-" reg_j in ssc[N]) {
                    seg_type = N ; break
                }
            }
    
            # Apply seg_type based on priority
            for (i=0; i<=1; i++) {
                if (ssp[seg_type] > ssp[ss[chain][resi+i]]) {
                    ss[chain][resi+i] = seg_type
                }
            }
        }
    
        for (chain in crd) for (resi in crd[chain]) {
            if ( !("CA" in crd[chain][resi] && "C" in elm[chain][resi]) ) continue
            printf "%s %s %s\n", chain, resi, ss[chain][resi]
        }
    
        delete tlc ; delete bfa ; delete crd ; delete reg ; delete ss
    }
    
    function check_coordinates(crd,    ch, ri, at, x, y, z) {
        # Print out contents of coordinates multi-dimensional array
        # crd: coordinates array | crd[chain][resi][atom][1/2/3] = x/y/z
        for (ch in crd) for (ri in crd[ch]) for (at in crd[ch][ri]) {
            x = crd[ch][ri][at][1]
            y = crd[ch][ri][at][2]
            z = crd[ch][ri][at][3]
            printf "%-8s%-4s%-16s%-16s%-16s\n", ch ri, at, x, y, z
        }
    }
    
    function calc_dist(v1, v2,    dx, dy, dz) {
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        dx = v1[1]-v2[1] ; dy = v1[2]-v2[2] ; dz = v1[3]-v2[3]
        return sqrt( dx*dx + dy*dy + dz*dz )
    }
    
    function calc_dist_sq(v1, v2,    dx, dy, dz) {
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        dx = v1[1]-v2[1] ; dy = v1[2]-v2[2] ; dz = v1[3]-v2[3]
        return dx*dx + dy*dy + dz*dz
    }
    
    function calc_dihedral(v1, v2, v3, v4,    b1, b2, b3, n1, n2, b2n, m, x, y, rad) {
        # Calculate dihedral angle of v1-v2-v3-v4
        # v1/v2/v3/v4: input array | v[1/2/3] = x/y/z
        vec_sub(v2, v1, b1)
        vec_sub(v3, v2, b2)
        vec_sub(v4, v3, b3) 
        vec_cross(b1, b2, n1)
        vec_cross(b2, b3, n2)
        vec_scale(b2, 1 / vec_mag(b2), b2n)
        vec_cross(b2n, n1, m)
        x = vec_dot(n1, n2)
        y = vec_dot(m, n2)
        rad = atan2(y, x)
        return rad * 180 / PI
    }
    
    function vec_add(v1, v2, vf) {
        # Add two 3D vectors
        # v1: input array 1        | v1[1/2/3] = x/y/z
        # v2: input array 2        | v2[1/2/3] = x/y/z
        # vf: name of output array | vf[1/2/3] = x/y/z
        delete vf ; vf[1] = v1[1]+v2[1] ; vf[2] = v1[2]+v2[2] ; vf[3] = v1[3]+v2[3]
    }
    
    function vec_sub(v1, v2, vf) {
        # Subtract 3D vectors
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        # vf: output array  | vf[1/2/3] = x/y/z
        delete vf ; vf[1] = v1[1]-v2[1] ; vf[2] = v1[2]-v2[2] ; vf[3] = v1[3]-v2[3]
    }
    
    function vec_cross(v1, v2, vf) {
        # Cross product of two 3D vectors
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        # vf: output array  | vf[1/2/3] = x/y/z
        delete vf
        vf[1] = v1[2]*v2[3] - v1[3]*v2[2] 
        vf[2] = v1[3]*v2[1] - v1[1]*v2[3]
        vf[3] = v1[1]*v2[2] - v1[2]*v2[1]
    }
    
    function vec_dot(v1, v2) {
        # Dot product of two 3D vectors
        # v1: input array 1 | v1[1/2/3] = x/y/z
        # v2: input array 2 | v2[1/2/3] = x/y/z
        return v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3]
    }
    
    function vec_scale(v, scalar, vf) {
        # Scale a 3D vector
        # v: input array       | v[1/2/3] = x/y/z
        # scalar: scalar value | int
        # vf: output array     | vf[1/2/3] = x/y/z
        delete vf
        vf[1] = scalar * v[1]
        vf[2] = scalar * v[2]
        vf[3] = scalar * v[3]
    }
    
    function vec_mag(v,    a_sq, b_sq, c_sq) {
        # Get magnitude of 3D vector
        # v: input array | v[1/2/3] = x/y/z
        a_sq = v[1] * v[1]
        b_sq = v[2] * v[2]
        c_sq = v[3] * v[3]
        return sqrt( a_sq + b_sq + c_sq )
    }
    
    function vec_norm(v, vf,    mag) {
        # Normalize vector
        # v: input array   | v[1/2/3] = x/y/z
        # vf: output array | vf[1/2/3] = x/y/z
        delete vf ; mag = vec_mag(v) ; vec_scale(v, 1/mag, vf)
    }
EOF
}
# <<< dss_disicl_tab () <<< }}}

gen_contigs_unconditional () {
    local num=$(seq $mintot $((mintot + addtot)) | shuf -n1)
    echo $num-$num
}

# >>> gen_contigs_from_reslist () >>> {{{
gen_contigs_from_reslist () {
    local reslist=($@)
    local numfixed=$(sed 's/ /\n/g' <<< "${reslist[@]}" | awk '
        {
            split($0, a, "-") ; sub(/^[A-Z]/, "", a[1])
            if ( a[2] == "" ) a[2] = a[1]
            for (i=int(a[1]) ; i<=int(a[2]) ; i++) j++
        }

        END { print j }
    ')
    [[ $random_order == 1 ]] && reslist=($(shuf -e ${reslist[@]}))
    [[ $numfixed -gt $mintot ]] && mintot=$numfixed
    mintot=$((mintot - numfixed))
    sed 's/ /\n/g' <<< "${reslist[@]}" | awk -v mingap="$mingap" -v addgap="$addgap" '
    BEGIN {
        srand()
        while (!(sum > '$mintot' && sum < '$mintot' + '$addtot')) {
            sum = 0
            for (i=1; i<='${#reslist[@]}'+1; i++) {
                g[i] = int(mingap + rand() * (addgap + 1))
                sum += g[i]
            }
            if (sum < mintot) addgap += 1
            if (sum > mintot) addgap = (addgap > 0) ? addgap - 1 : addgap
            tries++ ; if (tries > 10000) { print "FAILED" ; exit}
        }
        i = 0
    }

    {
        split($0, a, "-"); ch = substr(a[1], 1, 1); sub(/^[A-Z]/, "",a[1])
        if (a[2] == "") a[2]=a[1]
        printf g[++i] "-" g[i] "/" ch a[1] "-" a[2] "/"
    }

    END {
        if (tries > 10000) exit
        print g[++i]"-"g[i]
    }
    '
}
# <<< gen_contigs_from_reslist () <<< }}}

# >>> gen_contigs_from_disicl () >>> {{{
gen_contigs_from_disicl () {
    local inpdb=$1
    local reslist=${@:2}
    local restab=$(dss_disicl_tab $inpdb)
    echo "$restab" |
    gawk -v reslist="$reslist" -v ss_trim="$ss_trim" \
         -v nterm_trim="$nterm_trim" -v cterm_trim="$cterm_trim" \
         -v helix_cap="$helix_cap" '
    BEGIN {
        n = split(reslist, a, " ")
        for (i=1; i<=n; i++) {
            ch_i = substr(a[i], 1, 1) ; sub(/^[A-Z]/, "", a[i])
            split(a[i], b, "-") ; if (!b[2]) b[2] = b[1]
            for (j=b[1]; j<=b[2]; j++) reskey[ch_i j] = 1
        }
    }
    ch[NR-1] && $1!=ch[NR-1] { ss_str = ss_str " " }
    { ch[NR]=$1 ; ri[NR]=$2 ; ss[NR]=$3 ; ss_str = ss_str ss[NR] }
    $1 $2 in reskey { key[NR] = 1 }
    END {

        # Make sure ss structures are not adjacent
        gsub(/EH/, "LL", ss_str)
        gsub(/HE/, "LL", ss_str)

        # Cap helix length
        if (helix_cap) {
            n = split(ss_str, a, "")
            i=1
            while (i<=n) {
                j = i ; while (a[j] == a[i]) j++

                rep = j - i

                if (a[i] == "H" && rep > helix_cap + 0) {
                    insert = "" ; for (k=1; k<=rep; k++) { insert = insert "L" }
                    ss_str = substr(ss_str, 1, i - 1) insert substr(ss_str, j)
                }

                i = j
            }
        }

        # Trim ss_str string
        for (i=1; i<=ss_trim; i++) {
            gsub(/[^L ]L/, "LL", ss_str) ; gsub(/L[^L ]/, "LL", ss_str)
        }

        # Trim N and C termini
        nterm_trim = nterm_trim + 0
        cterm_trim = cterm_trim + 0

        # N-term
        if (nterm_trim > 1) {
            term = "" ; for (i=1; i<=nterm_trim; i++) term = term "L"
            gsub("^.{" nterm_trim "}", term, ss_str)
            gsub(" .{" nterm_trim "}", " " term, ss_str)

        }

        # C-term
        if (cterm_trim > 1) {
            term = "" ; for (i=1; i<=cterm_trim; i++) term = term "L"
            gsub(".{" cterm_trim "}$", term, ss_str)
            gsub(".{" cterm_trim "} ", term " ", ss_str)
        }

        # Slice in key residues
        gsub(" ", "", ss_str)
        gsub(/[^L]/, "X", ss_str)
        for (i=1; i<=NR; i++) {
            if (i in key) {
                ss_str = substr(ss_str, 1, i-1) "X" substr(ss_str, i+1)
            }
        }

        # Remake ss array from ss_str
        delete ss ; split(ss_str, ss, "")
    
        # Make contigs
        c = 1    # index for contigs
        for (i=1; i<=NR; i++) {
            if (ch[i-1] && ch[i-1]!=ch[i]) { ctg[c] = ctg[c] "/0" ; c++ ; count = 0 }
            if (ss[i-1] && ss[i-1]!=ss[i]) { count = 0 } ; count++
            if (ss[i]=="L") {
                block = count "-" count
                if (!ss[i+1] || ss[i+1] != "L" || (ch[i+1] && ch[i+1] != ch[i]) ) {
                    ctg[c] = ctg[c] ? ctg[c] "/" block : block
                }
            } else {
                block = ch[i] (ri[i] - count + 1) "-" ri[i]
                if (!ss[i+1] || ss[i+1] == "L") {
                    ctg[c] = ctg[c] ? ctg[c] "/" block : block
                }
            }
        }
        for (i=1; i<=c; i++) contigs = contigs ? contigs " " ctg[i] : ctg[i]
        print contigs
    }
    '
}
# <<< gen_contigs_from_disicl () <<< }}}

edit_linkers () {
    local contigs=$@
    echo $contigs | sed 's|/|\n|g' | sed 's/ /\n/g' |
    awk -v dev=$vary_linkers -v inc=$inc_only -v dec=$dec_only '
    BEGIN { srand() }
    $0~/[A-Z]/ || $0 == "0" {print $0; next} {
        delete a; split($0, a, "-")
        d = 0 + int(rand()*(dev-0+1)) ; coin = rand()
        if (inc)              { n = a[1] + d }
        else if (dec)         { n = a[1] - d } 
        else if (coin < 0.5)  { n = a[1] + d }
        else if (coin >= 0.5) { n = a[1] - d }
        n = (n > 1) ? n : 1
        if (rand() < 0.75) { print $0 } else { print n"-"n }
    }' | paste -sd '/' | sed 's|/0/|/0 |'
}

reslist=${args[@]:1}

# Generate initial contigs
if [[ $ss_to_contigs == 1 ]] ; then
    [[ -z $inpdb ]] && { echo "Requires --inpdb if invoking --ss_to_contigs"; exit ;}
    contigs=$(gen_contigs_from_disicl $inpdb $reslist)
elif [[ -n $reslist ]] ; then
    contigs=$(gen_contigs_from_reslist $reslist)
else
    contigs=$(gen_contigs_unconditional)
fi

# Vary linkers
contigs=$(edit_linkers $contigs)

# Regap
[[ $ss_to_contigs == 1 && $regap == 1 ]] && {
    ss_reslist=$(echo $contigs | sed 's|/|\n|g' | grep ^[A-Z] | paste -sd ' ')
    contigs=$(gen_contigs_from_reslist $ss_reslist)
}

# Add to N and C terms
[[ $nterm_add -gt 0 || $cterm_add -gt 0 ]] && {
    
    # Check total residue count
    totres=$(sed 's/ /\n/g' <<< "$contigs" | sed 's|/|\n|g' | awk '
        $1 ~ /^[A-Z]/ {
            gsub(/[A-Z]/, "", $0)
            split($0, a, "-")
            a[2] = a[2] ? a[2] : a[1]
            for (i=a[1]+0; i<=a[2]+0; i++) { totres++ }
            next
        }

        $1 !~ /^[A-Z]/ { split($0, a, "-") ; totres += a[1] }

        END { print totres }
    ')

    # Only add residues if the final residue count does not exceed maximum
    [[ $((totres + nterm_add + cterm_add)) -lt $((mintot + addtot)) ]] && {
        contigs=$(
            sed 's/ /\n/g' <<< "$contigs" |
            awk -v nterm_add=$nterm_add -v cterm_add=$cterm_add -v OFS='/' -F '/' '
                $1 !~ /^[A-Z]/ {
                    n = split($1, a, "-")
                    $1 = a[1]+nterm_add "-" a[1]+nterm_add
                }

                $NF !~ /^[A-Z]/ && $NF != 0 {
                    n = split($NF, a, "-")
                    $NF = a[1]+cterm_add "-" a[1]+cterm_add
                }

                $(NF-1) !~ /^[A-Z]/ && $NF == 0 {
                    n = split($(NF-1), a, "-")
                    $(NF-1) = a[1]+cterm_add "-" a[1]+cterm_add
                }

                { print }
            ' | paste -sd ' '
        )
    }
}

# Output
echo $contigs
