#!/usr/bin/env bash
usage () { cat << EOF
Usage: dss_discl PDB [PDB...]
Uses DISICL algorithm to assign secondary structure to each residue

Parameters:
    PDB                     Input PDB

Options:
  --append                  Add secondary structure info directly to input PDB
                            Will have the following appended:
                            ss_str, perc_loop, perc_helix, perc_sheet
  --fixedres [str]          Fixed residues in chain+resi format. If specified,
                            will return the secondary structure composition for
                            the fixed residues +/- the value for --fixedrange.
                            For example, --fixedres A7 A11 --fixedrange 1 will
                            return the secondary structure composition for
                            A6 A7 A8 A10 A11 A12.
  --fixedrange [int]        Number of residues that flank each fixed residue
                            specified by --fixedres. See --fixedres.
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

val_opts=(fixedres fixedrange)
bool_opts=(append help)
fixedres=""
fixedrange="2"

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

# >>> dss_disicl () >>> {{{
dss_disicl () {

    # Check append flag
    if [[ $append == 1 ]] ; then
        local inplace_flag="-i inplace"
    else
        local inplace_flag=""
    fi

    gawk ${inplace_flag} -v fixedres="$fixedres" -v fixedrange="$fixedrange" \
    -f - << 'EOF' $@
    BEGIN {
        # DISICL was made by Nagy and Oostenbrink, J. Chem. Inf. Model. 2013
        # This is my attempt to reimplement a GNU awk version of DISICL
    
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

        # If fixedres was specified, then make list of chain/resi of fixed residues
        if (fixedres != "") {
            n = split(fixedres, a, " ")
            for (i=1; i<=n; i++) {
                chain = substr(a[i], 1, 1)
                resi = substr(a[i], 2)
                for (j=0-fixedrange; j<=0+fixedrange; j++) {
                    fixcheck[chain][resi + j] = 1
                }
            }
        }
    }
    
    $1 ~ /^(ATOM|HETATM)$/ {
        tlc[$5][$6] = $4         # residue name (three-letter code)
        bfa[$5][$6][$3] = $11    # b-factor / pLDDT
        crd[$5][$6][$3][1] = $7  # x coordinate
        crd[$5][$6][$3][2] = $8  # y coordinate
        crd[$5][$6][$3][3] = $9  # z coordinate
    }
    
    $1 !~ /(ss_str|perc_loop|perc_sheet|perc_helix|totres)$/ {
        lin[++l] = $0
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
    
        # Output original file if inplace was specified
        if (inplace::enable) {
            for (i=1; i<=l; i++) { print lin[i] }
        }

        # Get ss_str
        ss_str = ""
        for (chain in crd) for (resi in crd[chain]) {
            ss_str = ss_str ss[chain][resi]
        }
    
        # Calculate percent helix/sheet/loop
        totres = split(ss_str, a, "")
        for ( i=1; i<=totres; i++ ) {
            ss_count[a[i]]++
        }
        for ( N in ssp ) {
            ss_perc[N] = 100 * (ss_count[N]+0) / totres
        }
    
        # Output secondary structure composition
        print "perc_helix " ss_perc["H"]
        print "perc_sheet " ss_perc["E"]
        print "perc_loop " ss_perc["L"]
        print "ss_str " ss_str
        print "totres " totres
        delete ss_count ; delete ss_perc

        # If --fixedres was specified, then calculate / output ss composition for only fixed residues
        if (fixedres != "") {
            ss_str = ""
            for (chain in crd) for (resi in crd[chain]) {
                if (chain in fixcheck && resi in fixcheck[chain]) {
                    ss_str = ss_str ss[chain][resi]
                }
            }
            totres = split(ss_str, a, "")
            for ( i=1; i<=totres; i++ ) {
                ss_count[a[i]]++
            }
            for ( N in ssp ) {
                ss_perc[N] = 100 * (ss_count[N]+0) / totres
            }
            print "fixedres_perc_helix " ss_perc["H"]
            print "fixedres_perc_sheet " ss_perc["E"]
            print "fixedres_perc_loop " ss_perc["L"]
            print "fixedres_ss_str " ss_str
            print "fixedres_totres " totres
        }
    
        delete tlc ; delete bfa ; delete crd ; delete reg ; delete ss
        delete lin ; l = 0
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
# <<< dss_disicl () <<< }}}

inpdbs=(${args[@]:1})

dss_disicl ${inpdbs[@]}
