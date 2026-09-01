#!/usr/bin/env bash
usage () { cat << EOF
Usage: calc_scafscore INPDB... [OPTIONS]
Calculate scaffold-only metrics for input PDBs. The backbone atoms are used to
infer a CB position at each residue, resulting in a poly-ala scaffold. This
poly-ala scaffold is then used for the calculation various metrics meant to
assess the quality of the scaffold in a sequence-agnostic way.

If focus residues are provided (with --chri and/or --resn), then additional
metrics for those residues will also be provided. If designing a ligand binder
or a protein binder, then the focus residues would be the ligand or hotspot
residues, respectively. Additional explanation for the metrics are provided in
the Metrics section.

Parameters:
    INPDB                   Input PDB

Options:
  --chri [str]              Focus residues (e.g. ligands) using chain+resi
                            format.
                            E.g. --chri X1
  --resn [str]              Focus residues using residue name.
                            E.g. --resn IGP
  --nbr_cut [float]         Distance cutoff for counting neighbors. Used for
                            calculating certain metrics.
  --help                    Display this help and exit

Metrics:
    mp_dst                  Average distance between the midpoint of non-focus
                            atoms and all non-focus atoms. A proxy metric for
                            protein size.
    mp_dev                  Standard deviation of atomic distances calculated in
                            mp_dst. A proxy metric for a hollow cavity.
    fc_dst                  Average distance between focus atoms and all
                            non-focus atoms within [nbr_cut] angstroms.
    fc_dev                  Standard deviation of atomic distances calculated in
                            fc_dst.
    clash                   Number of atomic clashes between focus atoms and
                            the poly-Ala backbone. Clash defined as <3 angstrom
                            distance. Does not count same-chain clashes.
    rog_ala                 Approximate radius of gyration. Treats all poly-Ala
                            heavy atoms as equal mass.
    fc_compact              Applies ROG equation to each focus atom, then
                            average those values. Only non-focus atoms within
                            [nbr_cut] angstroms are included. Lower values mean
                            that the non-focus atoms are more compact around
                            focus atoms.
    mp_fc_dst               Average distance between the midpoint of non-focus
                            atoms and all focus atoms. A proxy metric for
                            the positioning of focus residues.
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

val_opts=(chri resn nbr_cut)
bool_opts=(help)

chri=""
resn=""
nbr_cut="100"

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

# Initialize
inpdbs=(${args[@]:1})
mets="mp_dst mp_dev rog_ala fc_dst fc_dev clash fc_compact mp_fc_dst"

# >>> gawk_script >>> {{{
gawk_script='
    BEGIN {
        if (!chri)    chri = "'"$chri"'"
        if (!resn)    resn = "'"$resn"'"
        if (!nbr_cut) nbr_cut = '"$nbr_cut"'
    
        if (chri != "" || resn != "") {
            chri = expand_range(chri)
            split(chri, a) ; for (i in a) foc_chri[a[i]] = 1
            split(resn, a) ; for (i in a) foc_resn[a[i]] = 1
        }
    }
    
    $1~/^('"${mets// /|}"')$/ { next } 1
    
    $1!~/^(ATOM|HETATM)$/ || $NF == "H" { next }
    
    $5 $6 in foc_chri || $4 in foc_resn {
        # Format: obj[chain][resi][atom][1/2/3] = x/y/z coord
        foc[$5][$6][$3][1] = $7    # x coord
        foc[$5][$6][$3][2] = $8    # y coord
        foc[$5][$6][$3][3] = $9    # z coord
        next
    }
    
    $3~/^(N|CA|C|O)$/ {
        nfc[$5][$6][$3][1] = $7
        nfc[$5][$6][$3][2] = $8
        nfc[$5][$6][$3][3] = $9
    }
    
    ENDFILE {
        # Add CB coordinates
        for (ch in nfc) for (ri in nfc[ch]) {
            if ("N" in nfc[ch][ri] && "CA" in nfc[ch][ri] && "C" in nfc[ch][ri]) {
                add_cb(nfc[ch][ri])
            }
        }
    
        midpoint(nfc, mp)         # Midpoint of non-focus residues
        mp_dst = pdst(mp, nfc)    # Mean midpoint distance
        mp_dev = pdev(mp, nfc)    # Mean midpoint deviation
        rog_ala = calc_arog(nfc)  # Output ROG for poly-alanine
        print "mp_dst " mp_dst
        print "mp_dev " mp_dev
        print "rog_ala " rog_ala
    
        # Exit if no focus residues were provided
        if (!(chri != "" || resn != "")) exit
        
        fc_dst = rdst(foc, nfc)           # Mean distance from focus residues
        fc_dev = rdev(foc, nfc, nbr_cut)  # Deviation for focus residues
        clash = clash_count(foc, nfc)     # Number of clashes
        fc_compact = calc_local_arog(foc, nfc, nbr_cut)  # Compactness around focus atoms
        mp_fc_dst = pdst(mp, foc)         # Mean distance between midpoint and focus atoms
        print "fc_dst " fc_dst
        print "fc_dev " fc_dev
        print "clash " clash
        print "fc_compact " fc_compact
        print "mp_fc_dst " mp_fc_dst
    
        delete foc ; delete nfc
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
    
    function clash_count(obj1, obj2,
                         objA, objB, nbr, ch, ri, at,
                         ch1, ri1, at1, ch2, ri2, at2, i, clashes)
    {
        # Count number of clashes between two objects
        # Excludes same-chain clashes
        clashes = 0
        for (ch1 in obj1) {
            for (ch2 in obj2) {
                if (ch1 == ch2) { continue }
                delete objA ; delete objB
                for (ri1 in obj1[ch1]) for (at1 in obj1[ch1][ri1]) for (i=1; i<=3; i++) {
                    objA[ch1][ri1][at1][i] = obj1[ch1][ri1][at1][i]
                }
                for (ri2 in obj2[ch2]) for (at2 in obj2[ch2][ri2]) for (i=1; i<=3; i++) {
                    objB[ch2][ri2][at2][i] = obj2[ch2][ri2][at2][i]
                }
                find_nbr(objA, objB, 3, nbr)
                for (ch in nbr) for (ri in nbr[ch]) for (at in nbr[ch][ri]) clashes++
            }
        }
        return clashes
    }
    
    function pdev(p, obj,    d, ch, ri, at, d_tmp, dist, n, sum, mean, i, dev, sum_sq_dev) {
        # Return standard deviation of atomic distances from specified point
        # p: point | p[1/2/3] = x/y/z
        # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z
        # d: distance | float
        for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
            d_tmp = calc_dist(p, obj[ch][ri][at])
            if (d && d_tmp > d) continue
            dist[++n] = d_tmp ; sum += d_tmp
        }
        mean = sum / n
        for (i in dist) { dev = dist[i] - mean ; sum_sq_dev += dev * dev }
        if (n > 0) { return sqrt(sum_sq_dev / n) } else { return 9999 }
    }
    
    function rdev(obj1, obj2,    d, ch, ri, at, sum, n) {
        # Return the average standard deviation of atomic distances
        # Basically, loop pdev for every atom found in obj1, then average
        # Excludes same-chain comparisons
        # obj1: object 1 | obj1[chain][resi][atom][1/2/3] = x/y/z
        # obj2: object 2 | obj2[chain][resi][atom][1/2/3] = x/y/z
        for (ch in obj1) {
            if (ch in obj2) continue
            for (ri in obj1[ch]) for (at in obj1[ch][ri]) {
                sum += pdev(obj1[ch][ri][at], obj2,    d) ; n++
            }
        }
        if (n > 0) { return sum / n } else { return 9999 }
    }
    
    function check_object(obj,    ch, ri, at, x, y, z) {
        # Print out contents of object
        # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z
        for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
            x = obj[ch][ri][at][1]
            y = obj[ch][ri][at][2]
            z = obj[ch][ri][at][3]
            printf "%-8s%-4s%-16s%-16s%-16s\n", ch ri, at, x, y, z
        }
    }
    
    function calc_arog(obj,    mp, i, j, d, sum, n) {
        # Calculate approximate ROG: treats all atoms as identical mass
        # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z
        # mp: midpoint | mp[1/2/3] = x/y/z
        if (typeof(mp) != "array") midpoint(obj, mp)
        for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
            d = calc_dist(mp, obj[ch][ri][at]) ; sum+=d*d ; n++
        }
        if (n > 0) { return sqrt(sum / n) } else { return 9999 }
    }
    
    function calc_local_arog(obj1, obj2, dcut,    ch, ri, at, i, dot, nbr, sum, n) {
        # Calculate arog values for obj2 using every obj1 atom as midpoint, then average
        # obj1: object 1 | obj[chain][resi][atom][1/2/3] = x/y/z
        # obj2: object 2 | obj[chain][resi][atom][1/2/3] = x/y/z
        # dcut: distance cutoff | float
        for (ch in obj1) for (ri in obj1[ch]) for (at in obj1[ch][ri]) {
            for (i=1; i<=3; i++) dot[ch][ri][at][i] = obj1[ch][ri][at][i]
            find_nbr(dot, obj2, dcut, nbr)
            sum += calc_arog(nbr, dot[ch][ri][at]) ; n++
            delete dot; delete nbr
        }
        if (n > 0) { return sum / n } else { return 9999 }
    }
    
    function pdst(p, obj,    ch, ri, at, sum, n) {
        # Return average distance between point p and all atoms in obj
        # p: point | p[1/2/3] = x/y/z
        # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z
        for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
            sum += calc_dist(p, obj[ch][ri][at]) ; n++
        }
        if (n > 0) { return sum / n } else { return 9999 }
    }
    
    function rdst(obj1, obj2,    ch, ri, at, sum, n) {
        # Return average distance between all atoms in obj1 and all atoms in obj2
        # obj1: object 1 | obj[chain][resi][atom][1/2/3] = x/y/z
        # obj2: object 2 | obj[chain][resi][atom][1/2/3] = x/y/z
        for (ch in obj1) for (ri in obj1[ch]) for (at in obj1[ch][ri]) {
            sum += pdst(obj1[ch][ri][at], obj2) ; n++
        }
        if (n > 0) { return sum / n } else { return 9999 }
    }
    
    function midpoint(obj, mp,    ch, ri, at, x, y, z, n) {
        # Get midpoint of all atoms in obj
        # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z coord
        # mp: midpoint | mp[1/2/3] = x/y/z
        delete mp
        for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
            x+=obj[ch][ri][at][1] ; y+=obj[ch][ri][at][2] ; z+=obj[ch][ri][at][3] ; n++
        }
        if ( n > 0 ) { mp[1] = x / n ; mp[2] = y / n ; mp[3] = z / n }
    }
    
    function add_cb(res,    n, ca, c, v, v1, v2, v3, vf, norm_n, norm_c, norm_k, norm_b, vec_cb) {
        # Adds CB atom based on positions of N, CA, and C
        # res: res[atom][1/2/3] = x/y/z
        # res will be modified to have res["CB"][1/2/3] = x/y/z
    
        for (i in res["N"])  n[i]  = res["N"][i]
        for (i in res["CA"]) ca[i] = res["CA"][i]
        for (i in res["C"])  c[i]  = res["C"][i]
    
        vec_sub(n, ca, v)            ; vec_norm(v, norm_n)
        vec_sub(c, ca, v)            ; vec_norm(v, norm_c)
        vec_cross(norm_n, norm_c, v) ; vec_norm(v, norm_k)
        vec_add(norm_n, norm_c, v)   ; vec_norm(v, norm_b)
    
        vec_scale(norm_b, -0.57735, v1)
        vec_scale(norm_k,  0.81649, v2)
        vec_add(v1, v2, v3)
        vec_scale(v3, 1.52, vec_cb)
    
        vec_add(ca, vec_cb, vf)
    
        for (i=1; i<=3; i++) res["CB"][i] = vf[i]
    }
    
    function calc_dist(a, b,    dx, dy, dz) {
        # a: a[1/2/3] = x/y/z
        # b: b[1/2/3] = x/y/z
        dx = a[1]-b[1] ; dy = a[2]-b[2] ; dz = a[3]-b[3]
        return sqrt( dx*dx + dy*dy + dz*dz )
    }
    
    function calc_dist_sq(a, b,    dx, dy, dz) {
        # a: a[1/2/3] = x/y/z
        # b: b[1/2/3] = x/y/z
        dx = a[1]-b[1] ; dy = a[2]-b[2] ; dz = a[3]-b[3]
        return dx*dx + dy*dy + dz*dz
    }
    
    function vec_add(v1, v2, vf) {
        delete vf ; vf[1] = v1[1]+v2[1] ; vf[2] = v1[2]+v2[2] ; vf[3] = v1[3]+v2[3]
    }
    
    function vec_sub(v1, v2, vf) {
        delete vf ; vf[1] = v1[1]-v2[1] ; vf[2] = v1[2]-v2[2] ; vf[3] = v1[3]-v2[3]
    }
    
    function vec_cross(v1, v2, vf) {
        delete vf
        vf[1] = v1[2]*v2[3] - v1[3]*v2[2] 
        vf[2] = v1[3]*v2[1] - v1[1]*v2[3]
        vf[3] = v1[1]*v2[2] - v1[2]*v2[1]
    }
    
    function vec_scale(vi, scalar, vf) {
        delete vf
        vf[1] = scalar * vi[1]
        vf[2] = scalar * vi[2]
        vf[3] = scalar * vi[3]
    }
    
    function vec_mag(v,    a_sq, b_sq, c_sq) {
        a_sq = v[1] * v[1]
        b_sq = v[2] * v[2]
        c_sq = v[3] * v[3]
        return sqrt( a_sq + b_sq + c_sq )
    }
    
    function vec_norm(vi, vf,    mag) {
        delete vf ; mag = vec_mag(vi) ; vec_scale(vi, 1/mag, vf)
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
'
# <<< gawk_script <<< }}}

# Calculate scaffold scores
echo "Calculating scaffold scores for ${#inpdbs[@]} pdbs..."
gawk -i inplace -f <(echo "$gawk_script") ${inpdbs[@]}

# Report
awk '
$1~/^('"${mets// /|}"')$/ {
    f = FILENAME ; sub(/^.*\//, "", f) ; sub(/\..*$/, "", f) ; print f, $0
}' ${inpdbs[@]}
