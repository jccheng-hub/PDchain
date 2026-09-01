function def_aa_map() {
    delete aa_map
    aa_map["ALA"]="A" ; aa_map["ARG"]="R" ; aa_map["ASN"]="N" ; aa_map["ASP"]="D"
    aa_map["CYS"]="C" ; aa_map["GLU"]="E" ; aa_map["GLN"]="Q" ; aa_map["GLY"]="G"
    aa_map["HIS"]="H" ; aa_map["ILE"]="I" ; aa_map["LEU"]="L" ; aa_map["LYS"]="K"
    aa_map["MET"]="M" ; aa_map["PHE"]="F" ; aa_map["PRO"]="P" ; aa_map["SER"]="S"
    aa_map["THR"]="T" ; aa_map["TRP"]="W" ; aa_map["TYR"]="Y" ; aa_map["VAL"]="V"
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

function clash_count(obj1, obj2,    nbr, clashes) {
    # Count number of clashes between two objects
    find_nbr(obj1, obj2, 3, nbr) ; clashes = 0
    for (ch in nbr) for (ri in nbr[ch]) for (at in nbr[ch][ri]) clashes++
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
    return sqrt(sum_sq_dev / n)
}

function rdev(obj1, obj2,    d) {
    # Return the average standard deviation of atomic distances
    # Basically, loop pdev for every atom found in obj1, then average
    # obj1: object 1 | obj1[chain][resi][atom][1/2/3] = x/y/z
    # obj2: object 2 | obj2[chain][resi][atom][1/2/3] = x/y/z
    for (ch in obj1) for (ri in obj1[ch]) for (at in obj1[ch][ri]) {
        sum += pdev(obj1[ch][ri][at], obj2,    d) ; n++
    }
    return sum / n
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
    return sqrt(sum / n)
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
    return sum / n
}

function pdst(p, obj,    ch, ri, at, sum, n) {
    # Return average distance between point p and all atoms in obj
    # p: point | p[1/2/3] = x/y/z
    # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z
    for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
        sum += calc_dist(p, obj[ch][ri][at]) ; n++
    }
    return sum / n
}

function rdst(obj1, obj2,    ch, ri, at, sum, n) {
    # Return average distance between all atoms in obj1 and all atoms in obj2
    # obj1: object 1 | obj[chain][resi][atom][1/2/3] = x/y/z
    # obj2: object 2 | obj[chain][resi][atom][1/2/3] = x/y/z
    for (ch in obj1) for (ri in obj1[ch]) for (at in obj1[ch][ri]) {
        sum += pdst(obj1[ch][ri][at], obj2) ; n++
    }
    return sum / n
}

function midpoint(obj, mp,    ch, ri, at, x, y, z, n) {
    # Get midpoint of all atoms in obj
    # obj: object | obj[chain][resi][atom][1/2/3] = x/y/z coord
    # mp: midpoint | mp[1/2/3] = x/y/z
    delete mp
    for (ch in obj) for (ri in obj[ch]) for (at in obj[ch][ri]) {
        x+=obj[ch][ri][at][1] ; y+=obj[ch][ri][at][2] ; z+=obj[ch][ri][at][3] ; n++
    }
    mp[1] = x / n
    mp[2] = y / n
    mp[3] = z / n
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
