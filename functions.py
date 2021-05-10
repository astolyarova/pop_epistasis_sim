import sys
sys.path.append("/mnt/.local")
sys.path.append("/mnt/.local/lib/python2.7/site-packages")
import fisher
from scipy.stats.stats import pearsonr
from math import sqrt

def maf(s):
    l = list(set(list(s)))
    l = [s.count(x) for x in l]
    return min(l)

def read_biall(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ] = line
            [A, a] = list(set(col))
            col = list(col)
            col = [0 if x==A else 1 for x in col]
            out[pop][int(position)] = col
    return scaf, out


def read_biall_wo(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col] = line
            [A, a] = list(set(col))
            col = list(col)
            col = [0 if x==A else 1 for x in col]
            out[pop][int(position)] = col
    return scaf, out


def read_biall_ann(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e] = line
            pop = 'USA'
            [A, a] = list(set(col))
            col = list(col)
            if col.count(A) < col.count(a):
                [A, a] = [a, A]
            col = [0 if x==A else 1 for x in col]
            out[pop][int(position)] = (col, typ, rna, exnum, b, e)
    return scaf, out


def read_biall_annaa(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum,b1,e1, aa1,aa2] = line
            [A, a] = list(set(col))
            col = list(col)
            col = [0 if x==A else 1 for x in col]
            out[pop][int(position)] = (col, typ, rna, exnum, aa1,aa2)
    return scaf, out


def read_biall_ann_maf1(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e] = line
            [A, a] = list(set(col))
            col = list(col)
            col = [0 if x==A else 1 for x in col]
            if maf(col) > 1:
                out[pop][int(position)] = (col, typ, rna, exnum, b, e)
    return scaf, out


def read_biall_ann_c(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e] = line
            [A, a] = list(set(col))
            col = list(col)
            col = [0 if x==A else 1 for x in col]
            out[pop][int(position)] = (col, typ, rna, exnum, b, e, A, a)
    return scaf, out


def read_biall_ann_c1(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e] = line
            [A, a] = list(set(col))
            if col.count(A) < col.count(a):
                [A, a] = [a, A]
            col = list(col)
            col = [0 if x==A else 1 for x in col]
            out[pop][int(position)] = (col, typ, rna, exnum, b, e, A, a)
    return scaf, out


def read_biall_ann_der(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e, d] = line
            out[pop][int(position)] = (col, typ, rna, exnum, b, e, int(d))
    return scaf, out


def read_biall_ann_der_typ(fname, t):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e, d] = line
            if typ in t:
                out[pop][int(position)] = (col, typ, rna, exnum, b, e, int(d))
    return scaf, out


def read_biall_ann_dersc(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            line = line.split()
            [scaf, position, pop, col, typ, rna, exnum, b, e, d] = line
            out[pop][int(position)] = (scaf, col, typ, rna, exnum, b, e, int(d))
    return scaf, out


def count_r2(s1, s2):
    if True:
        c = [x + 10*y for x, y in zip(s1, s2)]
        x = []
        for n in [0, 1, 10, 11]:
            x.append(c.count(n))
            x1 = tuple(x)
        [AB, aB, Ab, ab] = [float(y) for y in x]
        r2 = (AB*ab - aB*Ab)**2 / ((AB+Ab) * (aB+ab) * (AB+aB) * (Ab+ab))
        return r2

def count_r(s1, s2):
    #[A, a] = list(set(s1))
    #s1 = [0 if x==A else 1 for x in s1]
    #[A, a] = list(set(s2))
    #s2 = [0 if x==A else 1 for x in s2]
    if True:
        c = [x + 10*y for x, y in zip(s1, s2)]
        x = []
        for n in [0, 1, 10, 11]:
            x.append(c.count(n))
            x1 = tuple(x)
        [AB, aB, Ab, ab] = [float(y) for y in x]
        r2 = (AB*ab - aB*Ab) / sqrt((AB+Ab) * (aB+ab) * (AB+aB) * (Ab+ab))
        return r2


def count_Dprime(s1, s2):
    if True: #  return corrcoef([code[x] for x in s1],[code[x] for x in s2])
        l = len(s1)
        c = [x + 10*y for x, y in zip(s1, s2)]
        A = 0
        a = 1
        B = 0
        b = 1
        AB = int(A+B*10)
        A = s1.count(A) /float(len(s1))
        B = s2.count(B) /float(len(s1))
        AB = float(c.count(AB)) / float(len(s1))
        D = AB - A*B
        if D < 0:
            Dmax = max(-A*B,-(1-A)*(1-B))
        elif D > 0:
            Dmax = min(A*(1-B),B*(1-A))
        else:
            Dmax = 1.0
    return D/Dmax #/abs(Dm)

def count_r2n(s1, s2):
    [A, a] = list(set(s1))
    s1 = [0 if x==A else 1 for x in s1]
    [A, a] = list(set(s2))
    s2 = [0 if x==A else 1 for x in s2]
    if True:
        c = [x + 10*y for x, y in zip(s1, s2)]
        x = []
        for n in [0, 1, 10, 11]:
            x.append(c.count(n))
            x1 = tuple(x)
        [AB, aB, Ab, ab] = [float(y) for y in x]
        #print '#', AB, aB, Ab, ab
        r2 = (AB*ab - aB*Ab)**2 / ((AB+Ab) * (aB+ab) * (AB+aB) * (Ab+ab))
        return r2


def pab(s1, s2):
    if True:
        c = [x + 10*y for x, y in zip(s1, s2)]
        x = []
        for n in [0, 1, 10, 11]:
            x.append(c.count(n))
            x1 = tuple(x)
        [AB, aB, Ab, ab] = [float(y) for y in x]
        return int(min([AB, aB, Ab, ab]))


def count_Fn(s1, s2):
    [A, a] = list(set(s1))
    s1 = [0 if x==A else 1 for x in s1]
    [A, a] = list(set(s2))
    s2 = [0 if x==A else 1 for x in s2]
    if True: #	return corrcoef([code[x] for x in s1],[code[x] for x in s2])
        l = len(s1)
        c = [x + 10*y for x, y in zip(s1, s2)]
        x = []
        for n in [0, 1, 10, 11]:
            x.append(c.count(n))
        x = [float(y) for y in x]
        F = fisher.pvalue(x[0], x[1], x[2], x[3]).two_tail
        return F


def count_F(s1, s2):
    if True: #	return corrcoef([code[x] for x in s1],[code[x] for x in s2])
        l = len(s1)
        c = [x + 10*y for x, y in zip(s1, s2)]
        x = []
        for n in [0, 1, 10, 11]:
            x.append(c.count(n))
        x = [float(y) for y in x]
        F = fisher.pvalue(x[0], x[1], x[2], x[3]).two_tail
        return F


def count_Dminor(s1, s2):
    if True: #	return corrcoef([code[x] for x in s1],[code[x] for x in s2])
        l = len(s1)
        s1 = [int(x) for x in s1]
        s2 = [int(x) for x in s2]
        c = [x + 10*y for x, y in zip(s1, s2)]
#        AB = int(A+B*10)
        A = 1
        B = 1
        AB = 11
        A = s1.count(A) /float(len(s1))
        B = s2.count(B) /float(len(s1))
        AB = float(c.count(AB)) / float(len(s1))
        D = AB - A*B
        return D#/abs(Dm)

def count_Ddaf(s1, s2,d1,d2):
    if True: #	return corrcoef([code[x] for x in s1],[code[x] for x in s2])
        l = len(s1)
        c = [x + 10*y for x, y in zip(s1, s2)]
        if s1.count(0) < s1.count(1):
            A = 0
            a = 1
        else:
            A = 1
            a = 0
        if s2.count(0) < s2.count(1):
            B = 0
            b = 1
        else:
            B = 1
            b = 0
        if d1 == 'c':
            if A == 0: A = 1
            else: A = 0
        if d2 == 'c':
            if B == 0: B = 1
            else: B = 0
        AB = int(A+B*10)
        A = s1.count(A) /float(len(s1))
        B = s2.count(B) /float(len(s1))
        AB = float(c.count(AB)) / float(len(s1))
        D = AB - A*B
        return D#/abs(Dm)
def count_D(s1, s2):
    if True: #	return corrcoef([code[x] for x in s1],[code[x] for x in s2])
        l = len(s1)
        c = [x + 10*y for x, y in zip(s1, s2)]
        if s1.count(0) < s1.count(1):
            A = 0
            a = 1
        else:
            A = 1
            a = 0
        if s2.count(0) < s2.count(1):
            B = 0
            b = 1
        else:
            B = 1
            b = 0
        AB = int(A+B*10)
        A = s1.count(A) /float(len(s1))
        B = s2.count(B) /float(len(s1))
        AB = float(c.count(AB)) / float(len(s1))
        D = AB - A*B
        return D#/abs(Dm)

def read_blocks(fname):
    out = {'RUS': {}, 'USA': {}}
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if line[0] != '#' and line.split()[1] == '.':
                line = line.split()
                [block, begin, end] = [line[0], int(line[5]), int(line[6])]
                out['USA'][block] = [begin, end]
    return out


def merge_sites_blocks(sites, blocks):
    out = {'RUS': {}, 'USA': {}}
    for pop in out:
        s = sites[pop]
        b = blocks[pop]
        for x in s:
            bl = False
            for y in b:
                if x >= b[y][0] and x <= b[y][1]:
                    bl = True
                    out[pop][y]  = out[pop].get(y, []) + [s[x]]
            if not bl:
                out[pop][x] = [s[x]]
    return out


def read_slim(fname, num):
    sites = {}
    mutations = {}
    with open(fname, 'r') as f:
        for line in f:
            line = line.split()
            if len(line) > 2:
                if len(line) == 9 and line[2][0] == 'm':
                    [mut_num, mut_id, mut_type, position, sel, dom, pop, generation, maf] = line
                    mutations[int(mut_num)] = (int(position), mut_type)
                if line[0][:2] == 'p1':
                    s = line[2:]
                    sample = int(line[0][3:])
                    for site in s:
                        [site, mut_type] = mutations[int(site)]
                        sites[mut_type] = sites.get(mut_type, {})
                        if site not in sites[mut_type]:
                            sites[mut_type][site] = [0] * num
                        sites[mut_type][site][sample] = 1
    for m in sites:
        for site in sites[m].keys():
            if sites[m][site].count(1) == 0 or sites[m][site].count(0) == 0:
                del sites[m][site]
    return(sites)
