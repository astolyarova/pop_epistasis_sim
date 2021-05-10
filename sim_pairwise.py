import sys
import numpy as np
import FFPopSim as h
import functions as fn
import argparse
import scipy
import scipy.cluster.hierarchy as sch



def init_pop(add_sel, epi_values):
    pop = h.haploid_highd(L)
    pop.carrying_capacity = N
    pop.recombination_model = h.CROSSOVERS
    pop.mutation_rate = mu       # mutation rate per site per generation
    pop.outcrossing_rate = 1    # probability of sexual reproduction per gen
    pop.crossover_rate = rec      # probability of crossover per site per gen

    pop.set_fitness_additive(add_sel)
    for pair in epi_values:
        pop.add_fitness_coefficient(epi_values[pair], [pair[0],pair[1]])
    #g = np.random.choice([False,True],L,[0.5,0.5])
    pop.set_wildtype(N)

    pop.evolve(0) # it's crucial don't know why
    return pop


def gen_random_value(mean,scale):
    if scale == 0:
        return mean
    else:
        return np.random.gamma(mean/scale,scale,1)[0]


def pop_snps(pop,size):

    gclones = pop.get_genotypes()
    nclones = pop.get_clone_sizes()
    sclones = np.random.choice(len(gclones), size, p=nclones/float(sum(nclones)), replace=True)
    gclones = [gclones[i] for i in sclones]
    nclones = [1] * len(sclones)

    # extract SNPs
    snps = []
    gclones = np.array(gclones)
    for i in range(L):
        #if args.mode not in ['sweep'] or i not in [sweep_site,sweep_site2]:
        x = gclones[:, i]
        maf = fn.maf(list(x))
        if maf > size * freq and maf < size * (1 - freq):
            snps.append(i)
    snps_cols = {}
    for i in snps:
        col = []
        for clone in range(len(nclones)):
            col += [int(gclones[clone][i])] * nclones[clone]
        snps_cols[i] = col
    return gclones, nclones, sclones, snps, snps_cols


def pop_stats(pop,sample):
    [gclones, nclones, sclones, snps, snps_cols] = pop_snps(pop, sample)

    # count SNP stats (number, LD)

    LD = {'s': [],'n': [],'all': []}
    LDr2 = {'s': [],'n': [],'all': []}
    mafs = {'s':[],'n':[]}
    n_number = 0
    s_number = 0
    pi = {'s':[],'n':[],'all':[]}
    for i in range(L):
        for n1 in range(sample):
            for n2 in range(n1 + 1, sample):
                if gclones[n1][i] == gclones[n2][i]:
                    pi['all'].append(0)
                    if i in epis:
                        pi['n'].append(0)
                    else:
                        pi['s'].append(0)
                else:
                    pi['all'].append(1)
                    if i in epis:
                        pi['n'].append(1)
                    else:
                        pi['s'].append(1)
    for i in snps:
        if i in epis or i in nonsyn:
            n_number += 1
            mafs['n'].append(fn.maf(snps_cols[i]))
        else:
            s_number += 1
            mafs['s'].append(fn.maf(snps_cols[i]))
        for j in snps:
            if j > i:
                if (i not in epis and i not in nonsyn) and (j not in epis and j not in nonsyn):
                    t = 's'
                elif (i in epis or i in nonsyn) and (j in epis or j in nonsyn):
                    t = 'n'
                else:
                    t = 'sn'
                if t in ['s', 'n']:
                    LDr2[t].append(fn.count_r2(snps_cols[i],snps_cols[j]))
                    LD[t].append(fn.count_Dprime(snps_cols[i],snps_cols[j]))
                LDr2['all'].append(fn.count_r2(snps_cols[i],snps_cols[j]))
                LD['all'].append(fn.count_Dprime(snps_cols[i],snps_cols[j]))

    # count genotypes stats (fitness, epistatic fitness)

    fitness = []
    snp_epistatic_fitness = []  # makes sense only for pairwise epistasis model
    pfit = pop.get_fitness_statistics()
    for n in range(len(gclones)):
        fitness.append(pop.get_fitness(sclones[n]))
        a=sum((2 * gclones[n] - 1) * add_sel)
        E = 0
        for i in snps:
            for j in snps:
                if j > i:
                    E += (epi_values.get((i,j),0)+epi_values.get((j,i),0)) * ((gclones[n][i] == gclones[n][j]) * 2 - 1)
        snp_epistatic_fitness.append(E)

    return pfit.mean, pfit.variance, np.mean(snp_epistatic_fitness), s_number, n_number, LD, LDr2, np.mean(mafs['s']), np.mean(mafs['n']), np.mean(pi['all']),np.mean(pi['s']),np.mean(pi['n'])


def pop_distance(pop):
    [gclones, nclones, sclones, snps, snps_cols] = pop_snps(pop, sample)

    LD = {'s': {},'n': {},'sn': {}}
    for i in snps:
        for j in snps:
            if j > i:
                dist = ((j - i) / 50) * 50
                if i not in epis and j not in epis:
                    t = 's'
                elif i in epis and j in epis:
                    t = 'n'
                else:
                    t = 'sn'
                LD[t][dist] = LD[t].get(dist, [])
                LD[t][dist].append(fn.count_r2(snps_cols[i],snps_cols[j]))
    return LD


def pop_aln(pop):
    [gclones, nclones, sclones, snps, snps_cols] = pop_snps(pop, sample)
    Y = sch.linkage(gclones, method='weighted')
    Z = sch.dendrogram(Y, no_plot=True)['ivl']
    snps_cols1 = {}
    if sum(gclones[int(Z[0])]) < sum(gclones[int(Z[len(Z) - 1])]):
        Z = Z[::-1]
    for i in range(L):
        if i not in snps_cols:
            snps_cols1[i]= [int(gclones[0,i])] * len(gclones)
        else:
            snps_cols1[i] = [snps_cols[i][int(x)] for x in Z]
    return snps_cols1


def pop_fitness(pop):
    [gclones, nclones, sclones, snps, snps_cols] = pop_snps(pop, sample)
    efitness = []
    efitness_random = []
    snps1 = []
    comp_all = []
    #print epi_values1
    for i in snps:
        if i in epis:
            snps1.append(i)
    snps = snps1
   # random_snps = list(np.random.choice(epis,len(snps)))
    random_snps = {}
    for i in snps:
        random_snps[i] = np.random.permutation(snps_cols[i])

    ngen = []
    for n in gclones:
        ngen0 = 0
        for i in snps:
            ngen0 += int(n[i])
        ngen.append(ngen0)

    E = 0
    E_random = 0
    comp = {}
    pres = {}
    for i in snps:
        comp[i] = [0] * sample
        pres[i] = [0] * sample
        for n in range(sample):
            if snps_cols[i][n]: pres[i] = 1
        for j in snps:
            if j < i:
                for n in range(sample):
                    e = epi_values.get((j,i),0) * ((snps_cols[i][n] == snps_cols[j][n]) * 2 - 1)
                    E += e
                    E_random += epi_values.get((j,i),0) * ((random_snps[i][n] == random_snps[j][n]) * 2 - 1)
                    if e and snps_cols[j][n] and snps_cols[i][n]:
                        comp[i][n] = 1
                        comp[j][n] = 1
    for i in comp:
        comp_all.append(np.mean(comp[i])/np.mean(pres[i]))
    efitness.append(E)
    efitness_random.append(E_random)
    return efitness, efitness_random, np.mean(ngen),np.mean(comp_all) #len(snps)


def pop_arg(pop):
    ev = []
    for kk in range(100):
        rec_events = 0
        gclones = pop.get_genotypes()
        nclones = pop.get_clone_sizes()
        pop0 = h.haploid_highd(L)
        pop0.track_locus_genealogy(range(L))
        pop0.set_genotypes(gclones,nclones)
        pop0.set_fitness_additive(ff)
        for pair in epi_values:
            pop0.add_fitness_coefficient(epi_values[pair], [pair[0],pair[1]])
        pop0.recombination_model = h.CROSSOVERS
        pop0.mutation_rate = 0       # mutation rate per site per generation
        pop0.outcrossing_rate = 1    # probability of sexual reproduction per gen
        pop0.crossover_rate = rec      # probability of crossover per site per gen
        pop0.evolve(0)
        pop0.evolve(1)
        for i in range(L - 1):
            j = i + 1
            tree0 = pop0.genealogy.get_tree(i).print_newick()
            tree1 = pop0.genealogy.get_tree(j).print_newick()
            rec_events += int(tree0!=tree1)
        for kkk in range(1): ev.append(rec_events)
    return ev


def pop_out(pop):
    if args.fmt == 'stats':
        [fita, fite, fitsnps, s, n, LD, LDr2, mafs, mafn, pop_pi,ps,pn] = pop_stats(pop,sample)
        print xxx, sigma*2, L, N, mu,rec, epi_mode, is_epistasis, args.p, ii*gen_step, fita, fite, fitsnps, s, n, np.mean(LD['s']), np.mean(LD['n']),np.mean(LD['all']),np.mean(LDr2['s']),np.mean(LDr2['n']),np.mean(LDr2['all']), mafs, mafn, pop_pi,ps,pn
    elif args.fmt == 'alignment':
        snps_cols = pop_aln(pop)
        snps = snps_cols.keys()
        snps.sort()
        for i_rank, i in enumerate(snps):
            for g in range(len(snps_cols[i])):
                print xxx, sigma*2, L, N, mu,rec, epi_mode, is_epistasis, ii*gen_step, i, i_rank, g, snps_cols[i][g], int(i in epis)
    elif args.fmt == 'fitness':
        [E, E_random, nsnps, comp] = pop_fitness(pop)
        for g in range(len(E)):
            print xxx, sigma*2, L, N, mu,rec, epi_mode, is_epistasis, args.p,ii*gen_step, E[g], E_random[g], nsnps, comp
    elif args.fmt == 'distance':
        LD = pop_distance(pop)
        for allele_type in LD:
            for distance in LD[allele_type]:
                print xxx, sigma*2, L, N, mu,rec, epi_mode, is_epistasis, ii*gen_step, allele_type, distance, np.mean(LD[allele_type][distance]), len(LD[allele_type][distance])
    elif args.fmt == 'arg':
        rec_events = pop_arg(pop)
        for xx in range(1):
            for x in rec_events:
                print xxx, sigma*2, L, N, mu,rec, epi_mode, is_epistasis, ii*gen_step, x

parser = argparse.ArgumentParser()

parser.add_argument('-i', type=int, action='store', dest='i', help="number of simulations")
parser.add_argument('-N', action='store', dest='N', help="population size")
parser.add_argument('-L', type=int, action='store', dest='L', help="sequence length")
parser.add_argument('--mu', action='store', dest='mu', help="mutation rate(s)")
parser.add_argument('--epi_mode', choices=['isolated'],action='store', dest='epi_mode', help="epstasis mode; only 'isolated' is implemented so far")
parser.add_argument('--rec', action='store', dest='rec', help="recombination rate(s) population scaled")
parser.add_argument('--gen', action='store', dest='gen', help="generations to simulate, number of steps")
parser.add_argument('--sigma', action='store', dest='sigma', help="selection coefficient(s) of deleterious mutations (population-scaled)")
parser.add_argument('--sample', type=int,action='store', dest='sample', help="sample size")
parser.add_argument('--freq', type=float,action='store', dest='freq', help="frequency threshold")
parser.add_argument('--p', type=float,action='store', dest='p', help="scale coefficient of DFE")
parser.add_argument('--mode', action='store', dest='mode', help="population mode",choices=['panmixic','admixture','balancing','sweep'])
parser.add_argument('--outfmt', action='store', dest='fmt', help="output",choices=['stats','alignment','fitness','distance','arg'])
args = parser.parse_args()

L = args.L
epi_mode = args.epi_mode
sigma_list = [float(x)/2.0 for x in args.sigma.split(',')] # normalize for ffpopsim
rec_list = [float(x) for x in args.rec.split(',')]
mu_list = [float(x) for x in args.mu.split(',')]
freq = args.freq
sample = args.sample
N = int(args.N)

[gens, gen_ints] = [int(x) for x in args.gen.split(',')]
gen_step = gens / gen_ints

if args.fmt == 'stats':
    print 'xxx', 'sigma','L','N', 'mu','rec', 'epi_mode', 'epi', 'p','gen', 'fita', 'fite','fitsnps', 's', 'n', 'LDs','LDn', 'LD', 'r2s','r2n','r2', 'mafs','mafn','pi','ps','pn'
elif args.fmt == 'alignment':
    print 'xxx', 'sigma','L','N', 'mu','rec', 'epi_mode', 'epi', 'gen', 'position', 'position_rank', 'genotype', 'allele', 'allele_type'
elif args.fmt == 'fitness':
    print 'xxx', 'sigma','L','N', 'mu','rec', 'epi_mode', 'epi', 'p','gen', 'E', 'E_random', 'n', 'comp'
elif args.fmt == 'distance':
    print 'xxx', 'sigma','L','N', 'mu','rec', 'epi_mode', 'epi', 'gen', 'allele_type', 'distance', 'LD', 'n'
elif args.fmt == 'arg':
    print 'xxx', 'sigma','L','N', 'mu','rec', 'epi_mode', 'epi', 'gen', 'events'

mu = mu_list[0]
for jjj in range(args.i):
    for rec0 in rec_list:
        for sigma in sigma_list:
            for mu in mu_list:
                rec = rec0 /float(N) # normalization
                epi_list = [0.0, 1.0]
                for is_epistasis in epi_list:
                    xxx = np.random.random(1)[0]

                    epis = list(range(1, L, 3)) + list(range(2, L, 3))  # sites under selection
                    np.random.shuffle(epis)
                    k=int(round(len(epis)*(1.0-0.0)/2.0))
                    if k > 0:
                        epis1= epis[:k]
                        epis2 = epis[-k:][::-1]
                    else:
                        epis1 = []; epis2 = []
                    nonsyn=epis[k:-k]
                    epis=epis1+epis2
                    sel_site = epis[0]
                    epi_values = {}
                    if args.mode == 'balancing':
                        epis = epis[1:]
                    epis1 = epis[:len(epis)/2]
                    epis2 = epis[len(epis)/2:]
                    sigma_a = sigma * (1 - is_epistasis) / float(N)
                    sigma_e = sigma *  is_epistasis / float(N)
                    sigma_sum = sigma_a + sigma_e
                    if True: #sigma_e > 0:
                        if True:
                        #if epi_mode == 'abundant':
                        #    for i in range(len(epis)):
                        #        for j in range(len(epis)):
                        #            if j > i:
                        #                [ii,jj] = [min(epis[i],epis[j]),max(epis[i],epis[j])]
                        #                value = gen_random_value(sigma_sum/(len(epis)-1),args.p)
                        #                epi_values[(ii,jj)] = value
                        #elif epi_mode == 'isolated':
                            for i in range(len(epis1)):
                                [ii,jj] = [min(epis1[i],epis2[i]),max(epis1[i],epis2[i])]
                                value = gen_random_value(sigma_sum,args.p/float(N))
                                if args.mode not in ['sweep','balancing'] or i != 0:
                                    epi_values[(ii,jj)] = value

                    if sigma_e:
                        epi_val_use = epi_values
                    else:
                        epi_val_use = {}

                    add_sel = [0.0] * L
                    if sigma_a:
                        for i in epis:
                            value = gen_random_value(sigma_a,args.p/float(N))
                            add_sel[i] = -value

                    if args.mode == 'panmixic':
                        pop = init_pop(add_sel, epi_val_use)
                    elif args.mode == 'balancing':
                        pop = init_pop(add_sel, epi_val_use)
                        #pop.evolve(N*20)
                    elif args.mode == 'admixture':
                        pop1 = init_pop(add_sel, epi_val_use)
                        pop2 = init_pop(add_sel, epi_val_use)
                        pop1.evolve(N*20)
                        pop2.evolve(N*20)
                        pop1 = list(pop1.random_genomes(N/2))
                        pop2 = list(pop2.random_genomes(N/2))
                        pop = init_pop(add_sel, epi_val_use)
                        pop.set_genotypes(pop1 + pop2, [1] * N)
                        pop.evolve(0)

                    if args.mode == 'sweep':
                        sweep_site = sel_site #epis1[0]
                        add_sel[sweep_site] = -1e5
                        pop = init_pop(add_sel, epi_val_use)
                        pop.evolve(N*20)
                        sweep_site2 = sweep_site
                        add_sel[sweep_site] = 100 / float(N)
                        pop.set_fitness_additive(add_sel)

                    ii = 0
                    if args.mode == 'sweep':
                        gen_step = 1
                        sweep_freq = pop.get_allele_frequency(sweep_site)
                        while sweep_freq < .5 and ii < 100000:
                            sweep_freq = pop.get_allele_frequency(sweep_site)
                            sweep_freq1 = sweep_freq
                            if sweep_freq1 == 0:
                                genomes = list(pop.random_genomes(N))
                                genome_sweep = list(genomes[0][:sweep_site]) + [1] + list(genomes[0][sweep_site + 1:])
                                pop.set_genotypes([genome_sweep] + genomes[1:], [1] *N)
                            pop.evolve(1)
                            ii += gen_step

                    else:
                        while ii < gen_ints:
                            if args.fmt not in ['arg','fitness'] or ii > 0:
                                pop_out(pop)
                            if args.mode == 'balancing':
                                for ii_1 in range(gen_step/10):
                                    balanc_freq = pop.get_allele_frequencies()[sel_site]
                                    balanc_sel = 1.0 * (balanc_freq - 0.5) * 100 / float(N)
                                    add_sel1 = add_sel[:sel_site] + [-balanc_sel] + add_sel[sel_site+1:]
                                    pop.set_fitness_additive(add_sel1)
                                    pop.evolve(10)
                            else:
                                pop.evolve(gen_step)
                            ii += 1
                    pop_out(pop)
