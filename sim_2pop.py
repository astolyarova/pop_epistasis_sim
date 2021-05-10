import sys
sys.path.append('/home/avstolyarova/.local/lib/python2.7/site-packages')
#sys.path.append('/opt/rh/python27/root/usr/lib64/python2.7/site-packages')
#sys.path.append('~/tools/pylib/FFPopSim')
import numpy as np
import FFPopSim as h
#from FFPopSim import FFPopSim as h
import functions as fn
import argparse
import scipy
import scipy.cluster.hierarchy as sch


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
        x = gclones[:, i]
        maf = fn.maf(list(x))
        if maf > size * freq and maf < size * (1 - freq):
        #if sum(x) not in [0, len(x),1,len(x)-1]:
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

    LD = {'s': [],'n': [],'sn': []}
    LDr2 = {'s': [],'n': [],'sn': []}
    mafs = {'s':[],'n':[]}
    n_number = 0
    s_number = 0
    for i in snps:
        if i in epis:
            n_number += 1
            mafs['n'].append(fn.maf(snps_cols[i]))
        else:
            s_number += 1
            mafs['s'].append(fn.maf(snps_cols[i]))
        for j in snps:
            if j > i:
                if i not in epis and j not in epis:
                    t = 's'
                elif i in epis and j in epis:
                    t = 'n'
                else:
                    t = 'sn'
                LDr2[t].append(fn.count_r2(snps_cols[i],snps_cols[j]))
                LD[t].append(fn.count_Dprime(snps_cols[i],snps_cols[j]))

    # count genotypes stats (fitness, epistatic fitness)

    fitness = []
    #additive_fitness = []
    #epistatic_fitness = []
    snp_epistatic_fitness = []  # makes sense only for pairweise epistasis model
    #for i in epis:
    #    for j in epis:
    #        if i < j:
    #            c=[]
    #            c1=[]
    #            for n in range(len(gclones)):
    #                c.append(epi_values.get((i,j),0)+epi_values.get((j,i),0))# * ((gclones[n][i] == gclones[n][j]) * 2 - 1))
    #                c1.append((gclones[n][i] == gclones[n][j]) * 2 - 1)
                    #print i,j,c[-1],c1[-1]#np.mean(c),np.var(c),np.mean(c1),np.var(c1)
    #exit()
    pfit = pop.get_fitness_statistics()
    for n in range(len(gclones)):
        fitness.append(pop.get_fitness(sclones[n]))
        #print sum((2*gclones[n]-1)*ff),sum(gclones[n]),fitness
        #additive_fitness.append(sum((2 * gclones[n] - 1) * ff))
        a=sum((2 * gclones[n] - 1) * ff)
        #epistatic_fitness.append(fitness - additive_fitness[-1])
        E = 0
        for i in snps:#epis:#range(L):#snps:
            for j in snps:#epis:#range(L):#snps:
                if j > i:
                    E += (epi_values.get((i,j),0)+epi_values.get((j,i),0)) * ((gclones[n][i] == gclones[n][j]) * 2 - 1)
                    #print ii,is_epistasis, n,E#(gclones[n][i]==gclones[n][j])*2-1
        snp_epistatic_fitness.append(E)

    #genealogy
    #print additive_fitness[:10],np.mean(additive_fitness)
        #print xxx, sigma2, L, pi, rec, alpha, is_epistasis, ii*gen_step, n, fitness[-1], np.var(fitness), E,sum(gclones[n]), fitness[-1]-a,0,0, pop_entropy(pop), 0,0#mafs, mafn
    #return np.mean(fitness), np.var(fitness), np.mean(snp_epistatic_fitness), s_number, n_number, np.mean(LD['s']), np.mean(LD['n']), np.mean(LD['sn']), np.mean(mafs['s']), np.mean(mafs['n'])
    return pfit.mean, pfit.variance, np.mean(snp_epistatic_fitness), s_number, n_number, LD, LDr2, np.mean(mafs['s']), np.mean(mafs['n'])


def pop_distance(pop):
    [gclones, nclones, sclones, snps, snps_cols] = pop_snps(pop, sample)

    LD = {'s': {},'n': {},'sn': {}}
    for i in snps:
        for j in snps:
            if j > i:
                dist = ((j - i) / 10) * 10
                if i not in epis and j not in epis:
                    t = 's'
                elif i in epis and j in epis:
                    t = 'n'
                else:
                    t = 'sn'
                LD[t][dist] = LD[t].get(dist, [])
                LD[t][dist].append(fn.count_r2(snps_cols[i],snps_cols[j]))
    return LD


def pop_entropy(pop):
    S = 0
    for i in range(L):
        freq = pop.get_allele_frequency(i)
        if freq not in [0, 1]:
           S -= freq * np.log(freq) + (1.0 - freq) * np.log(1.0 - freq)
            #S1 = 2.0*N/(float(N)-1.0) * freq * (1.0-freq)
        #S.append(S1)
    #S = np.sum(S)
    return S


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
    #for i in snps_cols:
            snps_cols1[i] = [snps_cols[i][int(x)] for x in Z]
    return snps_cols1


def pop_fitness(piop):
    [gclones, nclones, sclones, snps, snps_cols] = pop_snps(pop, sample)
    efitness = []
    efitness_random = []
    for n in range(len(gclones)):
        dif_snps = {}
        genotype = []
        k = 0
        for i in range(len(gclones[n])):
            if i in L_dif and i in snps:
                genotype.append(gclones[n][i])
                dif_snps[i] = k
                k += 1
        E = 0
        for i in dif_snps:
            for j in dif_snps:
                if j > i:
                    E += (epi_values.get((i,j),0)+epi_values.get((j,i),0)) * ((genotype[dif_snps[i]] == genotype[dif_snps[j]]) * 2 - 1)
        random_genotype = genotype[:]
        np.random.shuffle(random_genotype)
        E_random = 0
        for i in dif_snps:
            for j in dif_snps:
                if j > i:
                    E_random += (epi_values.get((i,j),0)+epi_values.get((j,i),0)) * ((random_genotype[dif_snps[i]] == random_genotype[dif_snps[j]]) * 2 - 1)
        efitness.append(E)
        efitness_random.append(E_random)
    return efitness, efitness_random, len(dif_snps)


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
        [fita, fite, fitsnps, s, n, LD, LDr2, mafs, mafn] = pop_stats(pop,sample)
        print xxx, sigma, L, ppp, mu,rec, alpha, is_epistasis, ii*gen_step, fita, fite, fitsnps, s, n, np.mean(LD['s']), np.mean(LD['n']),np.mean(LD['sn']),np.mean(LDr2['s']),np.mean(LDr2['n']),np.mean(LDr2['sn']), pop_entropy(pop), mafs, mafn
    elif args.fmt == 'alignment':
        snps_cols = pop_aln(pop)
        snps = snps_cols.keys()
        snps.sort()
        for i_rank, i in enumerate(snps):
            for g in range(len(snps_cols[i])):
                print xxx, sigma0, L, ppp, mu,rec, alpha, is_epistasis, ii*gen_step, i, i_rank, g, snps_cols[i][g], int(i in epis)
    elif args.fmt == 'fitness':
        [E, E_random, nsnps] = pop_fitness(pop)
        for g in range(len(E)):
            print xxx, sigma0, L, ppp, mu,rec, alpha, is_epistasis, ii*gen_step, E[g], E_random[g], nsnps
    elif args.fmt == 'distance':
        LD = pop_distance(pop)
        for allele_type in LD:
            for distance in LD[allele_type]:
                print(xxx, sigma0, L, N, mu,rec, alpha, is_epistasis, ii*gen_step, allele_type, distance, np.mean(LD[allele_type][distance]), len(LD[allele_type][distance]))
    elif args.fmt == 'arg':
        rec_events = pop_arg(pop)
        for xx in range(1):
            for x in rec_events:
                print xxx, sigma0, L, ppp, mu,rec, alpha, is_epistasis, ii*gen_step, x

parser = argparse.ArgumentParser()

parser.add_argument('-i', type=int, action='store', dest='i', help="iterations")
parser.add_argument('-N', action='store', dest='N', help="population size(s)")
parser.add_argument('-L', type=int, action='store', dest='L', help="sequence length")
parser.add_argument('--mu', action='store', dest='mu', help="mutation rate(s) population scaled")
parser.add_argument('--alpha', action='store', dest='alpha', help="epistatic variance(s) from 0 to 1")
parser.add_argument('--rec', action='store', dest='rec', help="recombination rate(s) population scaled")
parser.add_argument('--gen', action='store', dest='gen', help="generations to simulate, number of steps")
parser.add_argument('--sigma', type=float,action='store', dest='sigma', help="fitness variation")
parser.add_argument('--sample', type=int,action='store', dest='sample', help="sample size for LD calculations")
parser.add_argument('--freq', type=float,action='store', dest='freq', help="initial allele frequencies")
parser.add_argument('--pi', type=float, action='store', dest='pi', help="genetic diversity")
parser.add_argument('--entropy', type=float,action='store', dest='entropy', help="final allelic entropy (if 0 - no limit)")
parser.add_argument('--epi', action='store', dest='epi', help="type of epistasis",choices=['random','pairwise','positive'])
parser.add_argument('--outfmt', action='store', dest='fmt', help="output",choices=['stats','alignment','fitness','distance','arg'])
args = parser.parse_args()

L = args.L
alpha_list = [float(x) for x in args.alpha.split(',')]
rec_list = [float(x) for x in args.rec.split(',')]
mu_list = [float(x) for x in args.mu.split(',')]
entropy = args.entropy
freq = args.freq
pi1 = args.pi
sample = args.sample
sigma0 = args.sigma
alpha=alpha_list[0]
#sigma0 = 1.0
#mu0 = args.mu
#N = args.N
#pi_list = [float(x) for x in args.pi.split(',')]
N_list = [int(x) for x in args.N.split(',')]

[gens, gen_ints] = [int(x) for x in args.gen.split(',')]
gen_step = gens / gen_ints
if entropy > 0:
    gen_step = 1
    gen_ints = 100000

#if args.fmt == 'stats':
#    print 'xxx', 'sigma','L','N', 'mu','rec', 'alpha', 'epi', 'gen', 'fita', 'fite','fitsnps', 's', 'n', 'LDs','LDn', 'LDsn', 'r2s','r2n','r2sn','entropy', 'mafs','mafn'
#elif args.fmt == 'alignment':
#    print 'xxx', 'sigma','L','N', 'mu','rec', 'alpha', 'epi', 'gen', 'position', 'position_rank', 'genotype', 'allele', 'allele_type'
#elif args.fmt == 'fitness':
#    print 'xxx', 'sigma','L','N', 'mu','rec', 'alpha', 'epi', 'gen', 'E', 'E_random', 'n'
#elif args.fmt == 'distance':
#    print('xxx', 'sigma','L','N', 'mu','rec', 'alpha', 'epi', 'gen', 'allele_type', 'distance', 'LD', 'n')
#elif args.fmt == 'arg':
#    print 'xxx', 'sigma','L','N', 'mu','rec', 'alpha', 'epi', 'gen', 'events'

for jjj in range(args.i):
    for rec0 in rec_list:
        for mu0 in mu_list:
            N = N_list[0]
            for ppp in [1.0,2.0,5.0]:
            #for N in N_list:
                #pi = 0.5 * (1.0 - np.sqrt(1.0 - 2.0 * pi1))
                pi0 = 1
                pi = 0.5
                #L = 102
                #pi0=0.5
                #pi = 0.5
                #pi0 = 0.5
                #pi=pi0
                #pi0=1.0
                #L = int(round(max(2,N*2*mu0 * 1000)))
                mu = mu0  #/ float(N) * 10000.0# normalisation
                #mu = N
                #N=10000
                #sigma0=0.005
                #N=1000
                rec = rec0 /float(N) # normalisation
                if args.fmt == 'fitness':
                    epi_list = [1.0]
                else:
                    epi_list = [0.0, 1.0]
                for is_epistasis in epi_list:
                    #sys.stderr.write('# ' + str(pi0) + ' ' + str(rec) + ' ' + str(alpha) + ' ' + str(is_epistasis) + '\n')
                    xxx = np.random.random(1)[0]
                    #L = int(round(pi0*L0*50))
                    ### initiate population
                    epis = list(range(1, L, 3)) + list(range(2, L, 3))  # sites under selection
                    #epis = range(L)
                    L1 = len(epis)
                    np.random.shuffle(epis)
                    L_dif = epis#[:len(epis)/2]
                    L1 = 0
                    for xxx1 in epis:
                        if xxx1 in L_dif:
                            L1 += 1
                    epi_values = {}
                    #L1 = float(L1) * ppp
                    #print L1
                    #sigma2 = pow(sigma0 / 1,2) *len(L_dif) * 10
                    #sigma2 = sigma0 #pow(sigma0 / 1,2) #*len(L_dif)
                    sigma2 = pow(sigma0/float(N), 2) * float(L1)#*float(L1-1) #* 4.0 * pi * (1.0 - pi)
                    Va = sigma2 * (1 - alpha * is_epistasis)
                    Vi = sigma2 * alpha * is_epistasis
                    #L_dif = epis
                    Va=sigma2
                    sigma=args.sigma
                    ff = []
                    #print sigma2, Va, Vi, np.sqrt(Va/float(L1)), L1
                    for i in range(L):
                        if i in epis and i in L_dif:
                            if i in L_dif:
                                if Va:
                                    #ff.append(-np.sqrt(Va/float(L1)))#/(4.0*pi*(1.0-pi))))#/float(N))
                                    ff.append(-1 * np.random.gamma(np.sqrt(Va/float(L1))*N,1.0,1)[0]/float(N))
                                else:
                                    ff.append(0)#-np.sqrt(sigma2))
                        else:
                            ff.append(0)
                    sys.stderr.write('# ' + str( len(L_dif)) + ' '+ str( is_epistasis)+' ' + str(Va)+' ' + str(Vi) + ' ' + str(np.sqrt(Va/float(L1))) + '\n')
                    epis1 = epis[:]
                    np.random.shuffle(epis1)
                    epis2 = epis1[:len(epis1)/2]
                    epis1 = epis1[len(epis1)/2:]
                    if args.epi == 'random':
                        pop.set_random_epistasis(np.sqrt(Vi))# * float(len(L_dif))))
                    elif args.epi == 'pairwise' and Vi > 0:
                        for kkk in range(len(epis1)):
                            if len(epis2) > kkk:
                                i = epis1[kkk]
                                j = epis2[kkk]
                                if j < i: [i,j] = [j,i]
                        #for i in epis:
                            #for j in epis:
                                if i in L_dif and j in L_dif and j > i and True: #np.random.random(1)[0] < ppp: #i in L_dif:
                                    #value = np.sqrt(Vi/float(L1))/float(L1-1)*float(N)
                                    #value = np.random.gamma(np.sqrt(Vi/float(L1))/float(L1-1)*N,1.0,1)[0]
                                    #value = np.random.normal(0,np.sqrt(Vi/float(L1))/float(L1-1),1)[0]
                                    value = np.random.normal(0,np.sqrt(Vi/float(L1)) * ppp,1)[0]
                                    #value = np.sqrt(Vi/float(L1))/float(L1-1)
                                    #value = value / float(N)
                                    epi_values[(i,j)] = value
                    pop1 = h.haploid_highd(L)
                    pop1.carrying_capacity = N
                    pop1.set_wildtype(N)
                    pop1.recombination_model = h.CROSSOVERS
                    pop1.mutation_rate = mu       # mutation rate per site per generation
                    pop1.outcrossing_rate = 1    # probability of sexual reproduction per gen
                    pop1.crossover_rate = rec      # probability of crossover per site per gen
                    pop1.set_fitness_additive(ff)
                    for pair in epi_values:
                        i = pair[0]
                        j = pair[1]
                        pop1.add_fitness_coefficient(epi_values[(i,j)], [i,j])
                    pop2 = h.haploid_highd(L)
                    pop2.carrying_capacity = N
                    pop2.set_wildtype(N)
                    pop2.recombination_model = h.CROSSOVERS
                    pop2.mutation_rate = mu       # mutation rate per site per generation
                    pop2.outcrossing_rate = 1    # probability of sexual reproduction per gen
                    pop2.crossover_rate = rec      # probability of crossover per site per gen
                    pop2.set_fitness_additive(ff)
                    for pair in epi_values:
                        i = pair[0]
                        j = pair[1]
                        pop2.add_fitness_coefficient(epi_values[(i,j)], [i,j])
                    pop1.evolve(0)
                    pop1.evolve(N*2)
                    pop2.evolve(0)
                    pop2.evolve(N*2)


                    [gclones, nclones, sclones, snps1, snps_cols1] = pop_snps(pop1, sample)
                    [gclones, nclones, sclones, snps2, snps_cols2] = pop_snps(pop2, sample)

                    snps = list(set(snps1).intersection(set(snps2)))
                    for snp1 in snps:
                        if snp1 in epis:
                            t1 = 1
                        else:
                            t1 = 0
                        for snp2 in snps:
                            if snp1 in epis:
                                t2 = 1
                            else:
                                t2 = 0
                            if snp2 > snp1 and t1 == t2:
                                r2_1 = fn.count_r2(snps_cols1[snp1],snps_cols1[snp2])
                                r2_2 = fn.count_r2(snps_cols2[snp1],snps_cols2[snp2])
                                print xxx, sigma, L, ppp, mu,rec, alpha, is_epistasis, t1, snp1,snp2,r2_1,r2_2
