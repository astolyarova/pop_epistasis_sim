
Scripts used to simulate evolution of populations with varying levels of genetic diversity with and without epistasis.

The script `run_simulations.sh` will reproduce simulations used in Fig. 1bc, Table 1 and Supplementary Fig. 14-19.
For large population size and mutation rate, simulations can take some time. To produce a lot of simulations used in the paper, the calculations were parallelized.

### Pairwise epistasis

We simulate pairwise compensatory epistasis using [FFPopSim](https://github.com/neherlab/ffpopsim) python package (Zanini and Neher, Bioinformatics 2012). 
We model two types of sites, depending on whether mutations in them are neutral (with selection coefficient s = 0) or weakly deleterious (s ≤ 0), 
representing synonymous and nonsynonymous sites correspondingly. There are twice as many nonsynonymous as synonymous sites. 
Under the non-epistatic model, selection coefficient of a nonsynonymous mutation is independent of the genetic background. 
Under pairwise epistasis, a mutation at one nonsynonymous site can be fully compensated by a mutation at another site.

The python script to perform such simulations is `sim_pairwise.py`.

```python
usage: sim_pairwise.py [-h] [-i I] [-N N] [-L L] [--mu MU]
                       [--epi_mode {isolated}] [--rec REC] [--gen GEN]
                       [--sigma SIGMA] [--sample SAMPLE] [--freq FREQ] [--p P]
                       [--mode {panmixic,admixture,balancing,sweep}]
                       [--outfmt {stats,alignment,fitness,distance,arg}]

optional arguments:
  -h, --help            show this help message and exit
  -i I                  number of simulations
  -N N                  population size
  -L L                  sequence length
  --mu MU               mutation rate(s)
  --epi_mode {isolated}
                        epstasis mode; only 'isolated' is implemented so far
  --rec REC             recombination rate(s) population scaled
  --gen GEN             generations to simulate, number of steps
  --sigma SIGMA         selection coefficient(s) of deleterious mutations
                        (population-scaled)
  --p P                 scale coefficient of DFE
  --sample SAMPLE       sample size
  --freq FREQ           frequency threshold
  --mode {panmixic,admixture,balancing,sweep}
                        population mode
  --outfmt {stats,alignment,fitness,distance,arg}
                        output

```


#####Detailed description:


If `-i` (the number of simulations) is equal to 1, a pair of simulations will be produced: one without epistasis and one with epistasis under the same population parameters
and simulation mode (see below). `I` means the number of such pairs of simulations. All simulations are independent from each other and are defined only with 
parameters of the population and the simulation mode.

`-N` indicates population size (the individuals are assumed to be haploid). Population size also defines the number of generations needed to be simulated to achieve mutation-selection equilibirum: we usually use 10N
generations.

`-L` parameters determines sequence length. 1/3 of these sites will be attributed as synonymous (neutral); 2/3 will be nonsynonymous (deleterious).

`--mu` is the mutation rate, measured as mutations per generation per site.

`--epi_mode` (is not a real parameter at the moment; should be "isolated")

`--rec` is population-scaled recombination rate.

`--gen` determines the simulation time, measured in generations, and the number of measurements output in the course of the simulation. The format is 
"TOTAL_GENERATIONS,STEPS". For example, `--gen 1000,1` means that the simulation will run for 1000 generations, the output will be given at generations 0 and 1000;
`--gen 10000,10` means that the simulation will run for 10000 generations and the output will be printed at generations 0,1000,2000,...,10000.

`--sigma` is the population-scaled (negative) selection coefficient of nonsynonymous mutations (`--sigma 1` means the fitness value of a single mutant is the fitness value of a wild type minus 1/N). If the `--p` argument is set to zero, every nonsynonymous mutation has the same selection coefficient; if not, sigma specifies their mean (see below).

`--p` defines the sparsity of fitness effects of nonsynonymous mutations. If set to zero, each nonsynonymous mutation has the selection coefficient of -sigma/N. If not, 
the selection coefficient values are sampled from gamma distribution with parameters mean=sigma/N and scale=p/N

`--sample` is the number of genotypes sampled from the population to calculate the output statistics (such as number of SNPs, LD, pi etc).

`--freq` is the low threshold for the frequency of SNPs used to calculate statistics; `--freq 0.05` means all SNPs with frequency less than 5% will be excludued.

`--mode` indicates the mode of a simulation. The options are:
 * `panmixic` -- the population is simulated starting from the homogeneous wild-type population, also having the most fit genotype, and then accumulates mutations, eventually reaching mutation-selection equilibrium. This is the *basic* simulation mode, assuming only negative selection, mutation and recombination acting on the population.
 * `balancing` -- in addition to the factors acting in the `panmixic` mode, one randomly chosen nonsynonymous site is assigned to be under frequency-dependent selection: it is positively selected at frequencies below 0.5, and negatively selected at frequencies above 0.5. If this mutation gets lost in the course of simulation due to genetic drift, it is re-introduced.
 * `sweep` -- the mode of ongoing selective sweep. First, the population is simulated for 10N generations in the presence of negative selection only (like in the `panmixic` mode). Than, a single positively selected mutation is introduced at a randomly chosen nonsynonymous site. If the mutation gets lost in the course of simulation, it is re-introduced. The simulation ends when this mutation achieves frequency 0.5 (or for 100000 generations, if it doesn't happen). The `--gen` parameter doesn't mean anything in this simulation mode.
 * `admixture` -- simulation of a sequence evolution after an admixture event. First, two separate populations are simulated under the same selection mode (either epistatic or non-epistatic) for 10N generations; for the epistatic mode, the same pairs of nonsynonymous sites are assumed to interact in both populations. After that, randomly sampled 50% of genotypes from each population are pooled to form a new population. This descendant population is then simulated for the number of generations stated in `--gen`. 

`--outfmt` defines the output format. There are several options available:
 * `stats` is the most basic mode; it returns multiple statistics calculated in a random population sample of size stated in `--sample`. These statistics are:
   * mean and variance of the sampled individuals' fitness
   * number of synonymous and nonsynonymous SNPs
   * linkage disequilibrium between all SNPs, between synonymous and nonsynonymous SNPs separately; is calculated as D and r<sup>2</sup>
   * average minor allele frequency (MAF) of synonymous and nonsynonymous SNPs
   * nucleotide diversity (the average number of differences between two genotypes); pn and ps
 * `distance` gives LD between pairs of SNPs of different type as a function of the distance between the SNPs, and a number of such pairs;
 * `fitness` returns statistics of the amount of mutational load reduced due to epistasis:
   * the average sum of epistatic coefficients between **polymorphic** positions (in the case of a non-epistatic simulation, it is calculated **as if** the compensatory epistasis existed; this can be used as a control of how often the "compensatory" mutations happen to co-occurr by chance);
   * the average sum of epistatic coefficients between polymorphic positions after random permutation of genotypes; can be again used as a control of a random distribution of alleles not caused by epistatic-driven linkage;
   * number of SNPs
   * the average (per genotype) mutation load compensated by epistasis
 * `alignment` gives the alignment-like representation of the sampled genotypes: for each genotype and each position, there is the allelic state of the position (either 0 or 1; 0 corresponds to the wild-type); the genotypes are sorted by pairwise identity.


The script can be easily adjusted for your purpose. The examples of how it could be run are in `run_simulations.sh` in the **pairwise epistasis** section.


### Global epistasis

In the global epistasis model  selection against each mutation depends only on the total number of such mutations present in this genotype. Under the global epistasis model, s depends on the total number of deleterious mutations in the genotype ndel. Specifically, s was set to equal 
`s(ndel) = -max[0, s0 * (1 - k * ndel)],`
where `s0` is the selection coefficient of a mutation on the wild-type background (set to -0.005), and k is the epistatic coefficient. Positive (antagonistic, or widening) epistasis between nonsynonymous mutations is modelled with `k > 0`, meaning that they become less deleterious in the context of other nonsynonymous mutations and become neutral if `ndel > 1/k`. Under negative (synergistic, or narrowing) epistasis (`k< 0`), the mutations become more deleterious in the presence of other mutations. Without epistasis `(k = 0)`, selection coefficient of each mutation equals s0. Global epistasis was simulated using [SLiM3.0](https://messerlab.org/slim/) (Haller and Messer 2019).

The `sim_global.slim` config file can be used to run the simulation. The parameters needed to be set are:
 * `_MUT_` - mutation rate (per generation per nucleotide)
 * `_REC_` - population-scaled recombination rate
 * `_S_` - selection coefficient of deleterious mutations (multiplied by N=2000)
 * `_K_` - epistatic coefficient *k*

The epistasis function and other parameters can be easily changed in the `.slim` script.


