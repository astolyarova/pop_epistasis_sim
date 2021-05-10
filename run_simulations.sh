
# admixture model (Supplementary Fig. 17)

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 > epild_adm.out
python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 > epild_adm.out
python2.7 sim_pairwise.py -i 2000 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 > epild_adm.out
python2.7 sim_pairwise.py -i 5000 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 > epild_adm.out
python2.7 sim_pairwise.py -i 10000 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 > epild_adm.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 2000,100 --outfmt stats --p 0 > epild_adm_long.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt distance --p 0 > epild_adm_dist.out

# balancing selection model (Supplementary Fig. 15)

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 > epild_bal.out
python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 > epild_bal.out
python2.7 sim_pairwise.py -i 2000 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 > epild_bal.out
python2.7 sim_pairwise.py -i 5000 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 > epild_bal.out
python2.7 sim_pairwise.py -i 10000 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 > epild_bal.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,100 --outfmt stats --p 0 > epild_bal_long.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt distance --p 0 > epild_bal_dist.out


# negative selection only model (Supplementary Fig. 14)

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 > epild_pan.out
python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 > epild_pan.out
python2.7 sim_pairwise.py -i 2000 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 > epild_pan.out
python2.7 sim_pairwise.py -i 5000 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 > epild_pan.out
python2.7 sim_pairwise.py -i 10000 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 > epild_pan.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,100 --outfmt stats --p 0 > epild_pan_long.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 > epild_pan_dist.out


# ongoing selective sweep model (Supplementary Fig. 16)

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 > epild_sweep.out
python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 > epild_sweep.out
python2.7 sim_pairwise.py -i 2000 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 > epild_sweep.out
python2.7 sim_pairwise.py -i 5000 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 > epild_sweep.out
python2.7 sim_pairwise.py -i 10000 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 > epild_sweep.out

#python2.7 sim_pairwise.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 2000,10 --outfmt stats --p 0 > epild_sweep_long.out

python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt distance --p 100 > epild_sweep_dist.out


# mutation rate scaling vs population size scaling (Supplementary Fig. 19)

python2.7 sim_pairwise.py -i 100 -N 2000 -L 300 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 > mut_vs_Ne.out
python2.7 sim_pairwise.py -i 500 -N 2000 -L 300 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 | grep -v 'gen' >> mut_vs_Ne.out
python2.7 sim_pairwise.py -i 1000 -N 2000 -L 300 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 | grep -v 'gen'  >> mut_vs_Ne.out
python2.7 sim_pairwise.py -i 500 -N 4000 -L 300 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 40000,1 --outfmt stats --p 0| grep -v 'gen'  >> mut_vs_Ne.out
python2.7 sim_pairwise.py -i 100 -N 20000 -L 300 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 200000,1 --outfmt stats --p 0 | grep -v 'gen' >> mut_vs_Ne.out

# epistasis efficiency in populations of different diversity (Fig. 1bc)

python2.7 sim_pairwise.py -i 100 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 > aa_fit4.out
python2.7 sim_pairwise.py -i 100 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen' >> aa_fit4.out
python2.7 sim_pairwise.py -i 200 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen' >> aa_fit4.out
python2.7 sim_pairwise.py -i 500 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen'  >> aa_fit4.out
python2.7 sim_pairwise.py -i 1000 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen'  >> aa_fit4.out
python2.7 sim_pairwise.py -i 2000 -N 2000 -L 999 --mu 1e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen'  >> aa_fit4.out
python2.7 sim_pairwise.py -i 5000 -N 2000 -L 999 --mu 5e-7 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen'  >> aa_fit4.out
python2.7 sim_pairwise.py -i 10000 -N 2000 -L 999 --mu 2.5e-7 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 40000,1 --outfmt fitness --p 0 | grep -v 'gen'  >> aa_fit4.out




