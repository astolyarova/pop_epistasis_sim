cd ~/schizo100/sim_ld

#rm schizorecr.out_$PBS_ARRAYID
#rm schizorec.out_$PBS_ARRAYID
#rm schizorects100.out_$PBS_ARRAYID
#rm schizorec100.out_$PBS_ARRAYID
rm sim2popdiv.out_${PBS_ARRAYID}
rm sim2popdiv1.out_${PBS_ARRAYID}

module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

for i in `seq 1 100`; do 
#while read -r gene rec; do
rec=`shuf -n 1 ../ldhat/genes_rho.out | awk '{print $2}'`
#python2.7 sim_2pop.py -i 1 -N 2000 -L 1000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 5 --sample 34,21 --freq 0 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2pop.out_${PBS_ARRAYID}
#python2.7 sim_2popdiv.py -i 1 -N 2000 -L 2000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 5 --sample 34,21 --freq 0.05 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2popdiv.out_${PBS_ARRAYID}
#python2.7 sim_2popdiv.py -i 1 -N 2000 -L 2000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 3 --sample 34,21 --freq 0.05 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2popdiv.out_${PBS_ARRAYID}
python2.7 sim_2popdiv.py -i 1 -N 2000 -L 2000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 10 --sample 34,21 --freq 0.05 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2popdiv.out_${PBS_ARRAYID}

#rec=0.05
#python2.7 sim_2popdiv.py -i 1 -N 2000 -L 2000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 5 --sample 34,21 --freq 0.05 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2popdiv1.out_${PBS_ARRAYID}
#python2.7 sim_2popdiv.py -i 1 -N 2000 -L 2000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 3 --sample 34,21 --freq 0.05 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2popdiv1.out_${PBS_ARRAYID}
#python2.7 sim_2popdiv.py -i 1 -N 2000 -L 2000 --mu 5e-5 --alpha 1.0 --rec $rec --sigma 10 --sample 34,21 --freq 0.05 --mode panmixic --gen 20000,1 --p 0 --outfmt stats >> sim2popdiv1.out_${PBS_ARRAYID}
done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

#rm adm_aa*.out_$PBS_ARRAYID


python2.7 pairwise_panmix1.py -i 100 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .0 --sigma 10 --sample 50 --freq 0.05 --mode admixture  --gen 2000,4 --outfmt distance --p 0 > dist.out ##adm_aa2.out_$PBS_ARRAYID
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm adm_aa_mut*.out_$PBS_ARRAYID

for i in `seq 1 10`; do

python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode admixture   --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode admixture    --gen 10000,10 --outfmt stats --p 0 | grep -v 'gen' >> adm_aa_mut3.out_$PBS_ARRAYID

done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm epild_adm*.out_$PBS_ARRAYID

for i in `seq 1 10`; do


python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 >> epild_adm.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 >> epild_adm.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 2 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 >> epild_adm.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 5 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 >> epild_adm.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 10 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt stats --p 0 >> epild_adm.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 2000,100 --outfmt stats --p 0 >> epild_adm_long.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode admixture --gen 200,1 --outfmt distance --p 0 >> epild_adm_dist.out_$PBS_ARRAYID

done
python sim_ld_dist2.py -i 1 -N 100000 -L 150 --mu 5e-7 --alpha 1.0 --rec 0.1,0.5,1.0 --sigma 3 --sample 100 --freq 0 --epi pairwise --gen 10000,10 --entropy 0 --outfmt alignment > aln.out
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm bal_aa_mut*.out_$PBS_ARRAYID

for i in `seq 1 10`; do

python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode balancing   --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode balancing    --gen 20000,10 --outfmt stats --p 10 | grep -v 'gen' >> bal_aa_mut3.out_$PBS_ARRAYID

done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm epild_bal*.out_$PBS_ARRAYID

for i in `seq 1 100`; do


python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 >> epild_bal.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 >> epild_bal.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 2 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 >> epild_bal.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 5 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 >> epild_bal.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 10 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt stats --p 0 >> epild_bal.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,100 --outfmt stats --p 0 >> epild_bal_long.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode balancing --gen 20000,1 --outfmt distance --p 0 >> epild_bal_dist.out_$PBS_ARRAYID

done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

#rm epild_genes.out #_$PBS_ARRAYID
#rm epild_rec.out_$PBS_ARRAYID
#rm epild_lowfreq.out_$PBS_ARRAYID
#rm epild_highfreq.out_$PBS_ARRAYID
rm epild_genes_st.out

for i in `seq 1 100`; do

#python2.7 pairwise_panmix0.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec 0.001 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_rec.out_$PBS_ARRAYID
#python2.7 pairwise_panmix0.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec 0.0 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_rec.out_$PBS_ARRAYID
#python2.7 pairwise_panmix0.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec 0.01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_rec.out_$PBS_ARRAYID
#python2.7 pairwise_panmixlow.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec 0.01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_lowfreq.out_$PBS_ARRAYID
#python2.7 pairwise_panmix0.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec 0.01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_highfreq.out_$PBS_ARRAYID
for i in `seq 1 10`; do
rec=`shuf -n 1 ../ldhat/genes_rho1.out | awk '{print $2}'`
#python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 10 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt distance --p 0 | grep -v 'gen' >> epild_genes.out #_$PBS_ARRAYID
python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 10 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 | grep -v 'gen' >> epild_genes_st.out #_$PBS_ARRAYID
python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 5 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 | grep -v 'gen' >> epild_genes_st.out #_$PBS_ARRAYID
python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 2 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 | grep -v 'gen' >> epild_genes_st.out #_$PBS_ARRAYID
python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 0 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 | grep -v 'gen' >> epild_genes_st.out #_$PBS_ARRAYID
#python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 0 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt distance --p 0 | grep -v 'gen' >> epild_genes.out #_$PBS_ARRAYID
#python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 1 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt distance --p 0 | grep -v 'gen' >> epild_genes.out #_$PBS_ARRAYID
#python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 2 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt distance --p 0 | grep -v 'gen' >> epild_genes.out #_$PBS_ARRAYID
#python2.7 pairwise_panmix_genes.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec $rec --sigma 5 --sample 32 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt distance --p 0 | grep -v 'gen' >> epild_genes.out #_$PBS_ARRAYID
done
done #done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm mut_vs_Ne*.out_$PBS_ARRAYID

for i in `seq 1 20`; do

python2.7 pairwise_panmix0.py -i 5 -N 2000 -L 300 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> mut_vs_Ne.out_${PBS_ARRAYID}
python2.7 pairwise_panmix0.py -i 5 -N 2000 -L 300 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> mut_vs_Ne.out_${PBS_ARRAYID}
python2.7 pairwise_panmix0.py -i 5 -N 2000 -L 300 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> mut_vs_Ne.out_${PBS_ARRAYID}
python2.7 pairwise_panmix0.py -i 5 -N 4000 -L 300 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 40000,1 --outfmt stats --p 0 >> mut_vs_Ne.out_${PBS_ARRAYID}
python2.7 pairwise_panmix0.py -i 5 -N 20000 -L 300 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 0 --sample 100 --freq 0.0 --mode panmixic --gen 200000,1 --outfmt stats --p 0 >> mut_vs_Ne.out_${PBS_ARRAYID}
#0
done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''


python2.7 pairwise_panmix.py -i 100 -N 2000 -L 300 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 2 --sample 50 --freq 0.05 --mode sweep --gen 10000,10  --outfmt stats|grep -v 'gen'  >> a.out
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm epild_pan*.out_$PBS_ARRAYID

for i in `seq 1 100`; do


python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 2 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 5 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 10 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,100 --outfmt stats --p 0 >> epild_pan_long.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_pan_dist.out_$PBS_ARRAYID

done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm epild_pan*.out_$PBS_ARRAYID

for i in `seq 1 100`; do


#python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
#python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
#python2.7 pairwise_panmix.py -i 2 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
#python2.7 pairwise_panmix.py -i 5 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID
#python2.7 pairwise_panmix.py -i 10 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt stats --p 0 >> epild_pan.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 1000,100 --outfmt stats --p 0 >> epild_pan_long.out_$PBS_ARRAYID

#python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .05 --sigma 10 --sample 100 --freq 0.05 --mode panmixic --gen 20000,1 --outfmt distance --p 0 >> epild_pan_dist.out_$PBS_ARRAYID

done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm sweep_aa_mut*.out_$PBS_ARRAYID

for i in `seq 1 10`; do

python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 1 -N 2000 -L 150 --mu 5e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 2 -N 2000 -L 150 --mu 1e-5 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 10 -N 2000 -L 150 --mu 5e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut3.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.0 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.03 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut1.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.05 --mode sweep   --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut2.out_$PBS_ARRAYID
python2.7 pairwise_panmix1.py -i 20 -N 2000 -L 150 --mu 1e-6 --alpha 1.0 --rec .01 --sigma 5 --sample 50 --freq 0.1 --mode sweep    --gen 200000,10 --outfmt stats --p 100 | grep -v 'gen' >> sweep_aa_mut3.out_$PBS_ARRAYID

done
cd ~/schizo100/sim_ld


module unload ScriptLang/python/3.8.3
module load ScriptLang/python/2.7
export PYTHONPATH=''

rm epild_sweep*.out_$PBS_ARRAYID

for i in `seq 1 100`; do


python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 >> epild_sweep.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 2.5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 >> epild_sweep.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 2 -N 2000 -L 999 --mu 1e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 >> epild_sweep.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 5 -N 2000 -L 999 --mu 5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 >> epild_sweep.out_$PBS_ARRAYID
python2.7 pairwise_panmix.py -i 10 -N 2000 -L 999 --mu 2.5e-6 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt stats --p 100 >> epild_sweep.out_$PBS_ARRAYID

#python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 2000,10 --outfmt stats --p 0 >> epild_sweep_long.out_$PBS_ARRAYID

python2.7 pairwise_panmix.py -i 1 -N 2000 -L 999 --mu 5e-5 --epi_mode isolated --rec .01 --sigma 10 --sample 100 --freq 0.05 --mode sweep --gen 500,1 --outfmt distance --p 100 >> epild_sweep_dist.out_$PBS_ARRAYID

done
