module load imp
module load python3/pyrmsd

export cl=1
export st=0

cp ../A_models_clust${cl}_${st}.txt ../scoresA.txt
cp ../B_models_clust${cl}_${st}.txt ../scoresB.txt
cp ../../density.txt density.txt
nohup python /home/ignacia/SOFTW/imp-sampcon_2022/pyext/src/exhaust.py \
       -n BORCS6 -p ../ -ra A_models_clust${cl}_${st}.rmf3 -rb B_models_clust${cl}_${st}.rmf3 -d density.txt \
       -m cpu_omp -c 8 -g 2.0 --align > clustering.log &
