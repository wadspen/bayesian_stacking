
for i in $(seq 1 50);
do
	sbatch dyn_forc2.sh "$i"
done

