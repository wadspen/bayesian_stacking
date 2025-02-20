
for i in $(seq 1 10);
do
	sbatch dyn_forc2.sh "$i"
done

