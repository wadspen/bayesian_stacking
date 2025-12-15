
#for i in $(seq 1 30);
for i in $(seq 1 2);
do
  for j in ".8" ".9" ".95" ".98" "1"
  do
      sbatch dyn_forc2.sh "$i" "$j"
  done
done

