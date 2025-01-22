# Loop through each pair of populations
for pair in "CCW CCE" "CCW CO" "CCE CO"; do
  # Split the pair into P1 and P2
  read P1 P2 <<<$(echo $pair)

  # Loop through each iteration
  for IT in {1..10}; do
    # Submit the job with sbatch
    sbatch -J "CROSSCOAL_${P1}_${P2}_${IT}" ~/merondun/cuculus_migration/msmc/4.Crosscoalescent_Iterative.sh ${P1} ${P2} ${IT}
  done
done
