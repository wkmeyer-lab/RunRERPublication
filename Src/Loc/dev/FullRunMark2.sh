

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$1 --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunRERAndCorrelation.R "$@"

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$1 --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunMergedAnalysis.R "$@"

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$1 --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT EnrichmentAnalysisNew.R "$@"

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$1 --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CategoricalCombinePairwiseResults.R "$@"



for i in `seq 1 2`; do sbatch -t 3300 --job-name=Alternates$1$i --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT ProcessAlternates.R "$@" i=$i; done

for i in `seq 1 2`; do sbatch -t 3300 --job-name=Alternates$1$i --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunAlternateEnrichments.R "$@" i=$i; done

deps=$(squeue -h -o "%A %j" | awk '$2 ~ /Alternates/ {print $1}' | paste -sd:)

sbatch -t 3300 --job-name=CombineAlts$1 --dependency=afterok:$deps -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CombineAlternates.R "$@"

deps2=$(squeue -h -o "%A %j" | awk '$2 ~ /Combine/ {print $1}' | paste -sd:)

sbatch -t 3300 --job-name=GeneAltsToCombo$1 --dependency=afterok:$deps2 -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT AddAlternatesToCombinations.R"$@"

