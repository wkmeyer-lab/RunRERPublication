j=$1
echo $j
useAlts=$2

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$j --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunRERAndCorrelation.R r=$j s=g m=Data/zoonomiaAllMammalsTrees.rds l=170

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$j --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunMergedAnalysis.R r=$j s=g m=Data/zoonomiaAllMammalsTrees.rds l=170

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$j --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT EnrichmentAnalysisNew.R r=$j m=Data/KeggReactome.gmt p=F

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$j --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CategoricalCombinePairwiseResults.R r=$j g=gene

sbatch -t 3300 --job-name=MainRerMergeEnrichCombineCombine$j --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CategoricalCombinePairwiseResults.R r=$j g=KeggReactome



for i in `seq 1 100`; do sbatch -t 3300 --job-name=Alternates$j$i --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT ProcessAlternates.R r=$j s=g m=Data/zoonomiaAllMammalsTrees.rds l=170 i=$i; done

for i in `seq 1 100`; do sbatch -t 3300 --job-name=Alternates$j$i --dependency singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunAlternateEnrichments.R r=$j m=Data/KeggReactome.gmt p=F i=$i; done

deps=$(squeue -h -o "%A %j" | awk '$2 ~ /Alternates/ {print $1}' | paste -sd:)

sbatch -t 3300 --job-name=CombineAlts$j --dependency=afterok:$deps -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CombineAlternates.R r=$j g=gene

deps2=$(squeue -h -o "%A %j" | awk '$2 ~ /Combine/ {print $1}' | paste -sd:)

sbatch -t 3300 --job-name=GeneAltsToCombo$j --dependency=afterok:$deps2 -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT AddAlternatesToCombinations.R r=$j g=gene

sbatch -t 3300 --job-name=GoAltsToCombo$j --dependency=afterok:$deps2 -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT AddAlternatesToCombinations.R r=$j g=KeggReactome

