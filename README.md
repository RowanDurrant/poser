# POSER: Phylogenetic Outbreak Size Estimation in R
A simple method of estimating rabies outbreak sizes from phylogenetic trees in R. 

To install the package:

```
devtools::install_github("RowanDurrant/poser")
```

Load a phylogenetic tree of sequences from an outbreak, and specify a per-generation substitution rate estimate for this outbreak:

```
library(ape)
tree = read.tree("rabies_outbreak.treefile")
sub_rate = 0.0002 * (28/365) #TempEst, BEAST or literature estimate multiplied by generation interval in years

estimate_Size(tree, subrate)
```

Alternatively, you can use a BEAST log file for your substitution rate and the uncertainty will be built into the credible intervals:

```
library(data.table)
BEAST_log = fread("rabies_outbreak_logfile.log.txt", 
                        select = "meanRate")$meanRate
BEAST_log_trimmed = BEAST_log[round(length(BEAST_log)*0.1):length(BEAST_log)] #remove burn-in
sub_rate = BEAST_log_trimmed * (28/365)

estimate_Size(tree, subrate)
```

Output:
```
tree_length  mean_estimate_p lower_estimate_p upper_estimate_p  mean_estimate_N lower_estimate_N 
      0.00072900       0.08380967       0.05903736       0.10858198     167.04516218     128.93483978 
upper_estimate_N 
    237.13796496
```
Where 'tree length' is the sum of branch lengths of your input tree, 'mean_estimate_p' is the mean estimate for the proportion of cases sequenced, 'mean_estimate_N' is the outbreak size estimate calculated using the mean value of p, and the 'upper' and 'lower' values are the 95% credible interval values.
