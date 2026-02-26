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
