# POSER: Phylogenetic Outbreak Size Estimation in R
A simple method of estimating rabies outbreak sizes from phylogenetic trees in R. 

## Usage
To install the package:

```
devtools::install_github("RowanDurrant/poser")
```

Load a phylogenetic tree of sequences from an outbreak (substitution-scaled, not a time-scaled tree), and specify a per-generation substitution rate estimate for this outbreak:

```
tree = poser::example_tree #maximum likelihood tree- must be an object of class "phylo"
subrate = 0.0002 * (28/365) #TempEst, BEAST or literature estimate multiplied by generation interval in years

estimate_Size(tree, subrate)
```

Output:
```
tree_length  mean_estimate_p lower_estimate_p upper_estimate_p  mean_estimate_N lower_estimate_N upper_estimate_N 
0.009730849      0.009200249      0.002316431      0.036540946   6412.869606511   1614.627053190  25470.214009399

```
Where 'tree length' is the sum of branch lengths of your input tree, 'mean_estimate_p' is the mean estimate for the proportion of cases sequenced, 'mean_estimate_N' is the outbreak size estimate calculated using the mean value of p, and the 'upper' and 'lower' values are the 95% credible interval values.

Alternatively, you can use a BEAST log file for your substitution rate and the uncertainty will be built into the credible intervals:

```
subrate = poser::mock_BEAST_clockrates * (28/365) #when using your own log file, remember to remove the burn-in

estimate_Size(tree, subrate)
```

Output:
```
tree_length  mean_estimate_p lower_estimate_p upper_estimate_p  mean_estimate_N lower_estimate_N upper_estimate_N 
0.009730849      0.009052342      0.002179221      0.037889803   6517.650160933   1557.147186719  27073.892440975
```

## Troubleshooting

Sometimes you might get an estimated value of p considerably higher than 1. This is usually caused either by the input substitution rate estimate being too high, or by the outbreak being very small with most of the cases being sequenced, leading to lots of identical sequences in the dataset. In this case, the sequences/tree should be inspected to consider the possibility of a small outbreak or extreme sampling bias, and the substitution rate used should be scrutinised.
