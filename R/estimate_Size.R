#' Estimate the size of an outbreak from a phylogenetic tree.
#'
#' @param tree A nucleotide-scaled phylogenetic tree as an object of class "phylo".
#' @param MR A numeric value for the per-generation substitution rate,
#'  with the units substitutions per site per generation, OR a vector of per-generation
#'  substitution rate estimates (i.e., substitution rate posterior estimates from a BEAST log
#'  file multiplied by the generation interval)
#' @param adjusted Boolean. True by default, use False if you do not want your estimates
#'  to be adjusted based on simulated accuracy data. WARNING: This may make your estimate
#'  less accurate! If you turn off adjustments and use only a point value of the substitution
#'  rate, you will receive a point value estimate of p and N without credible intervals.
#'
#' @returns A named vector. 'tree length' is the sum of the branch lengths of your tree;
#'  'mean_estimate_p' is the estimate for the proportion of cases sequenced, adjusted
#'  using our simulated accuracy data; 'mean_estimate_N' is the estimate of the outbreak
#'  size calculated from this proportion; 'upper' and 'lower' values give 95\% credible intervals.
#' @examples
#'   require(ape)
#'   tree = poser::example_tree
#'   per_gen = 0.0002 * (28/365) #0.0002 subs/site/year, 28 day generation interval
#'   estimate_Size(tree, per_gen)
#'   #output:
#'   #tree_length  mean_estimate_p lower_estimate_p upper_estimate_p  mean_estimate_N lower_estimate_N upper_estimate_N
#'   #9.730849e-03     8.807112e-03     3.006062e-03     1.460816e-02     6.699132e+03     4.038838e+03  1.962701e+04
#'
#' @export estimate_Size

estimate_Size = function(tree, MR, adjusted = T){
  if(adjusted == T){
    if(length(MR) == 1){
      treelength = sum(tree$edge.length)
      tips = length(tree$tip.label)

      prop_ests = poser::est_p(MR, treelength, tips)

      #bring in uncertainty from method itself
      gtfit <- glmmTMB::glmmTMB((1/accuracy_ratio)~splines::bs(estimate_seqd, 4)
                                , dispformula = ~splines::bs(estimate_seqd, 4)
                                , data = accuracy_data
      )

      df = data.frame(estimate_seqd = prop_ests)
      df$error_sd = predict(gtfit, df, type = "disp")
      df$error_mean = predict(gtfit, df)

      ci_2.5_50_97.5 = qnorm(c(0.025, 0.5, 0.975), mean = df$error_mean,
                             sd = df$error_sd)
      p_2.5_50_97.5 = df$estimate_seqd*ci_2.5_50_97.5
      N_2.5_50_97.5 = tips/p_2.5_50_97.5

      df2 = c(tree_length = treelength,
              mean_estimate_p = unname(p_2.5_50_97.5[2]),
              lower_estimate_p = unname(p_2.5_50_97.5[1]),
              upper_estimate_p = unname(p_2.5_50_97.5[3]),
              mean_estimate_N = unname(N_2.5_50_97.5[2]),
              lower_estimate_N = unname(N_2.5_50_97.5[3]),
              upper_estimate_N = unname(N_2.5_50_97.5[1]))

      return(df2)
    }
    else{
      if(length(MR) <= 500){
        warning("The supplied substitution rate is a short vector;
              is this correct?")
      }
      treelength = sum(tree$edge.length)
      tips = length(tree$tip.label)
      gtfit <- glmmTMB::glmmTMB((1/accuracy_ratio)~splines::bs(estimate_seqd, 4)
                                , dispformula = ~splines::bs(estimate_seqd, 4)
                                , data = accuracy_data
      )

      prop_ests = c()
      for(i in 1:length(MR)){
        prop_ests[i] = poser::est_p(MR[i], treelength, tips)
      }

      df = data.frame(estimate_seqd = prop_ests)
      df$error_sd = predict(gtfit, df, type = "disp")
      df$error_mean = predict(gtfit, df)

      combined_dists = c()
      for(j in 1:nrow(df)){
        combined_dists[j] = df$estimate_seqd[j]*rnorm(1, mean = df$error_mean[j],
                                                      sd = df$error_sd[j])
      }

      p_2.5_50_97.5 = quantile(combined_dists, c(0.025, 0.5, 0.975))
      N_2.5_50_97.5 = tips/p_2.5_50_97.5

      df2 = data.frame(tree_length = treelength,
                       mean_estimate_p = unname(p_2.5_50_97.5[2]),
                       lower_estimate_p = unname(p_2.5_50_97.5[1]),
                       upper_estimate_p = unname(p_2.5_50_97.5[3]),
                       mean_estimate_N = unname(N_2.5_50_97.5[2]),
                       lower_estimate_N = unname(N_2.5_50_97.5[3]),
                       upper_estimate_N = unname(N_2.5_50_97.5[1]))

      return(df2)
    }
  }
  else{
    if(length(MR) == 1){
      treelength = sum(tree$edge.length)
      tips = length(tree$tip.label)

      prop_ests = poser::est_p(MR, treelength, tips)

      size_est = tips/prop_ests

      df2 = c(tree_length = treelength,
              estimate_p = prop_ests,
              estimate_N = size_est)

      return(df2)
    }
    else{
      if(length(MR) <= 500){
        warning("The supplied substitution rate is a short vector;
              is this correct?")
      }
      treelength = sum(tree$edge.length)
      tips = length(tree$tip.label)

      prop_ests = c()
      for(i in 1:length(MR)){
        prop_ests[i] = poser::est_p(MR[i], treelength, tips)
      }

      p_2.5_50_97.5 = quantile(prop_ests, c(0.025, 0.5, 0.975))
      N_2.5_50_97.5 = tips/p_2.5_50_97.5

      df2 = data.frame(tree_length = treelength,
                       mean_estimate_p = unname(p_2.5_50_97.5[2]),
                       lower_estimate_p = unname(p_2.5_50_97.5[1]),
                       upper_estimate_p = unname(p_2.5_50_97.5[3]),
                       mean_estimate_N = unname(N_2.5_50_97.5[2]),
                       lower_estimate_N = unname(N_2.5_50_97.5[3]),
                       upper_estimate_N = unname(N_2.5_50_97.5[1]))

      return(df2)
    }
  }
}
