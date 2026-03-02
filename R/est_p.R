#' Estimate the size of an outbreak from a phylogenetic tree.
#'
#' @param mu A numeric value for the per-generation substitution rate,
#'  with the units substitutions per site per generation.
#' @param Tp  A numeric value for the length of the phylogeny in substitutions per site.
#' @param S An integer value of the number of sequences/tips in the phylogeny.
#'
#' @returns A number. Estimated value of p.
#' @export est_p


est_p <- function(mu, Tp, S){
  rho <- (-Tp + sqrt(Tp^2+4*mu^2*S))/(2*mu)
  return (rho^2)
}
