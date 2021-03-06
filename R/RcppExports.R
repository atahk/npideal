# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

computeCombinations <- function(n) {
    .Call('npideal_computeCombinations', PACKAGE = 'npideal', n)
}

votetest <- function(v1, v2, focus = as.integer( c()), numLegisShuffles = 10L, numBootstrapSims = 10L) {
    .Call('npideal_votetest', PACKAGE = 'npideal', v1, v2, focus, numLegisShuffles, numBootstrapSims)
}

rankAgreeR <- function(v, i, j) {
    .Call('npideal_rankAgreeR', PACKAGE = 'npideal', v, i, j)
}

voteest <- function(v, right = 0L, left = 0L, maxLegisSamples = 10000L) {
    .Call('npideal_voteest', PACKAGE = 'npideal', v, right, left, maxLegisSamples)
}

