\name{alleleSummary}
\alias{alleleSummary}
\title{
Allele sharing summary
}
\description{
This function summarises the allele flow for selected individuals, for a single genome simulation.
}
\usage{
alleleSummary(x, ids, ibd.status = FALSE)
}
\arguments{
  \item{x}{
A single genome simulation, i.e. a list of simulated chromosomes. Each chromosome is a list, with one entry for each individual. Each of these entries is a list of two matrices (one for each strand). The matrices have 2 columns (start position; allele) and one row for each segment unbroken by recombination.
}
  \item{ids}{
A vector of numerical IDs. If missing, all individuals are included.
}
  \item{ibd.status}{
A logical, meaningful only if \code{length(ids)==2}. If TRUE the IBD status (number of alleles shared IBD, either 0 1 or 2) of each segment is computed, as well as the breakdown of their parental origin.
}
}
\details{
This function is useful for downstream analysis of simulations produced by \code{\link{IBDsim}}.
}
\value{
A numerical matrix. Each row corresponds to a chromosomal segment. The first 4 columns describe the segment (chromosome, start, end, length), and are followed by two columns (paternal allele, maternal allele) for each of the selected individuals. If \code{ibd.status=TRUE} five more columns are added: \code{ibd}, \code{ibd_pp}, \code{ibd_pm}, \code{ibd_mp} and \code{ibd_mm}. The first of these indicate the IBD status (0, 1 or 2) in the segment, while the latter 4 give the parental breakdown of this number. For instance, \code{ibd_pm} is 1 if the _p_aternal allele of the first individual is IBD with the _m_aternal allele of the second individual, and 0 otherwise.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{IBDsim}}
}
\examples{
### Sibling simulation (3 sims of chromosomes 1 and 2)
x = nuclearPed(2)
sim = IBDsim(x, sims=3, chromosomes=1:2)

alleleSummary(sim[[1]]) # First sim, summary of all individuals
alleleSummary(sim[[1]], ids=3:4) # Summary of the siblings
alleleSummary(sim[[1]], ids=3:4, ibd.status=TRUE) # IBD breakdown of the siblings

# Trivial example: Summary of the father.
# Being the first founder, his alleles are denoted 1 and 2 everywhere.
fa = alleleSummary(sim[[1]], ids=1) 
stopifnot(all(fa[,c('1p', '1m')]==c(1,1,2,2)))

}
\keyword{ math }

