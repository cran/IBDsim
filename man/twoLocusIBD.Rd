\name{twoLocusIBD}
\alias{twoLocusIBD}

\title{
Joint IBD probabilities for two linked loci
}
\description{
Estimates by simulation the IBD matrix associated to a pair of non-inbred individuals at two linked loci. The entries of the matrix are the joint IBD probabilities 
\deqn{\kappa(i,j) = Pr(sharing i alleles IBD at locus 1, and j alleles IBD at locus 2).}
The loci can be autosomal (in which case \eqn{0 \le i,j \le 2}) or X-linked (in which case \eqn{0 \le i,j \le 2} if both individuals are female; 
otherwise \eqn{0 \le i,j \le 1}).  
}
\usage{
twoLocusIBD(x, ind1, ind2, rho=NULL, cM=NULL, Nsim, Xchrom=FALSE, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A pedigree in the form of a \code{\link{linkdat}} object.
}
  \item{ind1, ind2}{
Numeric ID labels of the two individuals.
}
  \item{rho}{
A number in the interval [0, 0.5]: the recombination fraction between the two loci. Should be NULL if \code{cM} is non-NULL.
}
  \item{cM}{
A non-negative number: the distance in centiMorgan between the two loci. Should be NULL if \code{rho} is non-NULL.
}
  \item{Nsim}{
The number of simulations to be performed.
}
  \item{Xchrom}{
A logical indicating if the markers are X-linked (if TRUE) or autosomal (FALSE).
}
  \item{verbose}{
A logical.
}
  \item{\dots}{
Further arguments to be passed on to \code{\link{IBDsim}}.
}
}
\details{
If the recombination fraction \code{rho} is given, it is converted to centiMorgan using Haldanes map function: \code{cM = -50*log(1-2*rho)}.
}
\value{
A numerical matrix. The entry in row \eqn{a} and column \eqn{b} is the estimate of \eqn{\kappa(a-1, b-1)} defined above. 
}

\author{
Magnus Dehli Vigeland
}

\examples{

#Example 1
x <- nuclearPed(2) 
rho <- 0.1 
cM <- -50*log(1-2*rho)
Nsim <- 100 # Should be increased to at least 1000 
r1 <- twoLocusIBD(x, ind1=3, ind2=4, rho=rho, Nsim=Nsim, seed=17)
r2 <- twoLocusIBD(x, ind1=3, ind2=4, cM=cM, Nsim=Nsim, seed=17)
stopifnot(identical(r1,r2))

# X chromosome
r1 <- twoLocusIBD(x, ind1=3, ind2=4, rho=rho, Nsim=Nsim, seed=17, Xchrom=TRUE)
r2 <- twoLocusIBD(x, ind1=3, ind2=4, cM=cM, Nsim=Nsim, seed=17, Xchrom=TRUE )
stopifnot(identical(r1,r2))

# Example 2
# Autosomal testing, example p. 182-183 Egeland, Kling and Mostad(2016, Book).
# Each simulations followed by exact counterpart.
x <- addSon(addSon(nuclearPed(2,1:2), 4), 5)
plot(x)
rho <- 0.25
Nsim <- 100 # Should be increased to at least 10000
twoLocusIBD(x, 1, 6, rho=rho, Nsim=Nsim, verbose=FALSE)[2,2];.5*(1-rho)
twoLocusIBD(x, 8, 6, rho=rho, Nsim=Nsim, verbose=FALSE)[2,2];.5*(rho^2+(1-rho)^2)
twoLocusIBD(x, 3, 6, rho=rho, Nsim=Nsim, verbose=FALSE)[2,2];.5*((1-rho)*(rho^2+(1-rho)^2) + rho/2)

# Example 3
# X chromosome testing, example granddaughter vs maternal grandfather.
y <- addDaughter(nuclearPed(1, sex=2), 3)
plot(y)
rho <- 0.25
Nsim <-100
twoLocusIBD(y, 1, 5, rho=rho, Nsim=Nsim, Xchrom=TRUE)
matrix(c(1-rho, rho, rho, 1-rho)/2, ncol=2) #exact

# Example 4 
# X chromosomes can help in incest cases: Incest by brother or father?
# Could alternatively consider one marker and use IBDsim:::oneLocusIBD
x <- nuclearPed(2, sex=c(1,2))
x <- addOffspring(x, father=3, mother=4, noff=1, sex=2)
par(mfcol=c(1,2))
plot(x, title="Brother")
Nsim <- 100
r1 <- twoLocusIBD(x, ind1=2, ind2=5, rho=0.2,  Nsim=Nsim, Xchrom=TRUE)
stopifnot(r1[1,1]==0) # grand mother must share IBD 

#  Incest by father
x <- nuclearPed(1,sex=2)
x <- addOffspring(x, father=1, mother=3, noff=1, sex=2)
plot(x, title="Father")
r2 <- twoLocusIBD(x, ind1=2, ind2=4, rho=0.2,  Nsim=Nsim, Xchrom=TRUE)
stopifnot(r2[1,1]>0) # grand mother need not share IBD

}
\keyword{ math }