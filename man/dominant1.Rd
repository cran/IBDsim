\name{dominant1}
\alias{dominant1}
\docType{data}
\title{
Autosomal dominant pedigree
}
\description{
A three generation pedigree affected with an autosomal dominant disorder.
}
\usage{data(twoloops)}
\format{
A data frame describing the pedigree in standard LINKAGE format (i.e. the columns are: Family ID, individual ID, father, mother, sex, affection status).
}
\examples{
data(dominant1)
x = linkdat(dominant1)
plot(x)
}
\keyword{datasets}
