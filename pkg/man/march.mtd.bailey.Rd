% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/march.AllGenerics.R
\name{march.mtd.bailey}
\alias{march.mtd.bailey}
\title{Bailey Confidence Intervals for a MTD model.}
\usage{
march.mtd.bailey(object, alpha)
}
\arguments{
\item{object}{the march.Model object on which compute the confidence intervals.}

\item{alpha}{the significance level.}
}
\value{
A list of half-length confidence intervals for each probability distribution of the MTD model.
}
\description{
Compute the confidence intervals using Bailey's formula on a march.Mtd
object. See Bailey BJR (1980) Large sample simultaneous confidence intervals
for the multinomial probabilities based ontransformation of the cell frequencies,
Technometrics 22:583–589, for details.
}
\examples{
# Compute the independence model for the pewee data.
Indep <- march.indep.construct(pewee)
# Display the model
print(Indep)
# Compute the half-length 95\% confidence interval for each element of the distribution.
march.indep.bailey(Indep,alpha=0.05)

# Compute a second-order MTDg model for the pewee data.
MTD2g <- march.mtd.construct(pewee,2,mtdg=TRUE)
# Display the model
print(MTD2g)
# Compute the half-length 95\% confidence interval for all parameters
# of the MTD2g model.
march.mtd.bailey(MTD2g,alpha=0.05)
}
\author{
Berchtold André
}
