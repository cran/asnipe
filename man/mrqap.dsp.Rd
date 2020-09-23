\name{mrqap.dsp}
\alias{mrqap.dsp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MRQAP with Double-Semi-Partialing (DSP)
}
\description{
Calculate MRQAP with Double-Semi-Partialing (DSP) from Dekker et al (2007)
}
\usage{
mrqap.dsp(formula, intercept = TRUE, directed = "undirected",
	diagonal = FALSE, test.statistic = "t-value", 
	tol = 1e-07, randomisations = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
input formula (e.g. y ~ x1 + x2), where y and each x are NxN matrices
}
  \item{intercept}{
calculate intercept (TRUE or FALSE value)
}
  \item{directed}{
whether the network is directed or undirected (enter either "directed" or "undirected")
}
  \item{diagonal}{
whether to include self-loop values (TRUE or FALSE)
}
  \item{test.statistic}{
what to calculate P-value, either t-statistic ("t-value") or regression coefficient ("beta")
}
  \item{tol}{
tolerance value for the qr function
}
  \item{randomisations}{
number of randomisations to perform for calculating P-value.
}
}
\details{
Calculate the regression coefficient for each input matrix using the DSP method in Dekker et al (2007). This method randomises the residuals from the regression on each independent variable (fixed effect) in order to calculate the P value. This is the same as testing whether y is related to x1 on y while controlling for x2. This differs from regular mrqap, where the dependent (y) value is randomised, testing for whether y is related to x1 and x2 together.
}
\value{
Returns a mrqap.dsp object containing the regression coefficient and P-values for each indendent matrix (x) and associated statistics
}
\references{
Dekker, D., Krackhard, D., Snijders, T.A.B (2007) Sensitivity of MRQAP tests to collinearity and autocorellation conditions. \emph{Psychometrika} 72(4): 563-581.
}
\author{
Damien R. Farine
}

\examples{

library(asnipe)
data("individuals")
data("group_by_individual")

# Generate network
network <- get_network(gbi)

# Create a species similarity matrix
species <- array(0,dim(network))

# Create a sex similarity matrix
sex <- array(0,dim(network))

# Fill each matrix with 1 (same) or 0 (different)
for (i in 1:nrow(network)) {
	species[i,-i] <- as.numeric(inds$SPECIES[i] == inds$SPECIES[-i])
	sex[i,-i] <- as.numeric(inds$SEX[i] == inds$SEX[-i])
}

# Run mrqap.dsp
# Note randomisations are limited to 10 to reduce runtime
reg <- mrqap.dsp(network ~ species + sex, randomisations=10)

# Look at results
reg

}
