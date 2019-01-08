\name{mod_regression}
\alias{mod_regression}
\title{moderated regression}
\description{moderated regression}
\usage{mod_regression(idv, dv, mod1, mod2 = NULL, multMOD1 = 1, multMOD2 = 1, ci = 95)}
\arguments{
  \item{idv}{
  \code{}A numeric vector representing participant scores on the independent/predictor variable.}
  \item{dv}{
  \code{}A numeric vector representing the dependent/outcome variable.}
  \item{mod1}{
  \code{}Either a numeric vector (for continuous moderators) or factor (for categorical moderators) representing the first moderator.}
  \item{mod2}{
  \code{}An optional numeric vector or factor representing the second moderator.}
  \item{multMOD1}{
  \code{}The number of SD above and below the mean of the first moderator variable, for which simple slope statistics will be calculated.}
  \item{multMOD2}{
  \code{}The number of SD above and below the mean of the second moderator variable, for which simple slope statistics will be calculated.}
  \item{ci}{
  \code{}The confidence level to be calculated. Default is 95.}
}
\value{ A list with the following elements:
 \item{lm}{The significance test for the multiple regression equation.}
 \item{beta}{The raw and standardized beta weights for the full equation.}
 \item{intercept}{The intercept.}
 \item{slopes}{The raw (unstandardized) simple slope coefficients for the DV on the IDV at three levels of the moderator.}
 \item{std_slopes}{The standardized simple slopes (beta, SE, and 95\% confidence intervals).}
 \item{s_slope}{The score on the moderator for which the simple slope for the DV on the IDV is zero (flat).}
 \item{s_regression}{The score on the independent variable at which the simple regression lines for mod = high and mod = low intersect.}
}
\author{Zakary A. Draper & Brian P. O'Connor }
\examples{
data(mtcars)
attach(mtcars)
mod_regression(hp, mpg, qsec)
}
