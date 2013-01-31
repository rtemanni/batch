\name{interpretZ}
\alias{interpretZ}
\title{Compute transformation of non-identifiable parameters}
\usage{
  interpretZ(X, Z)
}
\arguments{
  \item{X}{full non-identifiable design matrix}

  \item{Z}{reparameterized design matrix}
}
\value{
  A matrix indacting the identifiable linear combinations
  of the X cols that make up the Z cols
}
\description{
  Compute transformation of non-identifiable parameters
}

