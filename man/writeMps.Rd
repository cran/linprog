\name{writeMps}
\alias{writeMps}

\title{write MPS files}

\description{
  This function writes MPS files - the standard format for Linear Programming problems.
}

\usage{
writeMps( file, cvec, bvec, Amat, name="LP" )
}

\arguments{
  \item{file}{a character string naming the file to write.}
  \item{cvec}{vector \eqn{c}.}
  \item{bvec}{vector \eqn{b}.}
  \item{Amat}{matrix \eqn{A}.}
  \item{name}{an optional name for the Linear Programming problem.}
}

\details{
   The exported LP can be solved by running other software on this MPS file
   (e.g. \code{lp_solve}, see \url{ftp://ftp.es.ele.tue.nl/pub/lp_solve}).
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{solveLP}}, \code{\link{readMps}}
}


\examples{

## example of Steinhauser, Langbehn and Peters (1992)
\dontrun{library( linprog )}

## Production activities
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Milk","Bulls","Pigs")

## Constraints (quasi-fix factors)
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")

## Needs of Production activities
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )

## Write to MPS file
writeMps( "steinh.mps", cvec, bvec, Amat, "Steinhauser" )
}

\keyword{ optimize }



