\name{summary.solveLP}
\alias{summary.solveLP}
\title{summary.solveLP}

\description{
  This function returns the results of the Linear Programming function.
}

\usage{
   summary.solveLP(object,...)
}

\arguments{
  \item{object}{an object of type \code{solveLP}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.solveLP} returns an object of type \code{solveLP}.
}

\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}
}

\seealso{
   \code{\link{solveLP}}, \code{\link{print.solveLP}},
   \code{\link{readMps}}, \code{\link{writeMps}}
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

## Maximize the gross margin
res <- solveLP( cvec, bvec, Amat, TRUE )

## print the results
summary.solveLP( res )
}

\keyword{ optimize }



