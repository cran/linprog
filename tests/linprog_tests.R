library( linprog )

## Example 1
## Steinhauser, Langbehn and Peters (1992)
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Cows","Bulls","Pigs")
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )
result1a <- solveLP( cvec, bvec, Amat, TRUE )
print( result1a )

## Example 2
## example 1.1.3 of Witte, Deppe and Born (1975)
cvec <- c(2.5, 2 )  # prices of feed
names(cvec) <- c("Feed1","Feed2")
bvec <- c( -10, -1.5, 12)
names(bvec) <- c("Protein","Fat","Fibre")
Amat <- rbind( c(-1.6,-2.4 ),
               c(-0.5,-0.2 ),
               c( 2.0, 2.0 ) )
result2a <- solveLP( cvec, bvec, Amat )
print( result2a )
