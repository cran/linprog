solveLP <- function( cvec, bvec, Amat, maximum=FALSE, maxiter=1000,
                     zero=1e-10, lpSolve=FALSE, verbose=FALSE )
{
   if( lpSolve ) {
      library( lpSolve )
      if( maximum ) {
         direction <- "max"
      } else {
         direction <- "min"
      }
      lpres <- lp (direction = direction, cvec, Amat, rep("<=",length(bvec)), bvec )

      solution  <- lpres$solution
      objval    <- lpres$objval
      Tab       <- NULL
      iter1     <- NULL
      iter2     <- NULL
      allvar    <- NULL
      basvar    <- NULL
      con       <- NULL

   } else {
      nVar <- length(cvec)  # number of variables
      nCon <- length(bvec)  # number of constraints
      iter1 <- 0
      iter2 <- 0

      if(maximum) cvec <- -cvec

      ## Labels
      if( is.null(names(cvec))) {
         clab <- as.character(1:nVar)
      } else {
         clab <- names(cvec)
      }
      if( is.null(names(bvec))) {
         blab <- as.character(1:nCon)
      } else {
         blab <- names(bvec)
      }

      ## Slack Variables
      for(i in 1:nCon) clab <- c( clab, paste("S", blab[i] ) )
      cvec <- c( cvec, rep( 0, nCon ) )

      ## Tableau ( Basic Variables, Slack,Variables, P0, Z-C )
      Tab <- rbind( cbind( Amat, diag( 1, nCon, nCon ), bvec ),
                 c( cvec, 0 ) )
      rownames(Tab) <- c( blab, "Z-C" )
      colnames(Tab) <- c( clab, "P0" )
      if(verbose) print(Tab)

      ## searching for feasible solution for starting
      # basis: Zero Solution ( Basic Variables = Slack Variables )
      basis <- c( (nVar+1) : (nVar+nCon) )
      if(min(Tab[ 1:nCon, nVar+nCon+1]) < 0 ) {
         Tab2 <- Tab
         Tab2 <- rbind( Tab2, rep(0, ncol(Tab2) ) )
         rownames(Tab2)[nCon+2] <- "M Z-C"     # additional artificial 'Z-C' row
         basis2 <- basis
         nArt   <- 0       # number of artificial variables
         for(i in 1:nCon) {
            if(Tab[ i, nVar+nCon+1] < 0 ) {
               Tab2[ i, ] <- -Tab2[ i, ]
               Tab2 <- cbind( Tab2[ , 1:(nVar+nCon+nArt) ], rep(0,nCon+2),
                              Tab2[ , (nVar+nCon+nArt+1)] )
               nArt <- nArt + 1
               colnames(Tab2)[ nVar+nCon+nArt ] <- paste("M", rownames(Tab2)[i] )
               Tab2[ i, nVar+nCon+nArt ] <- 1
               Tab2[ nCon+2, nVar+nCon+nArt ] <- 1
               # put artificial variables in basis
               rownames(Tab2)[ i ] <- paste("M", rownames(Tab2)[i] )
               basis2[i] <- nVar+nCon+nArt
            }
         }
         for(i in 1:nCon) {    # artificial objective function (Z-C)
            if(Tab[ i, nVar+nCon+1] < 0 ) {
               Tab2[nCon+2, 1:(nVar+nCon+nArt)] <- Tab2[nCon+2, 1:(nVar+nCon+nArt)] -
                                                   Tab2[ i , 1:(nVar+nCon+nArt)]
            }
         }
         for(i in 1:nCon) {    # value of artificial objective function
            Tab2[nCon+2, nVar+nCon+nArt+1 ] <- Tab2[nCon+2, nVar+nCon+nArt+1 ] -
                         Tab2[i, nVar+nCon+nArt+1] * Tab2[nCon+2, basis[i] ]
         }
         colnames(Tab2)[ nVar+nCon+nArt+1 ] <- "P0"
         if(verbose) print(Tab2)

         ## Simplex algorithm (Phase 1)
         while( min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ] ) < 0 & iter1 < maxiter) {
            iter1 <- iter1 + 1
            ## Pivot
            Tab[ abs(Tab) < zero ] <- 0
            pcolumn <- which.min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ]) # Pivot column
            rwerte  <- Tab2[ 1:nCon , nVar+nCon+nArt+1 ] / Tab2[ 1:nCon , pcolumn ] # R-values
            rwerte[ Tab2[1:nCon, pcolumn ] <= 0 ] <- max(rwerte)+1
            prow  <- which.min( rwerte )    # Pivot row
            if(verbose) {
               cat( paste( "\nPivot Column:", as.character(pcolumn),
                           "(",colnames(Tab2)[pcolumn],")\n" ) )
               cat( paste( "Pivot Row:", as.character( prow ) ,"(",rownames(Tab2)[prow],")\n\n") )
            }

            ## New Basis
            basis[prow] <- pcolumn
            rownames(Tab2)[prow] <- colnames(Tab2)[pcolumn]

            ## New Tableau
            Tab2[ prow, ] <- Tab2[ prow, ] / Tab2[ prow, pcolumn ]
            for( i in 1:(nCon+2) ) {
               if( i != prow ) {
                  Tab2[ i, ] <- Tab2[ i, ] - Tab2[ prow, ] * Tab2[ i, pcolumn ]
               }
            }
            if(verbose) print(Tab2)
         }
         if(iter1 >= maxiter ) warning("Simplex algorithm (phase 1) did reach optimum.")
         Tab <- cbind( Tab2[ 1:(nCon+1), 1:(nCon+nVar) ], Tab2[ 1:(nCon+1), nVar+nCon+nArt+1 ] )
         if(verbose) print(Tab)
      }
      ## Simplex algorithm (Phase 2)
      while( min( Tab[ nCon+1, 1:(nVar+nCon) ] ) < 0  & iter2 < maxiter) {
         iter2 <- iter2 + 1
         ## Pivot
         Tab[ abs(Tab) < zero ] <- 0
         pcolumn <- which.min( Tab[ nCon+1, 1:(nVar+nCon) ]) # Pivot column
         rwerte  <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon , pcolumn ]     # R-values
         rwerte[ Tab[1:nCon, pcolumn ] <= 0 ] <- max(rwerte)+1
         prow  <- which.min( rwerte )    # Pivot row
         if(verbose) {
            cat( paste( "\nPivot Column:", as.character(pcolumn),
                        "(",colnames(Tab)[pcolumn],")\n" ) )
            cat( paste( "Pivot Row:", as.character( prow ) ,"(",rownames(Tab)[prow],")\n\n") )
         }

         ## New Basis
         basis[prow] <- pcolumn
         rownames(Tab)[prow] <- colnames(Tab)[pcolumn]

         ## New Tableau
         Tab[ prow, ] <- Tab[ prow, ] / Tab[ prow, pcolumn ]
         for( i in 1:(nCon+1) ) {
            if( i != prow ) {
               Tab[ i, ] <- Tab[ i, ] - Tab[ prow, ] * Tab[ i, pcolumn ]
            }
         }
         if(verbose) print(Tab)
      }
      if(iter2 >= maxiter ) warning("Simplex algorithm (phase 2) did reach optimum.")
      ## Results: Basic Variables
      basvar <- matrix( NA, nCon, 1 )
      colnames(basvar) <- c("opt")
      rownames(basvar) <- rep("a",nCon)
      for( i in 1:nCon ) {
         rownames(basvar)[i] <- clab[sort(basis)[i]]
         basvar[i,1] <- Tab[ which(basis==sort(basis)[i]), nVar+nCon+1 ]
      }

      ## Results: All Variables (Including Slack Variables)
      allvar <- matrix( NA, nVar+nCon, 6 )
      colnames(allvar) <- c("opt", "c", "min c","max c","marg.", "marg.reg." )
      rownames(allvar) <- clab
      for( i in 1:(nVar+nCon) ) {
         if(i %in% basis ) {
            allvar[i,1] <- Tab[ which(basis==i), nVar+nCon+1 ]
            ## Stability of Basic Variables
            quot <- Tab[ nCon+1, 1:(nVar+nCon) ] / Tab[ which(basis==i), 1:(nVar+nCon) ]
            if(maximum) {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  op <- options()
                  options(warn=-1)
                  allvar[i,3] <- -cvec[ i ] - min(quot[quot>0 & !is.na(quot)])
                  options(op)
               }
               if(min(quot[!is.na(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ) {
                     allvar[i,4] <- -cvec[ i ] - max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar[i,4] <- Inf
                  }
               } else {
                  allvar[i,4] <- Inf
               }
            } else {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  allvar[i,4] <- cvec[ i ] + min(quot[quot>0 & !is.na(quot)])
               }
               if(min(quot[!is.na(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ){
                     allvar[i,3] <- cvec[ i ] + max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar[i,3] <- -Inf
                  }
               } else {
                  allvar[i,3] <- -Inf
               }
            }
         } else {
            allvar[i,1] <- 0
            if( i <= nVar ) {
               if(maximum) {
                  allvar[i,3] <- -Inf
                  allvar[i,4] <- Tab[ nCon+1, i ] - cvec[i]
               } else {
                  allvar[i,3] <- -Tab[ nCon+1, i ] + cvec[i]
                  allvar[i,4] <- Inf
               }
            }
         }
         if(maximum) {
            allvar[i,2] <- -cvec[i]
         } else {
            allvar[i,2] <- cvec[i]
         }
         # marginal contribution to objective function (Shadow prices)
         if( !( ( i %in% basis ) & ( i <= nVar ) ) ) {
            allvar[ i, 5] <- Tab[ nCon+1, i ] * (-1)^maximum
            quot <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon, i ]
            op <- options()
            options(warn=-1)
            if(min(quot[quot>0 & !is.na(quot)]) != Inf ) {
               allvar[i,6] <- min(quot[quot>0 & !is.na(quot)])
            } else {
               allvar[i,6] <- NA
            }
            options(op)
         }
      }
      allvar[ (allvar[,3] >  1e16), 3 ] <-  Inf
      allvar[ (allvar[,3] < -1e16), 3 ] <- -Inf

      ## Results: Constraints
      con <- matrix( NA, nCon, 5 )
      colnames(con) <- c("max", "actual", "diff", "dual price","dual.reg")
      rownames(con) <- blab
      con[ , 1 ] <- bvec
      for(i in 1: nCon) {
         if( (i+nVar) %in% basis ) {
            con[ i, 2 ] <- bvec[i] - Tab[ which((i+nVar)==basis), nVar+nCon+1 ]
         } else {
            con[ i, 2 ] <- bvec[i]
         }
         con[ i, 4 ] <- -allvar[ i+nVar , 5 ]
         con[ i, 5 ] <-  allvar[ i+nVar , 6 ]
      }
      con[ , 3 ] <- con[ , 1 ] - con[ , 2 ]

      objval   <- -Tab[ nCon+1, nCon+nVar+1 ] * (-1)^maximum
      basvar   <- round( basvar, digits=10 )
      con      <- round( con, digits=10 )
      allvar   <- round( allvar, digits=10 )
      solution <- allvar[ 1 : nVar, 1 ]
   }

   ## List of Results
   result <- list()
   result$opt      <- round( objval, digits=10 )
   result$iter1    <- ifelse( iter1 < maxiter, iter1, -iter1 )
   result$iter2    <- ifelse( iter2 < maxiter, iter2, -iter2 )
   result$solution <- solution
   result$basvar   <- basvar
   result$con      <- con
   result$allvar   <- allvar
   result$maximum  <- maximum
   result$Tab      <- Tab
   result$lpSolve  <- lpSolve
   class(result)   <- "solveLP"
   result
}

## print the (summary) results
summary.solveLP <- function(object,...) {
  summary.solveLP <- object
  summary.solveLP
}

## print the results
print.solveLP <- function( x, digits=6,... ) {
   object <- x

   save.digits <- unlist(options(digits=digits))
   on.exit(options(digits=save.digits))

   cat("\n\nResults of Linear Programming / Linear Optimization\n")
   cat("\nObjective function")
   if( object$maximum ) {
      cat(" (Maximum):\n")
   } else {
      cat(" (Minimum):\n")
   }
   print( object$opt )
   cat("\nIterations in phase 1: ")
   cat(abs(object$iter1))
   if( !object$lpSolve ) {
      if( object$iter1 < 0 ) cat(" (equals 'maxiter' !!!)")
   }
   cat("\nIterations in phase 2: ")
   cat(abs(object$iter2))
   if( !object$lpSolve ) {
      if( object$iter2 < 0 ) cat(" (equals 'maxiter' !!!)")
   }
   cat("\n\nBasic Variables\n")
   print( object$basvar )
   cat("\nConstraints\n")
   print( object$con )
   cat("\nAll Variables (including slack variables)\n")
   print( object$allvar )
   cat("\n")
}

readMps <- function( file, solve=FALSE, maximum=FALSE ) {

   mps <- readLines(file)
   i <- 1

   ## Name
   while( substr( mps[i], 1, 4 ) != "NAME" & i < length(mps) ) {
      i <- i + 1
   }
   if( substr( mps[i], 1, 4 ) == "NAME" ) {
      name <- substr( mps[i], 15, nchar( mps[i] ) )
   } else {
      stop( "MPS file must have a line starting with 'NAME'" )
   }

   ## Rows / Constraints
   while( substr( mps[i], 1, 4 ) != "ROWS" & i < length(mps) ) {
      i <- i + 1
   }
   if( substr( mps[i], 1, 4 ) != "ROWS" ) stop( "MPS file must have a line starting with 'ROWS'" )

   objname <- NULL
   bvec <- NULL   # constraints
   svec <- NULL   # sign of the constraints

   i <- i + 1
   while( substr( mps[i], 1, 7 ) != "COLUMNS" & i < length(mps)) {
      sign <- substr( mps[i], 2,2 )
      if( sign == "E" ) stop("Equaltity constraints are not implemented yet.")
      if( !( sign %in% c("E", "L", "G", "N" ) ) ) {
         sign <- substr( mps[i], 3,3 )
      }
      if( sign %in% c("E", "L", "G" ) ) {
         rname <- strsplit( mps[i], " " )[[1]][length( strsplit( mps[i], " " )[[1]] )]
         bvec <- c( bvec, 0 )
         svec <- c( svec, sign )
         names(bvec)[length( bvec ) ]  <-  rname
         names(svec)[length( svec ) ]  <-  rname
      } else {
         if( sign == "N" ) {
            if( is.null( objname ) ) {
               objname <- strsplit( mps[i], " " )[[1]][length( strsplit( mps[i], " " )[[1]] )]
            }
         } else {
            stop("the 2nd or 3rd column of the rows section must be 'N', 'E', 'L' or 'G'")
         }
      }
      i <- i + 1
   }

   ## Columns
   if( substr( mps[i], 1, 7 ) != "COLUMNS" )
      stop( "MPS file must have a line starting with 'COLUMS'" )
   cvec <- NULL
   Amat <- matrix(0, length(bvec), 0 )
   rownames(Amat) <- names(bvec)
   i <- i + 1
   while( substr( mps[i], 1, 3 ) != "RHS" & i < length(mps)) {
      temp <- strsplit( mps[i], " " )[[1]]
      temp <- temp[ temp != "" ]
      if( !(temp[1] %in% colnames(Amat) ) ) {
         cvec <- c( cvec, 0 )
         Amat <- cbind( Amat, rep( 0, nrow(Amat) ) )
         names(cvec)[length(cvec)]  <- temp[1]
         colnames(Amat)[ncol(Amat)] <- temp[1]
      }
      for( j in 1:((length(temp)-1)/2) ) {
         if( temp[ 2*j ] == objname ) {
            cvec[ temp[ 1 ] ] <- as.numeric( temp[ 2*j + 1 ] )
         } else {
            if( temp[ 2*j ] %in% names(bvec) ) {
               Amat[ temp[ 2*j ], temp[ 1 ] ] <- as.numeric( temp[ 2*j + 1 ] )
            } else {
               stop( paste( "Constraint name '",temp[ 2*j ],"' is not defined", sep="") )
            }
         }
      }
      i <- i + 1
   }

   ## Restriction values
   if( substr( mps[i], 1, 3 ) != "RHS" ) stop( "MPS file must have a line starting with 'RHS'" )
   i <- i + 1
   while( substr( mps[i], 1, 6 ) != "BOUNDS" & substr( mps[i], 1, 6 ) != "ENDATA" & i < length(mps)) {
      temp <- strsplit( mps[i], " " )[[1]]
      temp <- temp[ temp != "" ]
      for( j in 1:((length(temp)-1)/2) ) {
         if( temp[ 2*j ] %in% names(bvec) ) {
            bvec[ temp[ 2*j ] ] <- as.numeric( temp[ 2*j + 1 ] )
         } else {
            stop( paste( "Constraint name '",temp[ 2*j ],"' is not defined", sep="") )
         }
      }
      i <- i + 1
   }

   ## Bounds
   if( substr( mps[i], 1, 6 ) == "BOUNDS" ) {
      i <- i + 1
      while( substr( mps[i], 1, 6 ) != "ENDATA" & i <= length(mps)) {
         temp <- strsplit( mps[i], " " )[[1]]
         temp <- temp[ temp != "" ]
         if( temp[ 3 ] %in% colnames(Amat) ) {
            if(temp[1] == "UP") {
               svec <- c( svec, "L" )
               bvec <- c( bvec, as.numeric(temp[ 4 ]) )
               Amat <- rbind( Amat, rep( 0, ncol(Amat) ) )
               Amat[ nrow(Amat), temp[3] ] <- 1
               names( svec )[length(svec)] <- paste(temp[1], temp[3], sep="" )
               names( bvec )[length(bvec)] <- paste(temp[1], temp[3], sep="" )
               rownames( Amat )[nrow(Amat)] <- paste(temp[1], temp[3], sep="" )
            } else {
               if( temp[1] %in% c("LO","FX","FR") ) {
                  stop("'LO', 'FX', and 'FR' Bounds are not implemented yet")
               } else {
                  stop(" A 'BOUND' line must start with 'UP', 'LO', 'FX' or 'FR'")
               }
            }
         } else {
            stop( paste( "Variable name '",temp[ 3 ],"' is not defined", sep="") )
         }
         i <- i + 1
      }
   }

   if( substr( mps[i], 1, 6 ) != "ENDATA" )
      stop( "MPS file must have a line starting with 'ENDDATA'" )

   ## Changing 'Greater' constraints to 'Lower' constraints
   for( j in 1:length(svec) ) {
      if(svec[j] == "G" ) {
         bvec[ j ]  <- - bvec[ j ]
         Amat[ j, ] <- - Amat[ j, ]
      }
   }

   res <- NULL
   if(solve) {
      res <- solveLP(cvec,bvec,Amat,maximum)
   }
   return( name=name, cvec=cvec, bvec=bvec, Amat=Amat, res=res )
}

writeMps <- function( file, cvec, bvec, Amat, name="LP" ) {
   nCon <- length(bvec)
   nVar <- length(cvec)

   if( is.null( names(bvec) ) ) {
      blab <- rep("",nCon)
      for(i in 1:nCon) {
         blab[i] <- paste("R_",as.character(i))
      }
   } else {
      blab <- names(bvec)
      for(i in 1:nCon) {
         blab[i] <- gsub(" ","",blab[i])
         if( nchar( blab[i] ) > 8 ) {
            blab[i] <- substr( blab[i], 1, 8 )
         }
         j <- 2
         while( i>1 & blab[i] %in% blab[1:(i-1)])  {
            blab[i] <- paste( substr(blab[i], 1, 7-nchar(as.character(j))),
                              "_", as.character(j), sep="" )
            j <- j+1
         }
      }
   }

   if( is.null( names(cvec) ) ) {
      clab <- rep("",nVar)
      for(i in 1:nVar) {
         clab[i] <- paste("C_",as.character(i))
      }
   } else {
      clab <- names(cvec)
      for(i in 1:nVar) {
         clab[i] <- gsub(" ","",clab[i])
         if( nchar( clab[i] ) > 8 ) {
            clab[i] <- substr( clab[i], 1, 8 )
         }
         j <- 2
         while( i>1 & clab[i] %in% clab[1:(i-1)])  {
            clab[i] <- paste( substr(clab[i], 1, 7-nchar(as.character(j))),
                              "_", as.character(j), sep="" )
            j <- j+1
         }
      }
   }

   write( paste("NAME          ",name,sep=""), file )

   write( "ROWS", file, append=TRUE )
   write( " N  obj", file, append=TRUE )
   for(i in 1:nCon) {
      write( paste(" L  ",blab[i],sep="" ), file, append=TRUE )
   }

   write( "COLUMNS", file, append=TRUE )
   for(i in 1:nVar) {
      line <- paste("    ",clab[i], sep="" )
      line <- paste( line, paste( rep( " ", 14-nchar(line)), collapse=""), "obj", sep="")
      line <- paste( line, paste( rep( " ", 36-nchar(line) - nchar(
         as.character( signif( cvec[i], 10 )))), collapse=""),
         as.character( signif( cvec[i], 10)), sep="")
      write( line, file, append=TRUE )
      for(j in 1:nCon) {
         if( Amat[j,i] != 0 ) {
            line <- paste("    ",clab[i], sep="" )
            line <- paste( line, paste( rep( " ", 14-nchar(line)),collapse=""), blab[j], sep="")
            line <- paste( line, paste( rep( " ", 36 - nchar(line) - nchar(
                      as.character( signif( Amat[j,i], 10 )))), collapse=""),
                      as.character( signif( Amat[j,i], 10 )), sep="")
            write( line, file, append=TRUE )
         }
      }
   }

   write( "RHS", file, append=TRUE )
   for(i in 1:nCon) {
      line <- paste("    RHS       ",blab[i], sep="" )
      line <- paste( line, paste( rep( " ", 36-nchar(line) - nchar(
               as.character( signif( bvec[i], 10 )))), collapse=""),
               as.character( signif( bvec[i], 10 )), sep="")
      write( line, file, append=TRUE )
   }

   write( "ENDATA", file, append=TRUE )
}
