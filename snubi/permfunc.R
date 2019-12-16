#######################################
# permfunc.R
#######################################
# Functions for doing permutation tests
#######################################
# Utility function
#     returns binary representation of 1:(2^n)
binary.v <-
  function(n)
  {
    x <- 1:(2^n)
    mx <- max(x)
    digits <- floor(log2(mx))
    ans <- 0:(digits-1); lx <- length(x)
    x <- matrix(rep(x,rep(digits, lx)),ncol=lx)
    x <- (x %/% 2^ans) %% 2
  }

# Function to perform a paired permutation test
#     Input: differences, d
#            no. permutations
# Use paired.perm.test(d, n.perm=NULL) to 
#     do the exact test
paired.perm.test <-
  function(d, n.perm=1000, pval=FALSE)
  {
    n <- length(d)
    obs <- t.test(d)$statistic
    
    if(is.null(n.perm)) { # do exact test
      ind <- binary.v(n)
      allt <- apply(ind,2,function(x,y)
        t.test((2*x-1)*y)$statistic,d)
    }
    else { # do n.perm samples
      allt <- 1:n.perm
      for(i in 1:n.perm) 
        allt[i] <- t.test(d*sample(c(-1,1),n,repl=TRUE))$statistic
    }
    
    if(pval) return(mean(abs(allt) >= abs(obs)))
    allt
  }


# Function to perform permutation test
#     x, y = the two samples
#     n.perm = number of permutations
#     var.equal = passed to the function t.test
# Use perm.test(x, y, n.perm=NULL)
#     to get exact P-value
perm.test <-
  function(x, y, n.perm=1000, var.equal=TRUE, pval=FALSE)
  {
    # number of data points
    kx <- length(x)
    ky <- length(y)
    n <- kx + ky
    
    # observed statistic
    obs <- t.test(x,y, var.equal=var.equal)$statistic
    
    # Data re-compiled
    X <- c(x,y)
    z <- rep(1:0,c(kx,ky))
    
    if(is.null(n.perm)) { # do exact permutation test
      o <- binary.v(n)  # indicator of all possible samples
      o <- o[,apply(o,2,sum)==kx]  
      nc <- choose(n,kx)
      allt <- 1:nc
      for(i in 1:nc) {
        xn <- X[o[,i]==1]
        yn <- X[o[,i]==0]
        allt[i] <- t.test(xn,yn,var.equal=var.equal)$statistic
      }
    }
    else { # do 1000 permutations of the data
      allt <- 1:n.perm
      for(i in 1:n.perm) {
        z <- sample(z)
        xn <- X[z==1]
        yn <- X[z==0]
        allt[i] <- t.test(xn,yn,var.equal=var.equal)$statistic
      }
    }
    
    if(pval) return(mean(abs(allt) >= abs(obs)))
    allt
  }

# end of permfunc.R