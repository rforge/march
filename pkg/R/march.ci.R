# This file is part of March.
# It contains functions for the computation of confidence intervals.
# 
# March is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyrigth 2014-2020, Ogier Maitre, Kevin Emery, André Berchtold 
# andre.berchtold@unil.ch 


###############################################################################
# Thompson allows to compute confidence intervals according to a given model,
# using thompson's confidence interval method describe into: Thompson, S.K. (1987) 
# "Sample size for estimating multinomial proportions," American Statistician, 41, 42-46.
# Adaptation to markov models is described into : Berchtold, "Confidence Intervals for Markovian Models"
###############################################################################

#' Thompson Confidence Intervals for a march.Model.
#' 
#' Compute the confidence intervals using Thompson's formula on a march.Model
#' object. See Thompson SK (1987) Sample size for estimating multinomial proportions,
#' American Statistician 41:42-46, for details.
#' 
#' @param object the march.Model object on which compute the confidence intervals.
#' @param alpha the significance level among : 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.025, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001.
#' 
#' @return A list of half-length confidence intervals for each probability distribution of the considered model.
#' @author Ogier Maitre, Kevin Emery
#' @example tests/examples/march.thompson.example.R
#' @export 
march.thompson <- function(object,alpha){}

march.model.thompson <- function(object,alpha){
  warning("Confidence interval with thompson formula cannot be computed for abstract class \"model\", check the parameters of the call to march.thompson")
}

march.indep.thompson <- function(object,alpha){
  d2n <- march.ci.h.d2n(alpha)
  
  d <- sqrt(d2n/object@dsL)
  
  # Display
  cat("\nHalf CI for the independence model :\n")
  cat("------------------------------------\n")
  print(d)
  cat("\n")
}

march.mc.thompson <- function(object,alpha){
  d2n <- march.ci.h.d2n(alpha)
  d <- array(NA,object@y@K^object@order)
  
  for( i in 1:length(d)){
    s <- sum(object@RT[i,])
    if( s>0 ){
      d[i]<-sqrt(d2n/s)
    }
  }
  
  # Display
  cat("\nHalf CI for each row of the transition matrix :\n")
  cat("-----------------------------------------------\n")
  print(d)
  cat("\n") 

  }


march.mtd.thompson <- function(object, alpha){
  # TO DO: CI for the high-order matrix
  
  d2n <- march.ci.h.d2n(alpha)
  
  dphi <- d2n/object@dsL
  
  if(dim(object@Q)[1]>1){
    is_mtdg <- TRUE
  }else{
    is_mtdg <- FALSE
  }
  l <- march.mtd.h.n(object,object@y,is_mtdg)
  
  dQ <- list()
  if(is_mtdg==FALSE){
    dQ[[1]] <- d2n/rowSums(l$nki_0)
  }else{
    for(ord in 1:object@order){
      dQ[[ord]] <- d2n/rowSums(l$nki_0[ord,,])
    }
  }
  
  # Covariates
  dS <- list()
  if(sum(object@MCovar)>0){
    placeCovar <- which(object@MCovar==1)
    for(i in 1:sum(object@MCovar)){
      dS[[i]] <- d2n/rowSums(l$numcov[i,1:object@y@Kcov[placeCovar[i]],])
    }
  }
  
  # High-order transition matrix
  CQ <- l$nki_0
  
  if(is_mtdg==FALSE){
    # MTD model

    # Matrix of index
    INDEX <- BuildArrayCombinations(ncol(CQ),(object@order-1),0,0)
    
    # Matrix of number of data
    NHO <- matrix(0,nrow=(ncol(CQ)^object@order),ncol=ncol(CQ))
    
    for (i in 1:(ncol(CQ)^object@order)){
      for (j in 1:ncol(CQ)){
        for (k in 1:object@order){
          NHO[i,j] <- NHO[i,j]+object@phi[k]*CQ[INDEX[i,k],j]
        } 
      }
    }    
    
  }else{
    # MTDg model
    
    # Matrix of index
    INDEX <- BuildArrayCombinations(ncol(CQ[1,,]),(object@order-1),0,0)
    
    # Matrix of number of data
    NHO <- matrix(0,nrow=(ncol(CQ[1,,])^object@order),ncol=ncol(CQ))
    
    for (i in 1:(ncol(CQ[1,,])^object@order)){
      for (j in 1:ncol(CQ[1,,])){
        for (k in 1:object@order){
          NHO[i,j] <- NHO[i,j]+object@phi[k]*CQ[ord,INDEX[i,k],j]
        } 
      }
    }  
  }

  # CI for the high-order transition matrix
  dNHO <- d2n/rowSums(NHO)

  
  # Display
  cat("\nHalf CI for the vector of weights :\n")
  cat("-----------------------------------\n")
  print(dphi)
  cat("\n")
  
  if (dim(object@Q)[1]==1){
    # MTD
    cat("Half CI for each row of the transition matrix :\n")
    cat("---------------------------------------------\n")
    print(dQ[[1]])
    cat("\n")
  } else {
    # MTDg
    for (ord in 1:dim(object@Q)[1]){
      cat("Half CI for each row of the transition matrix of lag",ord,":\n")
      cat("  ------------------------------------------------------\n")
      print(dQ[[ord]])
      cat("\n")
    }
  }
 
  if(sum(object@MCovar)>0){
    # Covariates
    for (cov in 1:sum(object@MCovar)){
      cat("Half CI for each row of the transition matrix of covariate",cov,":\n")
      cat("  ---------------------------------------------------------\n")
      print(dS[[cov]])
      cat("\n")
    }    
  }
    
    
  cat("Half CI for each row of the high-order transition matrix :\n")
  cat("----------------------------------------------------------\n")
  print(dNHO)
  cat("\n")
  
}


march.dcmm.thompson <- function(object,alpha){
  
  stop("This part is not yet implemented")
  # ys <- march.dataset.h.extractSequence(object@y,1)
  # alpha<-march.dcmm.forward(object,ys)
  # beta<-march.dcmm.backward(object,ys)
  # 
  # C
} 


#This part create the generic method and describe how a call to this generic
#has to be redirected to the rigth method, according to the considered object.
setGeneric(name="march.thompson",def=function(object,alpha)march.model.thompson(object,alpha))
setMethod(f="march.thompson",signature=signature("march.Indep",alpha="numeric"),definition=march.indep.thompson)
setMethod(f="march.thompson",signature=signature("march.Mc",alpha="numeric"),definition=march.mc.thompson)
setMethod(f="march.thompson",signature=signature("march.Mtd",alpha="numeric"),definition=march.mtd.thompson)
setMethod(f="march.thompson",signature=signature("march.Dcmm",alpha="numeric"),definition=march.dcmm.thompson)


###############################################################################
# Bailey allows to compute confidence intervals according to a given model,
# using Bailey's confidence interval method describe in: Bailey, B.J.R. (1980) 
# "Large sample simultaneous confidence intervals for the multinomial probabilities based on
# transformation of the cell frequencies." Technometrics 1980, 22, 583–589.
# Adaptation to Markov models is described into : Berchtold, "Confidence Intervals for Markovian Models"
###############################################################################

#' Bailey Confidence Intervals for a march.Model.
#' 
#' Compute the confidence intervals using Bailey's formula on a march.Model
#' object. See Bailey BJR (1980) Large sample simultaneous confidence intervals
#' for the multinomial probabilities based ontransformation of the cell frequencies,
#' Technometrics 22:583–589, for details.
#' 
#' @param object the march.Model object on which compute the confidence intervals.
#' @param alpha the significance level.
#' 
#' @return A list of half-length confidence intervals for each probability distribution of the considered model.
#' @author Berchtold André
#' @example tests/examples/march.bailey.example.R
#' @export 
march.bailey <- function(object,alpha){}

march.model.bailey <- function(object,alpha){
  warning("Confidence interval with Bailey formula cannot be computed for abstract class \"model\", check the parameters of the call to march.bailey")
}


march.indep.bailey <- function(object,alpha){
  n <- sum(object@dsL)
  
  p <- array(NA,c(2,object@y@K))
  colnames(p) <- object@y@dictionary
  rownames(p) <- c("p-","p+")
  for( i in 1:object@y@K){
    ni <- object@indC[i]
    A <- march.ci.h.A(n,ni)
    B <- march.ci.h.B(n,ni)	
    C <- march.ci.h.C(alpha,object@y@K,n)
    
    p["p-",i] = (sqrt(A)-sqrt(C*(C+1-A)))^2/(C+1)^2
    p["p+",i] = (sqrt(B)+sqrt(C*(C+1-B)))^2/(C+1)^2
  }
  
  # Display
  cat("\nCI for the independence model :\n")
  cat("-------------------------------\n")
  cat("Lower bound :")
  prmatrix(t(p["p-",]),collab=rep("",object@y@K))
  cat("\nUpper bound :")
  prmatrix(t(p["p+",]),collab=rep("",object@y@K))
  cat("\n")

}


march.mc.bailey <- function(object,alpha){
  
  NHO <- object@RT
  rNHO <- rowSums(NHO)  
  
  NHO.l <- matrix(NA,nrow=nrow(NHO),ncol=ncol(NHO))
  NHO.u <- matrix(NA,nrow=nrow(NHO),ncol=ncol(NHO))
  for (i in 1:nrow(NHO)){
    for (j in 1:ncol(NHO)){
      res <- march.bailey.ci(NHO[i,j],rNHO[i],alpha,ncol(NHO))
      NHO.l[i,j] <- res[[1]]
      NHO.u[i,j] <- res[[2]]
    }
  }
  
  # Display
  cat("\nCI for the transition matrix :\n")
  cat("------------------------------\n")
  cat("Lower bound :")
  prmatrix(NHO.l,collab=rep("",object@y@K))
  cat("\nUpper bound :")
  prmatrix(NHO.u,collab=rep("",object@y@K))
  cat("\n") 
}


march.mtd.bailey <- function(object,alpha){
  
  # lag weights
  N <- object@dsL
  ni <- N*object@phi
  
  phi.l <- matrix(NA,nrow=1,ncol=object@order)
  phi.u <- matrix(NA,nrow=1,ncol=object@order)
  for (i in 1:object@order){
    res <- march.bailey.ci(ni[i],N,alpha,object@order)
    phi.l[i] <- res[[1]]
    phi.u[i] <- res[[2]]
  }
  
  # transition probabilities

  if (dim(object@Q)[1]==1){
    
    # A. MTD model
    CQ <- march.mtd.h.n(object,object@y,is_mtdg=F)
    CQ <- CQ$`nki_0`
    rCQ <- rowSums(CQ)
    
    Q.l <- matrix(NA,nrow=nrow(CQ),ncol=ncol(CQ))
    Q.u <- matrix(NA,nrow=nrow(CQ),ncol=ncol(CQ))
    for (i in 1:nrow(CQ)){
      for (j in 1:ncol(CQ)){
        res <- march.bailey.ci(CQ[i,j],rCQ[i],alpha,ncol(CQ))
        Q.l[i,j] <- res[[1]]
        Q.u[i,j] <- res[[2]]
      }
    }
    
    # High-order matrix
    
    # Matrix of index
    INDEX <- BuildArrayCombinations(ncol(CQ),(object@order-1),0,0)
    
    # Matrix of number of data
    NHO <- matrix(0,nrow=(ncol(CQ)^object@order),ncol=ncol(CQ))
    
    for (i in 1:(ncol(CQ)^object@order)){
      for (j in 1:ncol(CQ)){
        for (k in 1:object@order){
          NHO[i,j] <- NHO[i,j]+object@phi[k]*CQ[INDEX[i,k],j]
        } 
      }
    }    
  } else {
    
    # B. MTDg model
    Q.l <- array(NA,dim=c(object@order,object@y@K,object@y@K))
    Q.u <- array(NA,dim=c(object@order,object@y@K,object@y@K))
    
    CQ <- march.mtd.h.n(object,object@y,is_mtdg=T)
    CQ <- CQ$`nki_0`
    
    for (ord in 1:object@order){
      QCQ <- CQ[ord,,]
      rQCQ <- rowSums(QCQ)

      for (i in 1:nrow(QCQ)){
        for (j in 1:ncol(QCQ)){
          res <- march.bailey.ci(QCQ[i,j],rQCQ[i],alpha,ncol(QCQ))
          Q.l[ord,i,j] <- res[[1]]
          Q.u[ord,i,j] <- res[[2]]
        }
      }
    }
    
    # High-order matrix
    
    # Matrix of index
    INDEX <- BuildArrayCombinations(ncol(QCQ),(object@order-1),0,0)
    
    # Matrix of number of data
    NHO <- matrix(0,nrow=(ncol(QCQ)^object@order),ncol=ncol(QCQ))
    
    for (i in 1:(ncol(QCQ)^object@order)){
      for (j in 1:ncol(QCQ)){
        for (k in 1:object@order){
          NHO[i,j] <- NHO[i,j]+object@phi[k]*CQ[k,INDEX[i,k],j]
        } 
      }
    }    
  }

  # CI for the high-order transition matrix
  rNHO <- rowSums(NHO)  
  
  NHO.l <- matrix(NA,nrow=nrow(NHO),ncol=ncol(NHO))
  NHO.u <- matrix(NA,nrow=nrow(NHO),ncol=ncol(NHO))
  for (i in 1:nrow(NHO)){
    for (j in 1:ncol(NHO)){
      res <- march.bailey.ci(NHO[i,j],rNHO[i],alpha,ncol(NHO))
      NHO.l[i,j] <- res[[1]]
      NHO.u[i,j] <- res[[2]]
    }
  }
  
  # Display
  cat("\nCI for the vector of weights :\n")
  cat("------------------------------\n")
  cat("Lower bound :\n")
  print(phi.l[1,])
  cat("\nUpper bound :\n")
  print(phi.u[1,])
  cat("\n")
  
  if (dim(object@Q)[1]==1){
    # MTD
    cat("CI for the transition matrix :\n")
    cat("------------------------------\n")
    cat("Lower bound :")
    prmatrix(Q.l,collab=rep("",object@y@K))
    cat("\nUpper bound :")
    prmatrix(Q.u,collab=rep("",object@y@K))
    cat("\n")
  } else {
    # MTDg
    for (ord in 1:dim(object@Q)[1]){
      cat("CI for the transition matrix of lag",ord,":\n")
      cat("---------------------------------------\n")
      cat("Lower bound :")
      prmatrix(Q.l[ord,,],collab=rep("",object@y@K))
      cat("\nUpper bound :")
      prmatrix(Q.u[ord,,],collab=rep("",object@y@K))
      cat("\n")
    }
  }
  
  cat("CI for the high-order transition matrix :\n")
  cat("-----------------------------------------\n")
  cat("Lower bound :")
  prmatrix(NHO.l,collab=rep("",object@y@K))
  cat("\nUpper bound :")
  prmatrix(NHO.u,collab=rep("",object@y@K))
  cat("\n")
  
  #list(phi.l,phi.u,Q.l,Q.u,NHO.l,NHO.u)
}


march.dcmm.bailey <- function(object,alpha){
  
  stop("This part is not yet implemented")
  # ys <- march.dataset.h.extractSequence(object@y,1)
  # alpha<-march.dcmm.forward(object,ys)
  # beta<-march.dcmm.backward(object,ys)
  # 
  # C
} 


#This part create the generic method and describe how a call to this generic
#has to be redirected to the rigth method, according to the considered object.
setGeneric(name="march.bailey",def=function(object,alpha)march.model.bailey(object,alpha))
setMethod(f="march.bailey",signature=signature("march.Indep",alpha="numeric"),definition=march.indep.bailey)
setMethod(f="march.bailey",signature=signature("march.Mc",alpha="numeric"),definition=march.mc.bailey)
setMethod(f="march.bailey",signature=signature("march.Mtd",alpha="numeric"),definition=march.mtd.bailey)
setMethod(f="march.bailey",signature=signature("march.Dcmm",alpha="numeric"),definition=march.dcmm.bailey)



##############################################################################
##############################################################################
# Tool functions used fo CIs

march.bailey.ci <- function(ni,N,alpha,k){
  
  if (ni>0){
    pMinus <- max((ni-(1/8))/(N+(1/8)),0)
    pPlus <- min((ni+(7/8))/(N+(1/8)),1)
    
    B <- qchisq(1-(alpha/k),1)
    C <- B/(4*N)
    
    t1 <- (sqrt(pMinus)-sqrt(C*(C+1-pMinus)))^2
    t2 <- (sqrt(pPlus)+sqrt(C*(C+1-pPlus)))^2
    t3 <- (C+1)^2
    
    LowerBd <- t1/t3
    UpperBd <- t2/t3
    
    list(LowerBd,UpperBd) 
  } else {
    list(0,0)
  }
  
}


march.ci.h.A <- function(n,ni){
  (ni-1/8)/(n+1/8)
}

march.ci.h.B <- function(n,ni){
  (ni+7/8)/(n+1/8)
}

march.ci.h.C <- function(alpha,K,n){
  qchisq(1-alpha/K,df=1)/(4*n)
}

# Get the d^2*n value as defined in Thompson 1987 p43 paper 
# ("Sample Size for Estimating Multinomial Proportions").
# The value is available for some value of alpha (as defined in the vector a).
# If the alpha value is not defined, the default value (for alpha=0.05) is
# returned.
#
# Parameters : 
#   alpha : the significance level (alpha)
#
# Returns :
#  the d^2*n associated with the alpha value, if defined, or the one
#	associated with alpha=0.05.
#
march.ci.h.d2n <- function(alpha){
  a <- c( 0.5,0.4,0.3,0.2,0.1,0.05,0.025,0.02,0.01,0.005,0.001,0.0005,0.0001 )
  d2n <- c( 0.44129,0.50729,0.60123,0.74739,1.00635,1.27359,1.55963,1.65872,1.96986,2.28514,3.02892,3.33530,4.11209)
  
  id <- which(a==alpha)
  
  if( !length(id) ){
    warning("alpha value has not been found, using 0.05 instead.",call.=FALSE)
    id <- which(a==0.05)
  }
  
  d2n[id]
}

# the expectation of Z_t(g)
# march.mtd.h.z <- function(mtd,y,t,g){
#   s <- 0
#   for( k in 1:mtd@order ){
#     s <- s+mtd@phi[k]*mtd@Q[1,y@y[t-k],y@y[t]]
#   }
#   if(y@Ncov>0){
#     for(j in 1:y@Ncov)
#       s <- s+mtd@phi[order+j]*mtd@S[[j]][cov[n,t,j],y@y[t]]
#   }
#   mtd@phi[g]*mtd@Q[1,y@y[t-g],y@y[t]]/s
# }
# 
# 
# # the weight coefficient for the g-th lag 
# march.mtd.h.l <- function(mtd,y,g){
#   s <- 0
#   for( i in 1:y@N ){
#     ys <- march.dataset.h.extractSequence(y,i)
#     for( t in march.h.seq(mtd@order+1,ys@N)){
#       s <- s+ march.mtd.h.z(mtd,ys,t,g)
#     }
#   }
#   s
# }

# Matrix containing the estimations of the number of data used to compute each element of Q
# march.mtd.h.n <- function(mtd,y){
#   nki_0 <- array(0,c(mtd@y@K,mtd@y@K))
#   for( k in 1:mtd@y@K ){
#     for( i0 in 1:mtd@y@K){
#       s <- 0
#       for( i in 1:y@N ){
#         ys <- march.dataset.h.extractSequence(y,i)
#         for( t in march.h.seq(mtd@order+1,ys@N) ){
#           if( ys@y[t]==i0 ){
#             for( g in 1:mtd@order ){
#               if( ys@y[t-g]==k)
#                 s <- s+march.mtd.h.z(mtd,ys,t,g,n)
#             }
#           }      
#         }
#       }
#       nki_0[k,i0] <- s
#     }
#   }
#   nki_0
# }

#Matrices containing the estimations of the number of data used to compute each element of Q
march.mtd.h.n <- function(mtd,y,is_mtdg){
  
  #Initialization
  if(is_mtdg==FALSE){
    nki_0 <- array(0,c(mtd@y@K,mtd@y@K))
  }else{
    nki_0 <- array(0,c(mtd@order,mtd@y@K,mtd@y@K))
  }
  numcov <- 0
  if(sum(mtd@MCovar)>0){
    placeCovar <- which(mtd@MCovar==1)
    numcov <- array(0,c(sum(mtd@MCovar),max(y@Kcov[placeCovar]),mtd@y@K))
  }
  
  #Computation of the matrices
  for(n in 1:y@N){
    ys <- march.dataset.h.extractSequence(y,n)
    for( t in march.h.seq(mtd@order+1,ys@N)){
      #Computation of the denominator (see p.9 Confidence Intervals for Markovian Models, Berchtold)
      tot <- march.mtd.h.z.tot(mtd,ys,t,n,is_mtdg)
      for(ord in 1:mtd@order){
        if(is_mtdg==FALSE){
          nki_0[ys@y[t-ord],ys@y[t]] <- nki_0[ys@y[t-ord],ys@y[t]]+mtd@phi[ord]*mtd@Q[1,ys@y[t-ord],ys@y[t]]/tot
        }else{
          nki_0[ord,ys@y[t-ord],ys@y[t]] <- nki_0[ord,ys@y[t-ord],ys@y[t]]+mtd@phi[ord]*mtd@Q[ord,ys@y[t-ord],ys@y[t]]/tot
        }
      }
      if(sum(mtd@MCovar)>0){
        for(i in 1:sum(mtd@MCovar)){
          numcov[i,y@cov[n,t,placeCovar[i]],ys@y[t]] <- numcov[i,y@cov[n,t,placeCovar[i]],ys@y[t]]+mtd@phi[mtd@order+i]*mtd@S[[i]][y@cov[n,t,placeCovar[i]],ys@y[t]]/tot
        }
      }
    }
  }
  list(nki_0=nki_0,numcov=numcov)
}

march.mtd.h.z.tot <- function(mtd,ys,t,n,is_mtdg){
  
  s <- 0
  placeCovar <- which(mtd@MCovar==1)
  
  if(is_mtdg==FALSE){
    for( k in 1:mtd@order ){
      s <- s+mtd@phi[k]*mtd@Q[1,ys@y[t-k],ys@y[t]]
    }
  }else{
    for( k in 1:mtd@order ){
      s <- s+mtd@phi[k]*mtd@Q[k,ys@y[t-k],ys@y[t]]
    }
  }
  
  if(sum(mtd@MCovar)>0){
    for(j in 1:sum(mtd@MCovar))
      s <- s+mtd@phi[mtd@order+j]*mtd@S[[j]][mtd@y@cov[n,t,placeCovar[j]],ys@y[t]]
  }
  s
}


alphat <- function(d,s){
  a <- array(0,c(s@N,2))
  
  a[1,] <- d@RB[,1,s@y[1]]*d@Pi[1,1,]
  for( t in 2:s@N ){
    for( g in 1:d@M ){
      for( i in 1:d@M ){
        a[t,g] <- a[t,g]+a[t-1,i]*d@A[i,g]*d@RB[g,1,s@y[t]]
      }
    }
  }
  a
}

betat <- function(d,s){
  b <- array(0,c(s@N,2))
  
  b[s@N,] <- c(1,1)
  for( t in (s@N-1):1 ){
    for( g in 1:d@M ){
      for( i in 1:d@M ){
        b[t,i] <- b[t,i]+d@A[i,g]*d@RB[g, 1 ,s@y[t+1]]*b[t+1,g]
      }
    }
  }
  b
}

march.dcmm.exp_z <- function(d,alpha,beta,t,g,n){
  L <- sum(alpha[n,])
  
  alpha[t,g]*beta[t,g]/L
}

#
# march.dcmm.test <- function(d){
#   PoA <-  array(0,c(d@M^d@orderHC))
#   PoCt <- array(0,c(d@M^d@orderHC))
#   for( i in 1:d@y@N ){
#     # number of point for A
#     s <- march.dataset.h.extractSequence(d@y,i)
#     
#     a <- march.dcmm.forward(d,s)
#     b <- march.dcmm.backward(d,s)
#     epsilon <- march.dcmm.epsilon(d,s,a$alpha,b$beta,a$l,b$l,a$LL)
#     gamma <- march.dcmm.gamma(d,s,a$alpha,b$beta,a$l,b$l,a$LL,epsilon)
#     
#     PoA <- PoA+colSums(gamma[d@orderHC:(d@y@T[i]),])
#     
#     if( d@M>2 ){
#       PoCt <- PoCt+colSums(gamma[1:(d@M-1),])+colSums(gamma[d@M:d@y@T[i],])  
#     }
#     else{
#       PoCt <- PoCt+gamma[1,]+colSums(gamma[d@M:d@y@T[i],])  
#     }
#   }
#   
#   PoC <- array(0,c(1,d@M))
#   for( i in 0:(d@M-1) ){
#     PoC[i+1] <- sum(PoCt[(i*d@orderHC+1):((i+1)*d@orderHC)])
#   } 
#   list(PoA,PoC)
# }