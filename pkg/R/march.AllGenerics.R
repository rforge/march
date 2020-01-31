# This file is part of March
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
#


###############################################################################
# Show methods here, are the implementation of generic methods to print
# models to the standard output. These methods are implemented for all
# models (indep, mc, mtd, dcmm).
###############################################################################

# show method for model (abstract) object
march.model.show <- function(object){
  cat("Log-likelihood : ",object@ll,"\n")
  cat("Number of data : ",object@dsL,"\n")
}

# show method for indep object
march.indep.show <- function(object){
  cat(march.name(object))
  cat("\n\nProbability distribution :\n")
  
  s <- array(0,c(1,length(object@indP)))
  colnames(s) <- 1:object@y@K
  
  s[1,] <- object@indP
  print(s)
  
  cat("\nFrequency distribution :\n")
  
  s[1,] <- object@indC
  print(s)
  cat("\n")
  
  march.model.show(object)
}

# show method for mc object
march.mc.show <- function(object){
  cat(march.name(object))
  cat("\n\nTransition matrix of order", object@order, ": \n")
  
  s <- march.h.mc.printableMatrix(object@RC,object@order,object@y@K)
  print(s)
  cat("\n")
  
  cat("Crosstable : \n")
  s <- march.h.mc.printableMatrix(object@RT,object@order,object@y@K)
  print(s)
  cat("\n")

  march.model.show(object)
}

# show method for mtd object
march.mtd.show <- function(object){
  placeCovar <- which(object@MCovar==1)
  
  cat(march.name(object))
  cat("\n\n")
  cat("High-order transition matrix : \n")
  s <- march.cov.h.mc.printableMatrix(object@RA,object@order,object@y@K,object@y@Kcov[placeCovar],sum(object@MCovar))
  print(s)
  cat("\n")
  
  # MTD model
  cat("Vector of weights : \n")
  print(object@phi)
  cat("\n")
  
  for (i in 1:(dim(object@Q)[1])){
    if (i==1){
      if (dim(object@Q)[1]==1){
        cat("Transition matrix : \n")
      } else {
        cat("Transition matrix, lag",i,": \n")        
      }
    } else {
      cat("Transition matrix, lag",i,": \n")        
    }
    print(object@Q[i,,])
    cat("\n")    
  }
  
  if (length(object@phi)>object@order){
    for (i in 1:sum(object@MCovar)){
      cat("Transition matrix for covariate",i,": \n")        
      print(object@S[i])
    }
  }
  
  march.model.show(object)
}

march.dcmm.show <- function(object){
	  placeACovar <- which(object@AMCovar==1)
	  placeCCovar <- which(object@CMCovar==1)
	  
	  AtmCovar <- 1
  	if(sum(object@AMCovar>0)){
  		AtmCovar <-prod(object@y@Kcov[placeACovar])
  	}

	  CtmCovar <- 1
	  if(sum(object@CMCovar>0)){
	    CtmCovar <-prod(object@y@Kcov[placeCCovar])
	  }
	  
  	cat(march.name(object),"\n\n")
  	
  	if (object@M==1 & sum(object@AMCovar)==0){
  	  cat("No hidden level \n\n")
  	}else{
  	  cat("Hidden level: \n")
  	  cat("=============  \n\n")
  	  
  	  for( i in 1:object@orderHC){
  	    cat("Distribution of hidden state",i,": \n")
  	    for( j in 1:(AtmCovar*object@M^(i-1))){
  	      cat(object@Pi[j,,i])
  	      cat("\n")
  	    }
  	    cat("\n")
  	  }
  	  
  	  cat("Hidden transition matrix : \n")
  	  if(object@orderHC==0){
  	    s <- march.cov.h.mc.printableMatrix(object@A[1:AtmCovar,,drop=FALSE],0,object@M,object@y@Kcov[placeACovar],sum(object@AMCovar))
  	  }else{
  	    s <- march.cov.h.mc.printableMatrix(march.dcmm.cov.h.compactA(object),object@orderHC,object@M,object@y@Kcov[placeACovar],sum(object@AMCovar))
  	  }
  	  print(s)
  	  cat("\n") 	  
  	  
  	  
  	  if(sum(object@AMCovar)>0 | (object@Amodel=="mtd" & object@orderHC>1) | (object@Amodel=="mtdg" & object@orderHC>1)){
  	    # Display the modeling of the hidden transition matrix
  	    # Rem: Covariates are allowed at the hidden level only when the hidden order is >0
  	    #      and there are at least two hidden states
  	    
  	    cat("Modeling of the hidden transition matrix : \n\n")
  	    
  	    cat("Vector of weights : \n")
  	    print(object@APhi)
  	    cat("\n")
  	    
  	    for (i in 1:(dim(object@AQ)[1])){
  	      if (i==1){
  	        if (dim(object@AQ)[1]==1){
  	          cat("Transition matrix : \n")
  	        } else {
  	          cat("Transition matrix, lag",i,": \n")        
  	        }
  	      } else {
  	        cat("Transition matrix, lag",i,": \n")        
  	      }
  	      print(object@AQ[i,,])
  	      cat("\n")    
  	    }
  	    
  	    if (sum(object@AMCovar)>0){
  	      for (i in 1:sum(object@AMCovar)){
  	        cat("Transition matrix for covariate",i,": \n")        
  	        print(object@ATCovar[i])
  	        cat("\n")  
  	      }
  	    }
  	    
  	  }
  	  
  	  
  	}
  	
  	cat("Visible level : \n")
  	cat("=============== \n\n")
  	
   	for( i in 1:object@M ){
   	  if (object@M>1){
   	    cat("Model corresponding to hidden state", i, ":\n")
   	    cat("------------------------------------- \n\n")
   	  }
   	  
   	  if(object@orderVC==0){
   	    s <- march.cov.h.mc.printableMatrix.orderVC0(matrix(object@RB[,,i],CtmCovar,object@y@K),object@y@K,object@y@Kcov[placeCCovar],sum(object@CMCovar))
   	  }else{
        s <- march.cov.h.mc.printableMatrix(object@RB[,,i],object@orderVC,object@y@K,object@y@Kcov[placeCCovar],sum(object@CMCovar))
   	  }
   	  print(s)
   	  cat("\n")
   	  
   	  if(sum(object@CMCovar)>0 | (object@Cmodel=="mtd" & object@orderVC>1) | (object@Cmodel=="mtdg" & object@orderVC>1)){
   	    # Display the modeling of the visible transition matrix

   	    cat("Modeling of the visible transition matrix : \n\n")
   	    
   	    cat("Vector of weights : \n")
   	    print(object@CPhi[1,,i])
   	    cat("\n")
   	    
   	    for (j in 1:(dim(object@CQ)[1])){
   	      if (j==1){
   	        if (dim(object@CQ)[1]==1){
   	          cat("Transition matrix : \n")
   	        } else {
   	          cat("Transition matrix, lag",j,": \n")        
   	        }
   	      } else {
   	        cat("Transition matrix, lag",j,": \n")        
   	      }
   	      print(object@CQ[j,,,i])
   	      cat("\n")    
   	    }
   	    
   	    if (sum(object@CMCovar)>0){
   	      for (j in 1:sum(object@CMCovar)){
   	        cat("Transition matrix for covariate",j,": \n")  
   	        tmpcov <- as.array(object@CTCovar[[j]])
   	        print(tmpcov[,,i])
   	        cat("\n")  
   	      }
   	    }
   	    
   	  }
   	  
   	  
   	} 
  	

  	march.model.show(object)
}

march.AIC.show <- function(object){
  cat("Number of parameters : ",object@nbParams,"\n")
  cat("AIC: ",object@AIC,"\n")
}

march.BIC.show <- function(object){
  cat("Number of parameters : ",object@nbParams,"\n")
  cat("BIC: ",object@BIC,"\n")
}

# this part describes how a call to show method (print) should be 
# redirected to the rigth method, depending on the object considered.
setMethod(f="show",signature="march.Indep",definition=march.indep.show)
setMethod(f="show",signature="march.Mc",definition=march.mc.show)
setMethod(f="show",signature="march.Mtd",definition=march.mtd.show)
setMethod(f="show",signature="march.Dcmm",definition=march.dcmm.show)
setMethod(f="show",signature="march.AIC",definition=march.AIC.show)
setMethod(f="show",signature="march.BIC",definition=march.BIC.show)


###############################################################################
# nbParams methods here, are the implementation of generic methods to obtain
# the number of parameters that have been used during the model construction.
# These methods are implemented for all models (indep, mc, mtd, dcmm).
###############################################################################

# nbParams method for model (abstract) object
march.model.nbParams <- function(object){
  0
}

# nbParams method for mc object
march.mc.nbParams <- function(object){
  object@y@K^object@order*(object@y@K-1)-object@nbZeros+length(which(rowSums(object@RT)==0))
}

# nbParams method for indep object
march.indep.nbParams <- function(object){
  object@y@K-1
}

# nbParams method for mtd object
march.mtd.nbParams <- function(object){
  nparam <- 0
  placeCovar <- which(object@MCovar==1)
  
  #Number of parameters of the matrix Q
  for(i in 1:dim(object@Q)[1]){
    if(object@phi[i]>0){
      nparam <- nparam+object@y@K*(object@y@K-1)+sum(rowSums(object@Q[i,,])==0)-sum(object@Q[i,,]==0)
    }
  }
  
  #Number of parameters in the vector of lags
  nparam <- nparam+length(object@phi)-1-sum(object@phi==0)
  
  #Number of parameters of the covariates transition matrices
  if(sum(object@MCovar)>0){
    for(i in 1:sum(object@MCovar)){
      if(object@phi[object@order+i]>0){
        nparam <- nparam+object@y@Kcov[placeCovar[i]]*(object@y@K-1)-sum(object@S[[i]]==0)+sum(rowSums(object@S[[i]])==0)
      }
    }
  }
  nparam
}

# nbParams method for dcmm object
march.dcmm.nbParams <- function(object){
  
  placeACovar <- which(object@AMCovar==1)
  placeCCovar <- which(object@CMCovar==1)
  AtmCovar <- 1
  if(sum(object@AMCovar>0)){
    for (i in 1:sum(object@AMCovar)){
      AtmCovar <- AtmCovar*object@y@Kcov[placeACovar[i]]
    }
  }
  #Pi
  nbparpi <- 0
  
  if(object@M>1 & object@orderHC>0){
    for(t in 1:object@orderHC){
      nbparpi <- nbparpi+object@M^(t-1)*AtmCovar*(object@M-1)
      for(i in 1:(object@M^(t-1)*AtmCovar)){
        for(j in 1:object@M){
          if(object@Pi[i,j,t]==0){
            nbparpi <- nbparpi-1
          }
        }
        if(sum(object@Pi[i,,t])==0){
          nbparpi <- nbparpi+1
        }
      }
    }
  }
  
  #A
  nbparA <- 0
  
  if(object@M>1){
    if(object@Amodel=="complete" & object@APhi[1,1]!=0){
      if(object@orderHC==0){
        nbparA <- object@M-1-sum(object@AQ[1,1,]==0)
      }else{
        nbparA <- object@M^object@orderHC*(object@M-1)-sum(object@AQ[1,,]==0)+sum(rowSums(object@AQ[1,,])==0)
      }
    }else if(object@Amodel=="mtd" & sum(object@APhi[1,1:object@orderHC])>0){
      nbparA <- object@M*(object@M-1)-sum(object@AQ[1,,]==0)+sum(rowSums(object@AQ[1,,])==0)
    }else if(object@Amodel=="mtdg" & sum(object@APhi[1,1:object@orderHC])>0){
      for(ord in 1:object@orderHC){
        if(object@APhi[1,ord]!=0){
          nbparA <- nbparA+object@M*(object@M-1)-sum(object@AQ[ord,,]==0)+sum(rowSums(object@AQ[ord,,])==0)
        }
      }
    }
    
    if(sum(object@AMCovar)>0){
      if(object@Amodel=="complete"){
        CCov <- 1
      }else{
        CCov <- object@orderHC
      }
      
      for(ord in 1:sum(object@AMCovar)){
        if(object@APhi[1,CCov+ord]!=0){
          KC <- object@y@Kcov[placeACovar[ord]]
          nbparA <- nbparA+KC*(object@M-1)-sum(object@ATCovar[[ord]]==0)+sum(rowSums(object@ATCovar[[ord]])==0)
        }
      }
    }
    
    #Independent parameters in APhi
    
    nbparA <- nbparA+length(object@APhi[1,])-1-sum(object@APhi[1,]==0)
  }
  
  #Number of parameters visible process
  nbparC <- 0
  for(state in 1:object@M){
    if(object@Cmodel=="complete" & object@CPhi[1,1,state]!=0){
      if(sum(object@CMCovar)==0){
        if(object@orderVC==0 & sum(object@CMCovar)==0){
          nbparC <- nbparC+object@y@K^object@orderVC*(object@y@K-1)-sum(object@RB[,,state]==0)+sum(sum(object@RB[,,state])==0)
        }else{
          nbparC <- nbparC+object@y@K^object@orderVC*(object@y@K-1)-sum(object@RB[,,state]==0)+sum(rowSums(object@RB[,,state])==0)
        }
      }else{
        nbparC <- nbparC+object@y@K^object@orderVC*(object@y@K-1)-sum(object@CQ[1,,,state]==0)+sum(rowSums(object@CQ[1,,,state])==0)
      }
    }else if(object@Cmodel=="mtd" & sum(object@CPhi[1,1:object@orderVC,state])!=0){
      nbparC <- nbparC+object@y@K*(object@y@K-1)-sum(object@CQ[1,,,state]==0)+sum(rowSums(object@CQ[1,,,state])==0)
    }else if(object@Cmodel=="mtdg" & sum(object@CPhi[1,1:object@orderVC,state])!=0){
      for(ord in 1:object@orderVC){
        if(object@CPhi[1,ord,state]!=0){
          nbparC <- nbparC+object@y@K*(object@y@K-1)-sum(object@CQ[ord,,,state]==0)+sum(rowSums(object@CQ[ord,,,state])==0)
        }
      }
    }
    
    if(sum(object@CMCovar)>0){
      if(object@Cmodel=="complete"){
        CCov <- 1
      }else{
        CCov <- object@orderVC
      }
      
      for(ord in 1:sum(object@CMCovar)){
        if(object@CPhi[1,CCov+ord,state]!=0){
          KC <- object@y@Kcov[placeCCovar[ord]]
          nbparC <- nbparC+KC*(object@y@K-1)-sum(object@CTCovar[[ord]][,,state]==0)+sum(rowSums(object@CTCovar[[ord]][,,state])==0)
        }
      }
    }
    
    #Independant parameters in CPhi
    nbparC <- nbparC+length(object@CPhi[1,,state])-1-sum(object@CPhi[1,,state]==0)
  }
  
  Nbparam <- nbparpi+nbparA+nbparC
  
  Nbparam
}
# This part create the generic method and describe how a call to this generic 
# has to be redirected to the rigth method, according to the considerd object.
setGeneric(name="march.nbParams",def=function(object)march.model.nbParams(object))
setMethod(f="march.nbParams",signature="march.Indep",definition=march.indep.nbParams)
setMethod(f="march.nbParams",signature="march.Mc",definition=march.mc.nbParams)
setMethod(f="march.nbParams",signature="march.Mtd",definition=march.mtd.nbParams)
setMethod(f="march.nbParams",signature="march.Dcmm",definition=march.dcmm.nbParams)



#' march.Model name.
#' 
#' Generate a name for the march model contained in the given \emph{object}.
#'
#' @param object contains the name of the model(Independence model, MTD,...). 
#' @author Ogier Maitre & Andre Berchtold
#' @example tests/examples/march.name.example.R
#' @export
march.name <- function(object){}

march.model.name <- function(object){ "abstract model" } # should never be called
march.indep.name <- function(object){ "Independence" }
march.mc.name <- function(object){ sprintf("MC(%d)",object@order) }
march.mtd.name <- function(object){
  if( dim(object@Q)[1]==1 ){
    sprintf("MTD(%d)",object@order)
  }
  else{
    sprintf("MTDg(%d)",object@order)
  }
  
}

march.dcmm.name <- function(object){
  if( object@orderVC==0 ){ sprintf("Hmm(%d)",object@orderHC) }
  else{ sprintf("Dcmm(%d,%d)",object@orderHC,object@orderVC) }
}


#This part create the generic method and describe how a call to this generic
#has to be redirected to the right method, according to the considered object.
setGeneric(name="march.name",def=function(object)march.model.name(object))
setMethod(f="march.name",signature=signature(object="march.Indep"),definition=march.indep.name)
setMethod(f="march.name",signature=signature(object="march.Mc"),definition=march.mc.name)
setMethod(f="march.name",signature=signature(object="march.Mtd"),definition=march.mtd.name)
setMethod(f="march.name",signature=signature(object="march.Dcmm"),definition=march.dcmm.name)

#quote=FALSE,digits = getOption("digits")

#' march.Model Summary.
#'
#' Print key values for the current model.
#' 
#' @param object can contain the results of any model computed  using march
#' @param ... should indicate any additional parameter passed to the function
#' @author Ogier Maitre & Andre Berchtold
#' @export
march.summary <- function(object,...){
  v <- array(NA,c(1,4))
  rownames(v) <- march.name(object)
  colnames(v) <- c(gettext("ll"),gettext("param"),gettext("BIC"),gettext("AIC"))
  v[1,]<-c(object@ll,march.nbParams(object),march.BIC(object),march.AIC(object))
  
  v
}

#setMethod('summary',"march.Model",march.summary)
#setMethod(f="summary",signature=signature(object="march.Mc"),definition=march.summary)
# setMethod(f="Summary",signature=signature(object="march.Mtd"),definition=march.summary)
# setMethod(f="Summary",signature=signature(object="march.Dcmm"),definition=march.summary)