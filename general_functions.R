### General Functions

# create function to create correlation matrix
Create.Cor.Mat4 <- function(No.Var.H, No.Var.Hsub, rho){
  corrmat <- matrix(NA,nrow=No.Var.H, ncol=No.Var.H)
  rho.vec <- vector()
  rho.vec[1] <- 1
  for (k in 1:(No.Var.H-1)){
    temp <- rho^(k)
    rho.vec <- c(rho.vec, temp)
  }
  rho.vec2 <- rep(rho.vec,times=(No.Var.H/2))
  for (m in 1:(No.Var.H)){
    corrmat[(m:No.Var.H),m] <- rho.vec2[1:((No.Var.H-(m-1)))]
  }
  corrmat[upper.tri(corrmat)] <- t(corrmat)[upper.tri(corrmat)]
  seq <- 1:No.Var.Hsub
  seq2 <- (No.Var.Hsub+1):No.Var.H
  vecd <- c(seq[1], No.Var.Hsub+1,No.Var.Hsub+2,No.Var.Hsub+3)
  for (i in vecd){
    for (j in vecd){
      corrmat[i,j] <- ifelse(i==j, 1, rho)
      corrmat[j,i] <- ifelse(i==j, 1, rho)
    }
  }
  corrmat[1,2] <- corrmat[2,1] <- corrmat[1,3] <- corrmat[3,1] <- corrmat[2,3] <- corrmat[3,2] <- 0
  return(corrmat)
}

## create function to create correlation matrix - simple exchangeable correlation structure
Create.Cor.Mat5 <- function(No.Var.H, No.Var.Hsub, rho){
  corrmat <- matrix(NA,nrow=No.Var.H, ncol=No.Var.H)
  rho.vec <- vector()
  rho.vec[1] <- 1
  for (k in 1:(No.Var.H-1)){
    temp <- rho
    rho.vec <- c(rho.vec, temp)
  }
  rho.vec2 <- rep(rho.vec,times=(No.Var.H/2))
  for (m in 1:(No.Var.H)){
    corrmat[(m:No.Var.H),m] <- rho.vec2[1:((No.Var.H-(m-1)))]
  }
  corrmat[upper.tri(corrmat)] <- t(corrmat)[upper.tri(corrmat)]
  return(corrmat)
}


### function to generate
A.sim<-function(matrix.pi){
  N<-nrow(matrix.pi) # sample size
  K<-ncol(matrix.pi) # treatment options
  if(N<=1 | K<=1) stop("Sample size or treatment options are insufficient!")
  if(min(matrix.pi)<0) stop("Treatment probabilities should not be negative!")
  pis<-apply(matrix.pi,1,sum)
  probs<-matrix(NA,N,K) 
  A<-rep(NA,N) 
  for(i in 1:N){
    probs[i,]<-matrix.pi[i,]/pis[i]
    A[i]<-sample(0:(K-1),1,prob=probs[i,])
  }
  A
}

### function to estimate propensity score
M.propen<-function(A,Xs){
  if(ncol(as.matrix(A))!=1) stop("Cannot handle multiple stages of treatments together!")
  if(length(A)!= nrow(as.matrix(Xs))) stop("A and Xs do not match in dimension!")
  if(length(unique(A))<=1) stop("Treament options are insufficient!")
  class.A<-sort(unique(A)) 
  A <- as.factor(A)
    require(nnet) 
  s.data<-data.frame(A,Xs) 
  model<-capture.output(mlogit<-multinom(A ~., data=s.data))
  s.p<-predict(mlogit,s.data,"probs") 
  if(length(class.A)==2){
    s.p<-cbind(1-s.p,s.p) 
  }
  colnames(s.p)<-paste("pi=",class.A,sep="")
  s.p 
}


mus.AIPW<-function(Y,A,pis.hat,mus.reg){
    ## input one dimensional vector outcome Y
    ## input one dimensional treatment vector A
    ## input a N x R propensity matrix
    ## mus.reg is the the estimated counterfactual mean outcomes (N x R) matrix
    class.A<-sort(unique(A))
    K<-length(class.A) ## represents the number of treatments
    N<-length(A) ## sample size, i.e, number of observations
    if(K<2 | N<2) stop("No multiple treatments or samples!")
    if(ncol(pis.hat)!=K | ncol(mus.reg)!=K | nrow(pis.hat)!=N | nrow(mus.reg)!=N) stop("Treatment, propensity or conditional means do not match!")
    #AIPW estimates
    mus.a<-matrix(NA,N,K) # sets of a blank matrix of N x (number of treaments)
    for(k in 1L:K){ ## computes AIPW estimates
        mus.a[,k]<-((A==class.A[k])*Y/pis.hat[,k])+(1-((A==class.A[k])/pis.hat[,k]))*mus.reg[,k]
    }
    mus.a ## outputs the AIPW estimates OF MU'S.
}


# choose the optimal treatment with given mu's matrix
Opt.A<-function(A,mus.hat){
    class.A<-sort(unique(A))
    if(length(class.A)==1){
        trt.opt1<-class.A
        Ey.opt1<-mean(mus.hat)
    } else{
        if(length(A)!= nrow(mus.hat) || length(class.A)!= ncol(mus.hat)){
            stop("Treatment options and mean matrix dimension do not match!")
        }
        
        # pick a single best treatment for all patients
        c.means<-apply(mus.hat,2,mean) ##
        Ey.opt1<-max(c.means)
        trt.opt1<-class.A[which(c.means==Ey.opt1)]
    }
    outs<-list(Ey.opt1,trt.opt1)
    names(outs)<-c("Ey.opt1","trt.opt1")
    outs
}

# combine matrices with different dimensions, add NA for additional rows/columns
combine.mat<-function(m1,m2,by="column"){
    if(is.matrix(m1)==F){
        m1a <- matrix(as.factor(m1),nrow=1, ncol=length(m1))}
    else {
        m1a <- m1
    }
    m2a <- matrix(as.character(m2), byrow=FALSE, nrow=nrow(m2), ncol=ncol(m2))
    nrow1<-nrow(m1a);ncol1<-ncol(m1a)
    nrow2<-nrow(m2);ncol2<-ncol(m2)
    if(by=="column"){
        combine<-matrix(NA,max(nrow1,nrow2),ncol1+ncol2)
        combine[1:nrow1,1:ncol1]<-m1a
        combine[1:nrow2,(ncol1+1):(ncol1+ncol2)]<-m2a
    }
    if(by=="row"){
        combine<-matrix(NA,nrow1+nrow2,max(ncol1,ncol2))
        combine[1:nrow1,1:ncol1]<-m1
        combine[(nrow1+1):(nrow1+nrow2),1:ncol2]<-m2
    }
    combine
}


# splix data by X to fit child nodes, calculate new means for each child node
Split.X.modified<-function(X,A,mus.hat,minsplit=20){
    n<-length(X)
    X.val<-unique(X)
    X.val <- sort(X.val)
    n.X.val<-length(X.val) ## caution number of levels in factor!!
    class.A<-sort(unique(A))
    
    if(n < 2*minsplit || n.X.val<2L || length(class.A)<2L) return(NULL)
    ## removed the n.X.val==2L b/c this appeared to be causing an issue
    if(is.numeric(X)==T || is.ordered(X)==T ){ #|| n.X.val==2L
        X.val<-sort(X.val)
        # reduce computation by using quantiles
        if(n.X.val>100L){
            X.val<-quantile(X,1:100/100)
            n.X.val<-100L
        }
        Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.val-1)
        for(i in 1L:(n.X.val-1)){
            left<-which(X<=X.val[i]) ## here's where we need to correct this for factors
            if(length(left)>=minsplit && length(left)<=n-minsplit){
                
                left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
                right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
                
                trt.left.X[i]<-left.est$trt.opt1 ## why is it always the same Ey.opt1 value??
                trt.right.X[i]<-right.est$trt.opt1
                Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1##
                ## the above is the purity measure, i.e., pg 10 in Yebin's paper
                ## just like a weighted average
            }
        }
        # pick the best split of X
        if(sum(!is.na(Ey.opt1.X))>0L){
            mEy.opt1<-max(Ey.opt1.X, na.rm=T)
            cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
            X.cutoff<-X.val[cutoff1]
            trt.L<-trt.left.X[cutoff1]
            trt.R<-trt.right.X[cutoff1]
            
            output<-data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
            names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
        } else{
            return(NULL)
        }
        
    }
    ## below is optin for is.numeric(X) == F!! check this out!!
    ## remove the n.X.val>2L condition here, at least temporarily.
    ## this appears to correct the factor issue.
    else if(is.numeric(X)==F && is.ordered(X)==F ){    # && n.X.val>2L
        n.X.combo<-2^(n.X.val-1)-1
        X.combo<-combn(X.val,1) ## does it matter if they are not ordered?
        X.combo <- sort(X.combo)
        if(n.X.val>3L && n.X.val%%2==1L){ ## could be a problem with combine.mat
            ## need to address this for factor levels >=7
            if(n.X.val>6){ ## changed from length(X.combo) to length(n.X.val)
                X.combo <- X.combo
                n.X.combo <- length(X.combo)
            }
            else if(n.X.val<=6){
                for(k in 2L:(n.X.val-1)/2) X.combo<-combine.mat(X.combo,combn(as.factor(X.val),k), by="column") #this line will probably need
            } else {
                return("error2")
            }
        }
        if(n.X.val>3L && n.X.val%%2==0L){
            if(n.X.val>6){ ## changed from length(X.combo) to length(n.X.val)
                X.combo <- X.combo
                n.X.combo <- length(X.combo)
            }
            else if(n.X.val<=6 & n.X.val>=5){
                
                for(k in 2L:(n.X.val/2)){
                    if(k<(n.X.val/2)) {
                        X.combo1<-combine.mat(m1=X.combo,m2=combn(X.val,k), by="column")
                    }##
                    else  ## problem here is that there are no observations for factor level 3
                    ## remove this from ifelse "if(k==(n.X.val/2))"
                    ## I added "as.factor" so as not to lose the factor labels when merging these
                    ## but I may also need to add a section for where it isn't a factor
                    {temp.mat<-combn(X.val[-1],k-1)
                        first.row<-rep(X.val[1],ncol(temp.mat))
                        #X.combo.mat <- matrix(as.factor(X.combo),nrow=1,ncol=length(X.combo)) ## added this row
                        ## NEED TO MAKE SURE IT STAYS AS A FACTOR LEVEL AND NOT AS NUMERIC START HERE!! 2:19PM 11/8/2018
                        X.combo2<-rbind(as.factor(first.row),as.factor(temp.mat)) ## why does rbind lose the labels??
                        #X.combo2<-combine.mat(X.combo1,rbind(as.factor(first.row),as.factor(temp.mat))) ## why does rbind lose the labels??
                        ## may need to bring X.combo2 outside of the bracket
                    }
                    
                }
                ## added the X.combo here
                X.combo <- combine.mat(X.combo1, X.combo2, by="column") ## this seems fine
            }
            else { ##  removed the else if(n.X.val==4) for just "else
                X.combo<-combine.mat(m1=X.combo,m2=combn(X.val,2), by="column")
            }
            #X.combo <- combine.mat(X.combo1, X.combo2, by="column")
        }
        
        Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.combo)
        ### now need to address whether n.X.val > 6
        if(n.X.val>6L){
            for(i in 1L:n.X.combo){
                left<-which(X %in% X.combo[i]) ## changed from X.combo[,i] to X.combo[i] ## changed back to X.combo[,i] back to X.combo[i]
                if(length(left)>=minsplit && length(left)<=n-minsplit){
                    
                    left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
                    right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
                    
                    trt.left.X[i]<-left.est$trt.opt1
                    trt.right.X[i]<-right.est$trt.opt1
                    Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
                }
            }
        }
        else if(n.X.val<=6L){
            for(i in 1L:n.X.combo){ ## problem getting the correct values for subset left
                if(is.matrix(X.combo)==F){
                    left<-which(X %in% X.combo[i]) ## changed from X.combo[,i] to X.combo[i] ## changed back to X.combo[,i] back to X.combo[i]
                    if(length(left)>=minsplit && length(left)<=n-minsplit){
                        
                        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
                        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
                        
                        trt.left.X[i]<-left.est$trt.opt1
                        trt.right.X[i]<-right.est$trt.opt1
                        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
                    }
                }
                else if(is.matrix(X.combo)==T){
                    left<-which(X %in% X.combo[,i]) ## changed from X.combo[,i] to X.combo[i] ## changed back to X.combo[,i] back to X.combo[i]
                    if(length(left)>=minsplit && length(left)<=n-minsplit){
                        
                        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
                        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
                        
                        trt.left.X[i]<-left.est$trt.opt1
                        trt.right.X[i]<-right.est$trt.opt1
                        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
                    }
                }
                else{
                    return("error4")
                }
            }
        } ## end n.X.val<=6
        else{
            return("error3")
        }
        # pick the best split of X
        if(sum(!is.na(Ey.opt1.X))>0L && length(Ey.opt1.X)>=2L){
            mEy.opt1<-max(Ey.opt1.X, na.rm=T)
            cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
            if(is.matrix(X.combo)==F){
                X.subset<-X.combo[cutoff1] ## changed from X.combo[,cutoff1] to X.combo[cutoff1]
                # change a vector into a single string while removing NA's
                X.subset<-paste(X.subset[!is.na(X.subset)], collapse=" ")
                trt.L<-trt.left.X[cutoff1]
                trt.R<-trt.right.X[cutoff1]
            }
            else if(is.matrix(X.combo)==T){
                X.subset<-X.combo[,cutoff1] ## changed from X.combo[,cutoff1] to X.combo[cutoff1]
                # change a vector into a single string while removing NA's
                X.subset<-paste(X.subset[!is.na(X.subset)], collapse=" ")
                trt.L<-trt.left.X[cutoff1] ## fix this
                trt.R<-trt.right.X[cutoff1]
            } else{
                return("error1")
            }
            
            output<-data.frame(X.subset, mEy.opt1, trt.L, trt.R)
            names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
        } else if(sum(!is.na(Ey.opt1.X))>0L && length(Ey.opt1.X)<2L){
            mEy.opt1<-max(Ey.opt1.X, na.rm=T)
            cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
            X.subset<-X.combo[cutoff1] ## changed to [cutoff1] from [,cutoff1]
            # change a vector into a single string while removing NA's
            X.subset<-paste(X.subset[!is.na(X.subset)], collapse=" ")
            trt.L<-trt.left.X[cutoff1]
            trt.R<-trt.right.X[cutoff1]
            
            output<-data.frame(X.subset, mEy.opt1, trt.L, trt.R)
            names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
        }
        else{
            return(NULL)
        }
        #} ## may need to add another bracket here - nope looks like it's okay
    } ## end bracket for "if numeric=false", i.e., factor
    else{
        return(NULL)
    }
    return(output)
}

#######################################################
# predit optimal treatment using the output from DTRtree
predict.DTR.modified<-function(treeout,newdata){
    n<-nrow(newdata)
    predicts<-rep(NA,n)
    
    # treeout is supposed to be a matrix
    # if there is no split
    #if(length(treeout)==5){ ## kas removed this and changed to dim(treeout)[1] == 1
    if(dim(treeout)[1]==1){
        predicts<-rep(treeout[1,5],n)## also changed this from treeout[5] to treeout[1,5]
    } else{ # if there are splits
        treeout<-as.data.frame(treeout)
        newdata<-as.data.frame(newdata)
        
        for(i in 1:n){
            nd<-1
            while(is.na(treeout$trt[treeout$node==nd])){
                ## add factor level here
                if(is.numeric(newdata[i,treeout$X[treeout$node==nd]])==F){
                    if(newdata[i,treeout$X[treeout$node==nd]] == treeout$cutoff[treeout$node==nd]){
                        nd=2*nd
                    } else{
                        nd=2*nd+1
                    }
                }
                if(is.numeric(newdata[i,treeout$X[treeout$node==nd]])==T){
                    if(newdata[i,treeout$X[treeout$node==nd]] <= treeout$cutoff[treeout$node==nd]){
                        nd=2*nd
                    } else{
                        nd=2*nd+1
                    }
                }
            } ## end here
            predicts[i]<-treeout$trt[treeout$node==nd]
        }
    }
    return(predicts)
}

predict.DTR.modified2<-function(treeout,newdata){
    n<-nrow(newdata)
    predicts<-rep(NA,n)
    
    # treeout is supposed to be a matrix
    # if there is no split
    #if(length(treeout)==5){ ## kas removed this and changed to dim(treeout)[1] == 1
    if(dim(treeout)[1]==1){
        predicts<-rep(treeout[1,5],n)## also changed this from treeout[5] to treeout[1,5]
    } else{ # if there are splits
        treeout<-as.data.frame(treeout)
        newdata<-as.data.frame(newdata)
        treeout$X <- as.numeric(treeout$X) ## added this in modified2 version.
        treeout$cutoff <- as.numeric(treeout$cutoff)
        
        for(i in 1:n){
            nd<-1
            while(is.na(treeout$trt[treeout$node==nd])){
                ## add factor level here
                # if(is.numeric(newdata[i,treeout$X[treeout$node==nd]])==F){ ### COMMENTED THIS OUT IN MODIFIED2 VERSION.
                #   if(newdata[i,treeout$X[treeout$node==nd]] == treeout$cutoff[treeout$node==nd]){
                #     nd=2*nd
                #   } else{
                #     nd=2*nd+1
                #   }
                # }
                #if(is.numeric(newdata[i,treeout$X[treeout$node==nd]])==T){
                if(newdata[i,treeout$X[treeout$node==nd]] <= treeout$cutoff[treeout$node==nd]){
                    nd=2*nd
                } else{
                    nd=2*nd+1
                }
                #}
            } ## end here
            predicts[i]<-treeout$trt[treeout$node==nd]
        }
    }
    return(predicts)
}
