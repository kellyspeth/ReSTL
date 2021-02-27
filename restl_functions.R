### REST-L Functions

Reg.mu.modified.Hsub<-function(Y,As,HsubC, Hsub){
  if(nrow(as.matrix(As))!=nrow(as.matrix(HsubC))) stop("Treatment and Covariates do not match in dimension!")
  Ts<-ncol(as.matrix(As)) # number of stages
  N<-nrow(as.matrix(As)) # number of observations
  if(Ts<0 | Ts>3) stop("Only support 1 to 3 stages!")
  
  if(Ts==1L){ ## if there is one stage only
    A1<-as.matrix(As)[,1] 
    A1<-as.factor(A1) 
    KT<-length(unique(A1)) 
    if(KT<2) stop("No multiple treatment options!")
    RegModel<-lm(Y ~ Hsub*A1 + HsubC) 
    mus.reg<-matrix(NA,N,KT) 
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(HsubC, Hsub,A1=factor(rep(sort(unique(A1))[k],N))))
  }
  if(Ts==2L){
    A1<-as.matrix(As)[,1];A2<-as.matrix(As)[,2]
    A1<-as.factor(A1);A2<-as.factor(A2)
    KT<-length(unique(A2))
    if(KT<2) stop("No multiple treatment options!")
    
    RegModel<-lm(Y ~ (Hsub + A1)*A2 + HsubC)
    
    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(Hsub,HsubC,A1,A2=factor(rep(sort(unique(A2))[k],N))))
  }
  if(Ts==3L){
    A1<-as.matrix(As)[,1];A2<-as.matrix(As)[,2];A3<-as.matrix(As)[,3]
    A1<-as.factor(A1);A2<-as.factor(A2);A3<-as.factor(A3)
    KT<-length(unique(A3))
    if(KT<2) stop("No multiple treatment options!")
    
    RegModel<-lm(Y ~ (H + A1 + A2)*A3)
    
    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1,A2,A3=factor(rep(sort(unique(A3))[k],N))))
  }
  
  output<-list(mus.reg, RegModel)
  names(output)<-c("mus.reg","RegModel")
  output 
}

### pick the best X for split
best.Hsub<-function(Hsub,A,mus.hat,minsplit=20){
    p<-ncol(Hsub)
    output<-as.data.frame(matrix(NA,p,5))
    output[,1]<-1:p
    colnames(output)<-c("X","X.subset","mEy.opt1","trt.L","trt.R")
    
    for(i in 1:p){
        split.i<-Split.X.modified(X=Hsub[,i],A=A,mus.hat=mus.hat,minsplit=minsplit) ## error with SITEID!! Many factors!!!
        
        if(!is.null(split.i)){
            if(is.numeric(split.i$X.subset)==T){
                output[i,-1]<-split.i ## another issue here with factors (converts to different value)
            } else if(is.numeric(split.i$X.subset)==F){
                #output[i,2] <- round(as.numeric(paste(split.i$X.subset)),0)
                output[i,2] <- paste(split.i$X.subset)
                output[i,3] <- split.i$mEy.opt1
                output[i,4] <- split.i$trt.L
                output[i,5] <- split.i$trt.R
            }
            else{
                return("error3a")
            }
        }
    }
    if(sum(!is.na(output$mEy.opt1))>0L){
        max.p<-which(output$mEy.opt1==max(output$mEy.opt1,na.rm=T))[1]
        opt.output<-output[max.p,]
        if(opt.output$trt.L==opt.output$trt.R){
            return(NULL)
        } else{
            return(opt.output)
        }
    } else{
        return(NULL)
    }
}


# DTRtree
# input: outcome Y, treatment A, covariate history H, propensity pis.hat
# lambda.pct is minimum purity improvement as a percent of the estimated counterfactaul mean at root node without splitting
# minsplit is minimum node size

DTRtree.modified.HSUB<-function(Y,A,HsubC, Hsub,pis.hat=NULL,m.method="AIPW",mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F){
  n<-length(Y)
  I.node<-rep(1,n)
  class.A<-sort(unique(A))
  output<-matrix(NA,1,5)
  output <- as.data.frame(output)
  colnames(output)<-c("node","X","cutoff","mEy","trt")
  
  if(m.method[1]=="AIPW"){
    if(is.null(pis.hat)){
      pis.hat<-M.propen(A=A,Xs=cbind(Hsub, HsubC))
      if(sum(pis.hat == 1) >= 1){
        stop("At least one Pis.hat == 1")
      }
    } 
    if(is.null(mus.reg)) mus.reg<-Reg.mu.modified.Hsub(Y=Y,As=A,HsubC=HsubC, Hsub=Hsub)$mus.reg
    mus.hat<-mus.AIPW(Y=Y,A=A,pis.hat=pis.hat,mus.reg=mus.reg)
  }
  else{
    stop("The method for estimating conditional means is not available!")
  }

  root<-Opt.A(A,mus.hat) 
  Ey0<-root$Ey.opt1
  lambda<-abs(Ey0)*lambda.pct
  
  for(k in 1L:depth){ 
    temp <- as.data.frame(matrix(NA,2^k,5))
    colnames(temp) <- c("node","X","cutoff","mEy","trt")
    output <- rbind.data.frame(output,temp)
    output[,1]<-1L:(2^(k+1)-1)
    if(k==1L){
      if(lookahead){
        best.H.1<-best.H.lh(H=cbind(Hsub),A=A,mus.hat=mus.hat,minsplit=0.15*n)
      } else{
        best.H.1<-best.Hsub(Hsub=Hsub,A=A,mus.hat=mus.hat,minsplit=minsplit) ##start
      }
      if(is.null(best.H.1)==F && best.H.1$mEy.opt1>Ey0+lambda){
        if(is.numeric(Hsub[,as.numeric(paste(best.H.1$X))])==T){ 
          output[k,-1]<-c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1, NA) 
          I.node[I.node==k & Hsub[,best.H.1$X] <= best.H.1$X.subset]<-2*k
          output[2*k,-1]<-c(NA,NA,NA,best.H.1$trt.L)
          I.node[I.node==k & Hsub[,best.H.1$X] > best.H.1$X.subset]<-2*k+1
          output[2*k+1,-1]<-c(NA,NA,NA,best.H.1$trt.R)
        } 
        else if(is.numeric(Hsub[,as.numeric(paste(best.H.1$X))])==F){
          output[k,2] <- round(as.numeric(paste(best.H.1$X)),0)
          output[k,3] <- paste(best.H.1$X.subset)#
          output[k,4] <- as.numeric(best.H.1$mEy.opt1)
          output[k,5] <- NA
          I.node[I.node==k & Hsub[,as.numeric(paste(best.H.1$X))] %in% unlist(strsplit(paste(best.H.1$X.subset),split= " "))]<-2*k
          output[2*k,-1]<-c(NA,NA,NA,best.H.1$trt.L)
          I.node[I.node==k & !(Hsub[,as.numeric(paste(best.H.1$X))] %in% unlist(strsplit(paste(best.H.1$X.subset),split= " ")))]<-2*k+1
          output[2*k+1,-1]<-c(NA,NA,NA,best.H.1$trt.R)
        } else{
          return("k=1 Unknown class of best.X")
        }
      } else{
        output[k,4:5]<-c(root$Ey.opt1,root$trt.opt1)
        break
      }
    } 
    else{
      for(j in (2^(k-1)):(2^k-1)){
        if(!is.na(output[trunc(j/2),2])){
          best.H.j<-best.Hsub(Hsub=Hsub[I.node==j,],A=A[I.node==j],mus.hat=mus.hat[I.node==j,],minsplit=minsplit)
          if(is.null(best.H.j)==F && best.H.j$mEy.opt1>as.numeric(output[trunc(j/2),4])+lambda){
            
            if(is.numeric(Hsub[,as.numeric(paste(best.H.j$X))])==T){ 
              output[j,-1]<-c(best.H.j$X, best.H.j$X.subset, best.H.j$mEy.opt1, NA)
              I.node[I.node==j & Hsub[,best.H.j$X] <= best.H.j$X.subset]<-2*j
              output[2*j,-1] <- c(NA,NA,NA,best.H.j$trt.L) 
              I.node[I.node==j & Hsub[,best.H.j$X] > best.H.j$X.subset]<-2*j+1
              output[2*j+1,-1]<-c(NA,NA,NA,best.H.j$trt.R) 
            } 
            else if(is.numeric(Hsub[,as.numeric(paste(best.H.j$X))])==F){ 
              output[j,2] <- round(as.numeric(paste(best.H.j$X)),0)
              output[j,3] <- paste(best.H.j$X.subset)
              output[j,4] <- as.numeric(best.H.j$mEy.opt1)
              output[j,5] <- NA
              I.node[I.node==j & Hsub[,as.numeric(paste(best.H.j$X))] %in% unlist(strsplit(paste(best.H.j$X.subset),split= " "))]<-2*j
              output[2*j,-1]<-c(NA,NA,NA,best.H.j$trt.L)
              I.node[I.node==j & !(Hsub[,as.numeric(paste(best.H.j$X))] %in% unlist(strsplit(paste(best.H.j$X.subset),split= " ")))]<-2*j+1
              output[2*j+1,-1]<-c(NA,NA,NA,best.H.j$trt.R)
            } else {
              return("error in k>=2")
            }
            
          }
        }  
      }
      if(sum(is.na(output[(2^(k-1)):(2^k-1),2]))==2^(k-1)) break ## means that you continue with the next k value
    }
  }
  output<-output[!is.na(output[,2]) | !is.na(output[,5]),]
  return(output)
}




