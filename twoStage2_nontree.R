## source functions
source("general_functions.R") ## general functions used in simulations
source("restl_functions.R") ## functions used to implement ReST-L
source("trl_functions.R") ## functions used to implement T-RL

library(rpart)
library(randomForest)

########################################################
# simulation - two stages - Section 4.1
########################################################

twoStage2_nontree <- function(N = 350, N2 = 1000, iter = 20, rho = 0.2, No.Var.H = 20, No.Var.Hsub = 7, Method = "REST-L", prop.model = "correct"){

select1<-select2<-selects<-rep(NA,iter) # percent of optimality
EYs<-rep(NA,iter) # estimated mean counterfactual outcome

for(i in 1:iter){
  set.seed(i)
  CorrMat <- Create.Cor.Mat5(No.Var.H, No.Var.Hsub, rho) ## K
  Z <- rbinom(N, 1, 0.4)
  X <- matrix((MASS::mvrnorm(n=N, mu=rep(0,No.Var.H), Sigma=CorrMat)),nrow=N) ## 20 normals
  if(No.Var.H == 20){
    x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
    x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,Z)
  } else if(No.Var.H == 50){
    x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
    x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
    x21<-X[,21]; x22<-X[,22]; x23<-X[,23]; x24<-X[,24]; x25<-X[,25]; x26<-X[,26]; x27<-X[,27]; x28<-X[,28]; x29<-X[,29]; x30<-X[,30]; 
    x31<-X[,31]; x32<-X[,32]; x33<-X[,33]; x34<-X[,34]; x35<-X[,35]; x36<-X[,36]; x37<-X[,37]; x38<-X[,38]; x39<-X[,39]; x40<-X[,40]; 
    x41<-X[,41]; x42<-X[,42]; x43<-X[,43]; x44<-X[,44]; x45<-X[,45]; x46<-X[,46]; x47<-X[,47]; x48<-X[,48]; x49<-X[,49]; x50<-X[,50]
    X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,Z)
  } else if(No.Var.H == 100){
    x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
    x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
    x21<-X[,21]; x22<-X[,22]; x23<-X[,23]; x24<-X[,24]; x25<-X[,25]; x26<-X[,26]; x27<-X[,27]; x28<-X[,28]; x29<-X[,29]; x30<-X[,30]; 
    x31<-X[,31]; x32<-X[,32]; x33<-X[,33]; x34<-X[,34]; x35<-X[,35]; x36<-X[,36]; x37<-X[,37]; x38<-X[,38]; x39<-X[,39]; x40<-X[,40]; 
    x41<-X[,41]; x42<-X[,42]; x43<-X[,43]; x44<-X[,44]; x45<-X[,45]; x46<-X[,46]; x47<-X[,47]; x48<-X[,48]; x49<-X[,49]; x50<-X[,50]; 
    x51<-X[,51]; x52<-X[,52]; x53<-X[,53]; x54<-X[,54]; x55<-X[,55]; x56<-X[,56]; x57<-X[,57]; x58<- X[,58]; x59<-X[,59]; x60<-X[,60]; 
    x61<-X[,61]; x62<-X[,62]; x63<-X[,63]; x64<-X[,64]; x65<-X[,65]; x66<-X[,66]; x67<-X[,67]; x68<- X[,68]; x69<-X[,69]; x70<-X[,70]; 
    x71<-X[,71]; x72<-X[,72]; x73<-X[,73]; x74<-X[,74]; x75<-X[,75]; x76<-X[,76]; x77<-X[,77]; x78<- X[,78]; x79<-X[,79]; x80<-X[,80]; 
    x81<-X[,81]; x82<-X[,82]; x83<-X[,83]; x84<-X[,84]; x85<-X[,85]; x86<-X[,86]; x87<-X[,87]; x88<- X[,88]; x89<-X[,89]; x90<-X[,90]; 
    x91<-X[,91]; x92<-X[,92]; x93<-X[,93]; x94<-X[,94]; x95<-X[,95]; x96<-X[,96]; x97<-X[,97]; x98<- X[,98]; x99<-X[,99]; x100<-X[,100]; 
    X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,
              x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100,Z)
  }
  
  ############### stage 1 data simulation ##############

    id1 <- No.Var.Hsub + 1
    id2 <- No.Var.Hsub + 2

    pi10<-rep(1,N); 
    pi11<-exp(0.5*X0[,id1]-0.5*x1 + 1*Z-0.5); 
    pi12<-exp(0.5*X0[,id2]+0.5*x1 - 1*Z)
    matrix.pi1<-cbind(pi10,pi11,pi12)
    
    A1<-A.sim(matrix.pi1)
    class.A1<-sort(unique(A1))
    
    # propensity stage 1
    if(prop.model == "correct"){
    pis1.hat<-M.propen(A1,cbind(x1,X0[,id1],X0[,id2],Z))
    } else if(prop.model == "incorrect"){
      pis1.hat<-M.propen(A1,cbind(X0))  
    } else(return("no propensity specification Stage 1"))
    
    # simulate stage 1 optimal g1.opt
    g1.opt <- (log2(abs(x1) + 1)  <= 2 )*((x2 < -0.25)) + (x2^2 > 0.35 ) ## about 65%
    # stage 1 outcome
    R1 <- exp(1.5+0.3*X0[,id1]-1.5*Z-abs(1.5*x1-2)*(A1-g1.opt)^2) + rnorm(N,0,1)

    ############### stage 2 data simulation ##############
    pi20<-rep(1.5,N); 
    pi21<-exp(0.2*R1+0.5-1*Z); 
    pi22<-exp(0.5*X0[,id2]+1*Z)
    matrix.pi2<-cbind(pi20,pi21,pi22)
    A2<-A.sim(matrix.pi2)
    class.A2<-sort(unique(A2))
    
    # propensity stage 2
    if(prop.model == "correct"){
      pis2.hat<-M.propen(A2,cbind(R1,X0[,id2],Z))
    } else if(prop.model == "incorrect"){
      pis2.hat<-M.propen(A2,cbind(R1,A1,X0))  
    } else(return("no propensity specification Stage 2"))
    
    # optimal g2.opt
    g2.opt <- (abs(x3)> 0.60)*(R1>1) + (R1^2>3)
    # stage 2 outcome R2
    R2<-exp(1.18+0.2*X0[,id2]-2*Z-abs(1.5*x3+2)*(A2-g2.opt)^2)+rnorm(N,0,1)

    ########### Using Backward induction #############################
    ############### stage 2 Estimation ###############################
  
if(Method == "REST-L"){
  # conditional mean model using linear regression 
  REG2<-Reg.mu.modified.Hsub(Y=R2,As=cbind(A1,A2),Hsub=cbind(X0[,1:No.Var.Hsub],R1), HsubC=cbind(X0[,(No.Var.Hsub+1):(No.Var.H+1)]))
  mus2.reg<-REG2$mus.reg

  tree2<-DTRtree.modified.HSUB(R2,A2,Hsub=cbind(X0[,1:No.Var.Hsub],A1,R1),HsubC=cbind(X0[,(No.Var.Hsub+1):(No.Var.H+1)]),pis.hat=pis2.hat,mus.reg=mus2.reg,lambda.pct=0.05,minsplit=max(0.05*N,20))

  ############### stage 1 Estimation ################################
  # calculate pseudo outcome (PO)

  # expected optimal stage 2 outcome 
  E.R2.tree<-rep(NA,N)
  
  # estimated optimal regime
  g2.tree<-predict.DTR.modified(tree2,newdata=data.frame(X0[,1:No.Var.Hsub],A1,R1))
  
  data <- data.frame(A2=A2,X0[,c(1:No.Var.Hsub)],A1=A1,R1=R1)
  RF2<-randomForest(R2~., data=data)
  mus2.RF<-matrix(NA,N,length(class.A2))
  for(d in 1L:length(class.A2)) mus2.RF[,d] <- predict(RF2,newdata=data.frame(A2=rep(class.A2[d],N),X0[,c(1:No.Var.Hsub)],A1=A1,R1=R1))
  
  for(m in 1:N){
    E.R2.tree[m]<-R2[m] + mus2.RF[m,g2.tree[m]+1]-mus2.RF[m,A2[m]+1] ## if the optimal is the same as the actual, E.R2.tree is just R2 ## removed the 
  }
  
  # pseudo outcomes
  PO.tree<-R1+E.R2.tree

  REG1<-Reg.mu.modified.Hsub(Y=PO.tree,As=A1,Hsub=cbind(X0[,1:No.Var.Hsub]), HsubC=cbind(X0[,(No.Var.Hsub+1):(No.Var.H+1)])) #
  mus1.reg<-REG1$mus.reg
  
  tree1<-DTRtree.modified.HSUB(PO.tree,A1,Hsub=cbind(X0[,1:No.Var.Hsub]),HsubC=cbind(X0[,(No.Var.Hsub+1):(No.Var.H+1)]),pis.hat=pis1.hat,mus.reg=mus1.reg, lambda.pct=0.05,minsplit=max(0.05*N,20)) ## 

} else if(Method == "T-RL"){ 

  REG2<-Reg.mu.modified(Y=R2,As=cbind(A1,A2),H=cbind(X0,R1))
  mus2.reg<-REG2$mus.reg
  
  tree2<-DTRtree.modified(R2,A2,H=cbind(X0,A1,R1),pis.hat=pis2.hat,mus.reg=mus2.reg,lambda.pct=0.05,minsplit=max(0.05*N,20))
  
  ############### stage 1 Estimation ################################
  # calculate pseudo outcome (PO)
  
  # expected optimal stage 2 outcome 
  E.R2.tree<-rep(NA,N)
  
  # estimated optimal regime
  g2.tree<-predict.DTR.modified(tree2,newdata=data.frame(X0,A1,R1))
  
  # random forest for the estimated mean
  A2 <- as.matrix(A2)
  RF2<-randomForest(R2~., data=data.frame(A2,X0,A1,R1))
  mus2.RF<-matrix(NA,N,length(class.A2))
  for(d in 1L:length(class.A2)) mus2.RF[,d]<-predict(RF2,newdata=data.frame(A2=rep(class.A2[d],N),X0,A1,R1))
  
  for(m in 1:N){
    E.R2.tree[m]<-R2[m] + mus2.RF[m,g2.tree[m]+1]-mus2.RF[m,A2[m]+1] ## if the optimal is the same as the actual, E.R2.tree is just R2 ## removed the 
  }
  
  # pseudo outcomes
  PO.tree<-R1+E.R2.tree

  tree1<-DTRtree.modified(PO.tree,A1,H=X0,pis.hat=pis1.hat,lambda.pct=0.05,minsplit=max(0.05*N,20))
  
} else if(Method == "Naive-T-RL"){ 

  REG2<-Reg.mu.modified(Y=R2,As=cbind(A1,A2),H=cbind(X0[,1:No.Var.Hsub],R1))
  mus2.reg<-REG2$mus.reg
  
  if(prop.model == "correct"){
    pis2.hat<-M.propen(A2,cbind(R1))
  } else if(prop.model == "incorrect"){
    pis2.hat<-M.propen(A2,cbind(R1,A1,X0[,1:No.Var.Hsub]))  
  } else(return("no propensity specification Stage 2"))
  
  tree2<-DTRtree.modified(R2,A2,H=cbind(X0[,1:No.Var.Hsub],A1,R1),pis.hat=pis2.hat,mus.reg=mus2.reg,lambda.pct=0.05,minsplit=max(0.05*N,20))
  
  ############### stage 1 Estimation ################################
  # calculate pseudo outcome (PO)
  
  # expected optimal stage 2 outcome 
  E.R2.tree<-rep(NA,N)
  
  # estimated optimal regime
  g2.tree<-predict.DTR.modified(tree2,newdata=data.frame(X0[,1:No.Var.Hsub],A1,R1))
  
  # random forest for the estimated mean
  A2 <- as.matrix(A2)
  RF2<-randomForest(R2~., data=data.frame(A2,X0[,1:No.Var.Hsub],A1,R1))
  mus2.RF<-matrix(NA,N,length(class.A2))
  for(d in 1L:length(class.A2)) mus2.RF[,d]<-predict(RF2,newdata=data.frame(A2=rep(class.A2[d],N),X0[,1:No.Var.Hsub],A1,R1))
  
  for(m in 1:N){
    E.R2.tree[m]<-R2[m] + mus2.RF[m,g2.tree[m]+1]-mus2.RF[m,A2[m]+1] ## if the optimal is the same as the actual, E.R2.tree is just R2 ## removed the 
  }
  
  # pseudo outcomes
  PO.tree<-R1+E.R2.tree
  
  # propensity stage 1
  if(prop.model == "correct"){
    pis1.hat<-M.propen(A1,cbind(x1))
  } else if(prop.model == "incorrect"){
    pis1.hat<-M.propen(A1,cbind(X0[,1:No.Var.Hsub]))  
  } else(return("no propensity specification Stage 1"))
  
  tree1<-DTRtree.modified(PO.tree,A1,H=X0[,1:No.Var.Hsub],pis.hat=pis1.hat,lambda.pct=0.05,minsplit=max(0.05*N,20))
  
} else if(Method == "Q-Linear"){
  
  A2 <- as.factor(A2)
  QLinear2 <- lm(R2 ~ A2*., data = data.frame(X0,A1,R1,A2,R2))

  ## predict
  newdata <- data.frame() 
  E.R2.tree<-rep(NA,N)

  # estimated optimal regime
  A2.vec <- rep(0,N)
  A2.vec1 <- rep(1,N)
  A2.vec2 <- rep(2,N)
  newdata1=data.frame(X0,A1,R1,A2=A2.vec)
  newdata2 <- data.frame(X0,A1,R1,A2=A2.vec1)
  newdata3 <- data.frame(X0,A1,R1,A2=A2.vec2)
  newdata <- rbind(newdata1,newdata2,newdata3)
  newdata$A2 = as.factor(newdata$A2)
  g2.tree<-matrix(predict(QLinear2,newdata=newdata), nrow=N, ncol=3)
  max.Y2.pred <- E.R2.tree <- apply(g2.tree, 1, max)
  A2.hat.opt <- ifelse(g2.tree[,1] > g2.tree[,2] & g2.tree[,1] > g2.tree[,3], 0,
                       ifelse(g2.tree[,2] > g2.tree[,1] & g2.tree[,2] > g2.tree[,3], 1,
                              ifelse(g2.tree[,3] > g2.tree[,1] & g2.tree[,3] > g2.tree[,2], 2,NA)))
    data2b <- data.frame(X = rbind(newdata), A2.opt = rep(A2.hat.opt, times=3), max.Y2.pred = rep(max.Y2.pred, times=3), A2.actual = rep(A2, times=3),R2 = rep(R2, times=3))
    colnames(data2b) <- c(colnames(newdata), "A2.opt", "max.Y2.pred","A2.actual", "R2")
    data2b$R2.tilde <- ifelse(data2b$A2.actual == data2b$A2.opt, data2b$R2, data2b$max.Y2.pred)
    data2b$R1.tilde <- data2b$R2.tilde + data2b$R1

    data3 <- data2b[data2b$A2 == data2b$A2.opt,]
    data3$A1 <- as.factor(data3$A1)
    QLinear1 <- lm(R1.tilde ~ A1*., data = data3[,c(1:(No.Var.H+2),No.Var.H+10)]) 
  
} else if(Method == "Q-Linear-REST"){
  
  A2 <- as.factor(A2)
  if(No.Var.Hsub == 7 & No.Var.H == 20){
    formula <- R2 ~ A2*(x1+x2+x3+x4+x5+x6+x7+A1+R1) + x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+Z
  } else if(No.Var.Hsub == 10 & No.Var.H == 50){
    formula <- R2 ~ A2*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+A1+R1) + x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+Z
  } else if(No.Var.Hsub == 20 & No.Var.H == 100){
    formula <- R2 ~ A2*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+A1+R1) + x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+
      x51+x52+x53+x54+x55+x56+x57+x58+x59+x60+x61+x62+x63+x64+x65+x66+x67+x68+x69+x70+x71+x72+x73+x74+x75+x76+x77+x78+x79+x80+x81+x82+x83+x84+x85+x86+x87+x88+x89+x90+x91+x92+x93+x94+x95+x96+x97+x98+x99+x100+Z
  } else if(No.Var.Hsub == 35 & No.Var.H == 50){
    formula <- R2 ~ A2*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+A1+R1)+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+Z
  } 
  QLinear2 <- lm(formula)
  
  ## predict
  newdata <- data.frame() 
  E.R2.tree<-rep(NA,N)
  
  # estimated optimal regime
  A2.vec <- rep(0,N)
  A2.vec1 <- rep(1,N)
  A2.vec2 <- rep(2,N)
  newdata1 <- data.frame(X0,A1,R1,A2=A2.vec)
  newdata2 <- data.frame(X0,A1,R1,A2=A2.vec1)
  newdata3 <- data.frame(X0,A1,R1,A2=A2.vec2)
  newdata <- rbind(newdata1,newdata2,newdata3)
  newdata$A2 = as.factor(newdata$A2)
  colnames(newdata) <- c(colnames(X0),"A1","R1","A2")
  g2.tree<-matrix(predict(QLinear2,newdata=cbind(newdata)), nrow=N, ncol=3)
  max.Y2.pred <- E.R2.tree <- apply(g2.tree, 1, max)
  A2.hat.opt <- ifelse(g2.tree[,1] > g2.tree[,2] & g2.tree[,1] > g2.tree[,3], 0, 
                       ifelse(g2.tree[,2] > g2.tree[,1] & g2.tree[,2] > g2.tree[,3], 1,
                              ifelse(g2.tree[,3] > g2.tree[,1] & g2.tree[,3] > g2.tree[,2], 2,NA)))  
  data2b <- data.frame(X = rbind(newdata), A2.opt = rep(A2.hat.opt, times=3), max.Y2.pred = rep(max.Y2.pred, times=3), A2.actual = rep(A2, times=3),R2 = rep(R2, times=3))
  colnames(data2b) <- c(colnames(newdata), "A2.opt", "max.Y2.pred","A2.actual", "R2")
  data2b$R2.tilde <- ifelse(data2b$A2.actual == data2b$A2.opt, data2b$R2, data2b$max.Y2.pred)
  data2b$R1.tilde <- data2b$R2.tilde + data2b$R1
  
  data3 <- data2b[data2b$A2 == data2b$A2.opt,]

  data3$A1 <- as.factor(data3$A1)
  if(No.Var.Hsub == 7 & No.Var.H == 20){
    formula2 <- R1.tilde ~ A1*(x1+x2+x3+x4+x5+x6+x7) + x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+Z
  } else if(No.Var.Hsub == 10 & No.Var.H == 50){
    formula2 <- R1.tilde ~ A1*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10) + x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+Z
  } else if(No.Var.Hsub == 20 & No.Var.H == 100){
    formula2 <- R1.tilde ~ A1*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20) + x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+
      x51+x52+x53+x54+x55+x56+x57+x58+x59+x60+x61+x62+x63+x64+x65+x66+x67+x68+x69+x70+x71+x72+x73+x74+x75+x76+x77+x78+x79+x80+x81+x82+x83+x84+x85+x86+x87+x88+x89+x90+x91+x92+x93+x94+x95+x96+x97+x98+x99+x100+Z
  } else if(No.Var.Hsub == 35 & No.Var.H == 50){
    formula2 <- R1.tilde ~ A1*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35)+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+Z
  } 
  QLinear1 <- lm(formula2, data = data3[,c(1:(No.Var.H+2),No.Var.H+10)]) 
  
} else if(Method == "Q-NP"){
  
  QRF2 <- randomForest(R2 ~ ., data = data.frame(X0,A1,R1,A2,R2))

  ## predict
  newdata <- data.frame() 
  E.R2.tree<-rep(NA,N)

  # estimated optimal regime
  A2.vec <- rep(0,N)
  A2.vec1 <- rep(1,N)
  A2.vec2 <- rep(2,N)
  newdata1=data.frame(X0,A1,R1,A2=A2.vec)
  newdata2 <- data.frame(X0,A1,R1,A2=A2.vec1)
  newdata3 <- data.frame(X0,A1,R1,A2=A2.vec2)
  newdata <- rbind(newdata1,newdata2,newdata3)
  g2.tree<-matrix(predict(QRF2,newdata=newdata), nrow=N, ncol=3)
  max.Y2.pred <- E.R2.tree <- apply(g2.tree, 1, max)
  A2.hat.opt <- ifelse(g2.tree[,1] > g2.tree[,2] & g2.tree[,1] > g2.tree[,3], 0,
                       ifelse(g2.tree[,2] > g2.tree[,1] & g2.tree[,2] > g2.tree[,3], 1,
                              ifelse(g2.tree[,3] > g2.tree[,1] & g2.tree[,3] > g2.tree[,2], 2,NA)))
  data2b <- data.frame(X = rbind(newdata), A2.opt = rep(A2.hat.opt, times=3), max.Y2.pred = rep(max.Y2.pred, times=3), A2.actual = rep(A2, times=3),R2 = rep(R2, times=3))
  colnames(data2b) <- c(colnames(newdata), "A2.opt", "max.Y2.pred","A2.actual", "R2")
  data2b$R2.tilde <- ifelse(data2b$A2.actual == data2b$A2.opt, data2b$R2, data2b$max.Y2.pred)
  data2b$R1.tilde <- data2b$R2.tilde + data2b$R1

  data3 <- data2b[data2b$A2 == data2b$A2.opt,]

  QRF1 <- randomForest(R1.tilde ~ ., data = data3[,c(1:(No.Var.H+2),No.Var.H+10)], na.action = na.omit) 
  
} else if(Method == "Q-NP-REST"){
  
  if(No.Var.Hsub == 7 ){
    formula3 <- R2 ~ A2+x1+x2+x3+x4+x5+x6+x7+A1+R1
  } else if(No.Var.Hsub == 10 ){
    formula3 <- R2 ~ A2+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+A1+R1
  } else if(No.Var.Hsub == 20 ){
    formula3 <- R2 ~ A2+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+A1+R1
  } else if(No.Var.Hsub == 35 ){
    formula3 <- R2 ~ A2+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+A1+R1
  }
  QRF2 <- randomForest(formula3)
  
  ## predict
  newdata <- data.frame() 
  E.R2.tree<-rep(NA,N)
  
  # estimated optimal regime
  A2.vec <- rep(0,N)
  A2.vec1 <- rep(1,N)
  A2.vec2 <- rep(2,N)
  newdata1 <- data.frame(X0,A1,R1,A2=A2.vec)
  newdata2 <- data.frame(X0,A1,R1,A2=A2.vec1)
  newdata3 <- data.frame(X0,A1,R1,A2=A2.vec2)
  newdata <- rbind(newdata1,newdata2,newdata3)
  g2.tree<-matrix(predict(QRF2,newdata=cbind(newdata)), nrow=N, ncol=3)
  max.Y2.pred <- E.R2.tree <- apply(g2.tree, 1, max)
  A2.hat.opt <- ifelse(g2.tree[,1] > g2.tree[,2] & g2.tree[,1] > g2.tree[,3], 0, 
                       ifelse(g2.tree[,2] > g2.tree[,1] & g2.tree[,2] > g2.tree[,3], 1,
                              ifelse(g2.tree[,3] > g2.tree[,1] & g2.tree[,3] > g2.tree[,2], 2,NA)))  
  data2b <- data.frame(X = rbind(newdata), A2.opt = rep(A2.hat.opt, times=3), max.Y2.pred = rep(max.Y2.pred, times=3), A2.actual = rep(A2, times=3),R2 = rep(R2, times=3))
  colnames(data2b) <- c(colnames(newdata), "A2.opt", "max.Y2.pred","A2.actual", "R2")
  data2b$R2.tilde <- ifelse(data2b$A2.actual == data2b$A2.opt, data2b$R2, data2b$max.Y2.pred)
  data2b$R1.tilde <- data2b$R2.tilde + data2b$R1
  
  data3 <- data2b[data2b$A2 == data2b$A2.opt,]

  if(No.Var.Hsub == 7 ){
    formula4 <- R1.tilde ~ A1+x1+x2+x3+x4+x5+x6+x7
  } else if(No.Var.Hsub == 10 ){
    formula4 <- R1.tilde ~ A1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10
  } else if(No.Var.Hsub == 20 ){
    formula4 <- R1.tilde ~ A1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20
  } else if(No.Var.Hsub == 35 ){
    formula4 <- R1.tilde ~ A1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35
  }
  QRF1 <- randomForest(formula4, data = data3[,c(1:(No.Var.H+2),No.Var.H+10)], na.action = na.omit)
  
} else{return("error-no method specified")}

############################################
# prediction using new data
############################################
set.seed(i+10000)
X <- matrix((MASS::mvrnorm(n=N2, mu=rep(0,No.Var.H), Sigma=CorrMat)),nrow=N2) ## 20 normals
Z <- rbinom(N2, 1, 0.4)
if(No.Var.H == 20){
  x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
  x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,Z)
} else if(No.Var.H == 50){
  x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
  x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  x21<-X[,21]; x22<-X[,22]; x23<-X[,23]; x24<-X[,24]; x25<-X[,25]; x26<-X[,26]; x27<-X[,27]; x28<-X[,28]; x29<-X[,29]; x30<-X[,30]; 
  x31<-X[,31]; x32<-X[,32]; x33<-X[,33]; x34<-X[,34]; x35<-X[,35]; x36<-X[,36]; x37<-X[,37]; x38<-X[,38]; x39<-X[,39]; x40<-X[,40]; 
  x41<-X[,41]; x42<-X[,42]; x43<-X[,43]; x44<-X[,44]; x45<-X[,45]; x46<-X[,46]; x47<-X[,47]; x48<-X[,48]; x49<-X[,49]; x50<-X[,50]
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,Z)
} else if(No.Var.H == 100){
  x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
  x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  x21<-X[,21]; x22<-X[,22]; x23<-X[,23]; x24<-X[,24]; x25<-X[,25]; x26<-X[,26]; x27<-X[,27]; x28<-X[,28]; x29<-X[,29]; x30<-X[,30]; 
  x31<-X[,31]; x32<-X[,32]; x33<-X[,33]; x34<-X[,34]; x35<-X[,35]; x36<-X[,36]; x37<-X[,37]; x38<-X[,38]; x39<-X[,39]; x40<-X[,40]; 
  x41<-X[,41]; x42<-X[,42]; x43<-X[,43]; x44<-X[,44]; x45<-X[,45]; x46<-X[,46]; x47<-X[,47]; x48<-X[,48]; x49<-X[,49]; x50<-X[,50]; 
  x51<-X[,51]; x52<-X[,52]; x53<-X[,53]; x54<-X[,54]; x55<-X[,55]; x56<-X[,56]; x57<-X[,57]; x58<- X[,58]; x59<-X[,59]; x60<-X[,60]; 
  x61<-X[,61]; x62<-X[,62]; x63<-X[,63]; x64<-X[,64]; x65<-X[,65]; x66<-X[,66]; x67<-X[,67]; x68<- X[,68]; x69<-X[,69]; x70<-X[,70]; 
  x71<-X[,71]; x72<-X[,72]; x73<-X[,73]; x74<-X[,74]; x75<-X[,75]; x76<-X[,76]; x77<-X[,77]; x78<- X[,78]; x79<-X[,79]; x80<-X[,80]; 
  x81<-X[,81]; x82<-X[,82]; x83<-X[,83]; x84<-X[,84]; x85<-X[,85]; x86<-X[,86]; x87<-X[,87]; x88<- X[,88]; x89<-X[,89]; x90<-X[,90]; 
  x91<-X[,91]; x92<-X[,92]; x93<-X[,93]; x94<-X[,94]; x95<-X[,95]; x96<-X[,96]; x97<-X[,97]; x98<- X[,98]; x99<-X[,99]; x100<-X[,100]; 
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,
            x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100,Z)
}

# true optimal for regime at stage 1
  g1.opt <- (log2(abs(x1) + 1)  <= 2 )*((x2 < -0.25)) + (x2^2 > 0.35 ) ## about 65%

  R1<-exp(1.5+0.3*X0[,id1]-1.5*Z)+rnorm(N2,0,1)
  
  R2<-exp(1.18+0.2*X0[,id2]-2*Z)+rnorm(N2,0,1)

  
####### stage 1 prediction #######

  # predict selection %

  if(Method == "REST-L"){

  g1.tree<-predict.DTR.modified2(tree1,newdata=data.frame(X0[,1:No.Var.Hsub])) 
  
  select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
  
  R1.tree<-exp(1.5+0.3*X0[,id1]-1.5*Z-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
  

####### stage 2 prediction #######

  g2.tree<-predict.DTR.modified2(tree2,newdata=data.frame(X0[,c(1:No.Var.Hsub)],A1=g1.tree,R1=R1.tree)) ## modified to include only values in Hsub instead of all X0
 # true optimal for regime at stage 2
  g2.opt.tree <- (abs(x3)> 0.60)*(R1.tree>1) + (R1.tree^2>3)

  select2[i]<-mean(g2.tree==g2.opt.tree, na.rm=TRUE)

  selects[i]<-mean(g1.tree==g1.opt & g2.tree==g2.opt.tree, na.rm=TRUE)

# predict R2
  R2.tree<-exp(1.18+0.2*X0[,id2]-2*Z-abs(1.5*x3+2)*(g2.tree-g2.opt.tree)^2)+rnorm(N2,0,1)
  
  
  EYs[i]<-mean(R1.tree+R2.tree, na.rm=TRUE)
    
  if(i%%50==0) print(i)
  } else if(Method == "T-RL"){
    
    # predict selection %
    
    g1.tree<-predict.DTR.modified2(tree1,newdata=data.frame(X0))
    
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)

    R1.tree<-exp(1.5+0.3*X0[,id1]-1.5*Z-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
  
    
    ####### stage 2 prediction #######
    g2.tree<-predict.DTR.modified2(tree2,newdata=data.frame(X0,A1=g1.tree,R1=R1.tree))
    # true optimal for regime at stage 2
    g2.opt.tree <- (abs(x3)> 0.60)*(R1.tree>1) + (R1.tree^2>3)
    
    select2[i]<-mean(g2.tree==g2.opt.tree, na.rm=TRUE)
    
    selects[i]<-mean(g1.tree==g1.opt & g2.tree==g2.opt.tree, na.rm=TRUE)
    
    # predict R2
    R2.tree<-exp(1.18+0.2*X0[,id2]-2*Z-abs(1.5*x3+2)*(g2.tree-g2.opt.tree)^2)+rnorm(N2,0,1)
  
    EYs[i]<-mean(R1.tree+R2.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
  } else if(Method == "Naive-T-RL"){
    
    # predict selection %
    
    g1.tree<-predict.DTR.modified2(tree1,newdata=data.frame(X0[,1:No.Var.Hsub]))
    
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
    
    R1.tree<-exp(1.5+0.3*X0[,id1]-1.5*Z-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
    
    ####### stage 2 prediction #######
    
    g2.tree<-predict.DTR.modified2(tree2,newdata=data.frame(X0[,1:No.Var.Hsub],A1=g1.tree,R1=R1.tree))
    # true optimal for regime at stage 2
    g2.opt.tree <- (abs(x3)> 0.60)*(R1.tree>1) + (R1.tree^2>3)
    
    select2[i]<-mean(g2.tree==g2.opt.tree, na.rm=TRUE)
    selects[i]<-mean(g1.tree==g1.opt & g2.tree==g2.opt.tree, na.rm=TRUE)
    
    # predict R2
    R2.tree<-exp(1.18+0.2*X0[,id2]-2*Z-abs(1.5*x3+2)*(g2.tree-g2.opt.tree)^2)+rnorm(N2,0,1)
    
    EYs[i]<-mean(R1.tree+R2.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
  } else if(Method == "Q-Linear" | Method == "Q-Linear-REST"){
    
    ### predict based on QLinear2 and QLinear1
    newdata <- data.frame() 
        # estimated optimal regime
    A1.vec <- rep(0,N2)
    A1.vec1 <- rep(1,N2)
    A1.vec2 <- rep(2,N2)
    newdata1=data.frame(X0,R1, g1.opt,A1=A1.vec)
    newdata2 <- data.frame(X0,R1, g1.opt,A1=A1.vec1)
    newdata3 <- data.frame(X0,R1, g1.opt,A1=A1.vec2)
    newdata <- rbind(newdata1,newdata2,newdata3)
    newdata$A1 <- as.factor(newdata$A1)
    g1.tree2<-matrix(predict(QLinear1,newdata=newdata), nrow=N2, ncol=3)
    max.Y1.pred <- E.R1.tree <- apply(g1.tree2, 1, max)
    g1.tree <- ifelse(g1.tree2[,1] > g1.tree2[,2] & g1.tree2[,1] > g1.tree2[,3], 0, 
                         ifelse(g1.tree2[,2] > g1.tree2[,1] & g1.tree2[,2] > g1.tree2[,3], 1,
                                ifelse(g1.tree2[,3] > g1.tree2[,1] & g1.tree2[,3] > g1.tree2[,2], 2,NA)))  
    
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
    
    R1.tree<-exp(1.5+0.3*X0[,id1]-1.5*Z-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
  
    ####### stage 2 prediction #######
    newdata <- data.frame() 
    # estimated optimal regime
    A2.vec <- rep(0,N2)
    A2.vec1 <- rep(1,N2)
    A2.vec2 <- rep(2,N2)
    newdata1=data.frame(X0,R1=R1.tree, A1 = g1.tree, A2=A2.vec)
    newdata2 <- data.frame(X0,R1=R1.tree, A1= g1.tree,A2=A2.vec1)
    newdata3 <- data.frame(X0,R1=R1.tree, A1=g1.tree,A2=A2.vec2)
    newdata <- rbind(newdata1,newdata2,newdata3)
    newdata$A2 <- as.factor(newdata$A2)
    g2.tree2<-matrix(predict(QLinear2,newdata=newdata), nrow=N2, ncol=3)
    max.Y2.pred <- E.R2.tree <- apply(g2.tree2, 1, max)
    g2.tree <- ifelse(g2.tree2[,1] > g2.tree2[,2] & g2.tree2[,1] > g2.tree2[,3], 0, 
                      ifelse(g2.tree2[,2] > g2.tree2[,1] & g2.tree2[,2] > g2.tree2[,3], 1,
                             ifelse(g2.tree2[,3] > g2.tree2[,1] & g2.tree2[,3] > g2.tree2[,2], 2,NA)))  

    # true optimal for regime at stage 2
    g2.opt.tree <- (abs(x3)> 0.60)*(R1.tree>1) + (R1.tree^2>3)
    
    select2[i]<-mean(g2.tree==g2.opt.tree, na.rm=TRUE)
    
    selects[i]<-mean(g1.tree==g1.opt & g2.tree==g2.opt.tree, na.rm=TRUE)
    
    # predict R2
    R2.tree<-exp(1.18+0.2*X0[,id2]-2*Z-abs(1.5*x3+2)*(g2.tree-g2.opt.tree)^2)+rnorm(N2,0,1)

    EYs[i]<-mean(R1.tree+R2.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
    
} else if(Method == "Q-NP" | Method == "Q-NP-REST"){
    
    ### predict based on QLinear2 and QLinear1
    newdata <- data.frame() 
    # estimated optimal regime
    A1.vec <- rep(0,N2)
    A1.vec1 <- rep(1,N2)
    A1.vec2 <- rep(2,N2)
    newdata1=data.frame(X0,R1, g1.opt,A1=A1.vec)
    newdata2 <- data.frame(X0,R1, g1.opt,A1=A1.vec1)
    newdata3 <- data.frame(X0,R1, g1.opt,A1=A1.vec2)
    newdata <- rbind(newdata1,newdata2,newdata3)
    g1.tree2<-matrix(predict(QRF1,newdata=newdata), nrow=N2, ncol=3) ## if there were missing data in the estimation it will still predict g1.tree2 so this shouldn't be a problem
    max.Y1.pred <- E.R1.tree <- apply(g1.tree2, 1, max)
    g1.tree <- ifelse(g1.tree2[,1] > g1.tree2[,2] & g1.tree2[,1] > g1.tree2[,3], 0, 
                      ifelse(g1.tree2[,2] > g1.tree2[,1] & g1.tree2[,2] > g1.tree2[,3], 1,
                             ifelse(g1.tree2[,3] > g1.tree2[,1] & g1.tree2[,3] > g1.tree2[,2], 2,NA)))  
    
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)

    R1.tree<-exp(1.5+0.3*X0[,id1]-1.5*Z-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)

    ####### stage 2 prediction #######
    newdata <- data.frame() 
    # estimated optimal regime
    A2.vec <- rep(0,N2)
    A2.vec1 <- rep(1,N2)
    A2.vec2 <- rep(2,N2)
    newdata1 <- data.frame(X0,R1=R1.tree, A1 = g1.tree, A2=A2.vec)
    newdata2 <- data.frame(X0,R1=R1.tree, A1= g1.tree,A2=A2.vec1)
    newdata3 <- data.frame(X0,R1=R1.tree, A1=g1.tree,A2=A2.vec2)
    newdata <- rbind(newdata1,newdata2,newdata3)
    g2.tree2<-matrix(predict(QRF2,newdata=newdata), nrow=N2, ncol=3)
    max.Y2.pred <- E.R2.tree <- apply(g2.tree2, 1, max)
    g2.tree <- ifelse(g2.tree2[,1] > g2.tree2[,2] & g2.tree2[,1] > g2.tree2[,3], 0, 
                      ifelse(g2.tree2[,2] > g2.tree2[,1] & g2.tree2[,2] > g2.tree2[,3], 1,
                             ifelse(g2.tree2[,3] > g2.tree2[,1] & g2.tree2[,3] > g2.tree2[,2], 2,NA)))  
    
    # true optimal for regime at stage 2
    g2.opt.tree <- (abs(x3)> 0.60)*(R1.tree>1) + (R1.tree^2>3)
    
    select2[i]<-mean(g2.tree==g2.opt.tree, na.rm=TRUE)
    
    selects[i]<-mean(g1.tree==g1.opt & g2.tree==g2.opt.tree, na.rm=TRUE)
    
    # predict R2
    R2.tree<-exp(1.18+0.2*X0[,id2]-2*Z-abs(1.5*x3+2)*(g2.tree-g2.opt.tree)^2)+rnorm(N2,0,1)
    
    EYs[i]<-mean(R1.tree+R2.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
    
  } else(return("error-2method not specified"))
  
  
}
return(c(No.Stages = 2, No.RXperStage = 3, rho=rho, No.Var.H=No.Var.H, No.Var.Hsub=No.Var.Hsub, n.Sample = N, n.TestSet = N2, propensity.model = prop.model, true.DTR = "nontree", B = iter, Method = Method,
         MedianOpt = round(median(selects, na.rm=TRUE),3), IQROpt = round(IQR(selects, na.rm=TRUE),3), MedianEYs = round(median(EYs, na.rm=TRUE),3), IQREYs = round(IQR(EYs, na.rm=TRUE),3)))
} ## end function


