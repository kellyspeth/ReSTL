## source functions
source("general_functions.R") ## general functions used in simulations
source("restl_functions.R") ## functions used to implement ReST-L
source("trl_functions.R") ## functions used to implement T-RL

library(rpart)
library(randomForest)

########################################################
# simulation - one stage - Section 4.2
########################################################

oneStage_nontree <- function(N = 350, N2 = 1000, iter = 20, rho = 0.2, No.Var.H = 20, No.Var.Hsub = 7, Method = "REST-L", prop.model = "correct"){

select1<-rep(NA,iter) 
EYs<-rep(NA,iter) 

for(i in 1:iter){
  set.seed(i)
  CorrMat <- Create.Cor.Mat4(No.Var.H, No.Var.Hsub, rho) 
  X <- matrix((MASS::mvrnorm(n=N, mu=rep(0,No.Var.H), Sigma=CorrMat)),nrow=N) 
  if(No.Var.H == 20){
    x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
    x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)
  } else if(No.Var.H == 50){
    x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
    x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
    x21<-X[,21]; x22<-X[,22]; x23<-X[,23]; x24<-X[,24]; x25<-X[,25]; x26<-X[,26]; x27<-X[,27]; x28<-X[,28]; x29<-X[,29]; x30<-X[,30]; 
    x31<-X[,31]; x32<-X[,32]; x33<-X[,33]; x34<-X[,34]; x35<-X[,35]; x36<-X[,36]; x37<-X[,37]; x38<-X[,38]; x39<-X[,39]; x40<-X[,40]; 
    x41<-X[,41]; x42<-X[,42]; x43<-X[,43]; x44<-X[,44]; x45<-X[,45]; x46<-X[,46]; x47<-X[,47]; x48<-X[,48]; x49<-X[,49]; x50<-X[,50]
    X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50)
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
              x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100)
  }

    id1 <- No.Var.Hsub + 1
    id2 <- No.Var.Hsub + 2

pi10<-rep(1,N); 
pi11<-exp(0.5*X0[,id1]+0.5*x1); 
pi12<-exp(0.5*X0[,id2]-0.5*x1)
# weights matrix
matrix.pi1<-cbind(pi10,pi11,pi12)

A1<-A.sim(matrix.pi1)
class.A1<-sort(unique(A1))

if(prop.model == "correct"){
pis1.hat<-M.propen(A1,cbind(x1,X0[,id1],X0[,id2]))
} else if(prop.model == "incorrect"){
  pis1.hat<-M.propen(A1,cbind(X0))  
} else(return("no propensity specification Stage 1"))

# simulate stage 1 optimal g1.opt
g1.opt <- (log2(abs(x1) + 1)  <= 2 )*((x2 < 0.25)) + (x2^2 <= 0.5 ) ## about 65%
# stage 1 outcome
R1 <- exp(1.5+0.3*X0[,id1]-abs(1.5*x1-2)*(A1-g1.opt)^2) + rnorm(N,0,1)

############### stage 1 Estimation ###############################
###########################################
  
if(Method == "REST-L"){
  REG1<-Reg.mu.modified.Hsub(Y=R1,As=cbind(A1),Hsub=cbind(X0[,1:No.Var.Hsub]), HsubC=cbind(X0[,(No.Var.Hsub+1):No.Var.H]))
  mus1.reg<-REG1$mus.reg
  
  tree1<-DTRtree.modified.HSUB(R1,A1,Hsub=cbind(X0[,1:No.Var.Hsub]),HsubC=cbind(X0[,(No.Var.Hsub+1):No.Var.H]),pis.hat=pis1.hat,mus.reg=mus1.reg,lambda.pct=0.05,minsplit=max(0.05*N,20))

} else if(Method == "T-RL"){ 
  REG1<-Reg.mu.modified(Y=R1,As=cbind(A1),H=cbind(X0))
  mus1.reg<-REG1$mus.reg
  
  tree1<-DTRtree.modified(R1,A1,H=cbind(X0),pis.hat=pis1.hat,mus.reg=mus1.reg,lambda.pct=0.05,minsplit=max(0.05*N,20))
  
} else if(Method == "Q-Linear"){
  
  QLinear1 <- lm(R1 ~ A1*., data = data.frame(X0,A1,R1))
  
} else if(Method == "Q-Linear-REST"){
  
  if(No.Var.Hsub == 7 & No.Var.H == 20){
    formula <- R1 ~ A1*(x1+x2+x3+x4+x5+x6+x7) + x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20
  } else if(No.Var.Hsub == 10 & No.Var.H == 50){
    formula <- R1 ~ A1*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10) + x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50
  } else if(No.Var.Hsub == 20 & No.Var.H == 100){
    formula <- R1 ~ A1*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20) + x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50+
      x51+x52+x53+x54+x55+x56+x57+x58+x59+x60+x61+x62+x63+x64+x65+x66+x67+x68+x69+x70+x71+x72+x73+x74+x75+x76+x77+x78+x79+x80+x81+x82+x83+x84+x85+x86+x87+x88+x89+x90+x91+x92+x93+x94+x95+x96+x97+x98+x99+x100
  } else if(No.Var.Hsub == 35 & No.Var.H == 50){
    formula <- R1 ~ A1*(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35)+x36+x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+x47+x48+x49+x50
  }
  QLinear1 <- lm(formula)
  
  
} else if(Method == "Q-NP"){
  
  QRF1 <- randomForest(R1 ~ ., data = data.frame(X0,A1,R1))
  
} else if(Method == "Q-NP-REST"){
  
  if(No.Var.Hsub == 7 ){
    formula3 <- R1 ~ A1+x1+x2+x3+x4+x5+x6+x7
  } else if(No.Var.Hsub == 10 ){
    formula3 <- R1 ~ A1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10
  } else if(No.Var.Hsub == 20 ){
    formula3 <- R1 ~ A1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20
  } else if(No.Var.Hsub == 35 ){
    formula3 <- R1 ~ A1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35
  }
  QRF1 <- randomForest(formula3)
  
} else{return("error-no method specified")}

############################################
# prediction - first generate new data
############################################

set.seed(i+10000)
X <- matrix((MASS::mvrnorm(n=N2, mu=rep(0,No.Var.H), Sigma=CorrMat)),nrow=N2) ## 20 normals
if(No.Var.H == 20){
  x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
  x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)
} else if(No.Var.H == 50){
  x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; x4<-X[,4]; x5<-X[,5]; x6<-X[,6]; x7<-X[,7]; x8<- X[,8]; x9<-X[,9]; x10<-X[,10]; 
  x11<-X[,11]; x12<-X[,12]; x13<-X[,13]; x14<-X[,14]; x15<-X[,15]; x16<-X[,16]; x17<-X[,17]; x18<-X[,18]; x19<-X[,19]; x20<-X[,20];
  x21<-X[,21]; x22<-X[,22]; x23<-X[,23]; x24<-X[,24]; x25<-X[,25]; x26<-X[,26]; x27<-X[,27]; x28<-X[,28]; x29<-X[,29]; x30<-X[,30]; 
  x31<-X[,31]; x32<-X[,32]; x33<-X[,33]; x34<-X[,34]; x35<-X[,35]; x36<-X[,36]; x37<-X[,37]; x38<-X[,38]; x39<-X[,39]; x40<-X[,40]; 
  x41<-X[,41]; x42<-X[,42]; x43<-X[,43]; x44<-X[,44]; x45<-X[,45]; x46<-X[,46]; x47<-X[,47]; x48<-X[,48]; x49<-X[,49]; x50<-X[,50]
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,
            x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50)
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
  X0<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,
            x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,
            x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,
            x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100)
}

# true optimal for regime at stage 1
  g1.opt <- (log2(abs(x1) + 1)  <= 2 )*((x2 < 0.25)) + (x2^2 <= 0.5 ) ## 

  R1<-exp(1.5+0.3*X0[,id1])+rnorm(N2,0,1)

  
####### stage 1 prediction #######

  if(Method == "REST-L"){

    g1.tree<-predict.DTR.modified2(tree1,newdata=data.frame(X0[,1:No.Var.Hsub])) 
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
    R1.tree<-exp(1.5+0.3*X0[,id1]-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
    EYs[i]<-mean(R1.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
  } else if(Method == "T-RL"){
    
    g1.tree<-predict.DTR.modified2(tree1,newdata=data.frame(X0))
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
    R1.tree<-exp(1.5+0.3*X0[,id1]-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
    EYs[i]<-mean(R1.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
  } else if(Method == "Q-Linear" | Method == "Q-Linear-REST"){
    
    newdata <- data.frame() 
    A1.vec <- rep(0,N2)
    A1.vec1 <- rep(1,N2)
    A1.vec2 <- rep(2,N2)
    newdata1=data.frame(X0,R1, g1.opt,A1=A1.vec)
    newdata2 <- data.frame(X0,R1, g1.opt,A1=A1.vec1)
    newdata3 <- data.frame(X0,R1, g1.opt,A1=A1.vec2)
    newdata <- rbind(newdata1,newdata2,newdata3)
    g1.tree2<-matrix(predict(QLinear1,newdata=newdata), nrow=N2, ncol=3)
    max.Y1.pred <- E.R1.tree <- apply(g1.tree2, 1, max)
    g1.tree <- ifelse(g1.tree2[,1] > g1.tree2[,2] & g1.tree2[,1] > g1.tree2[,3], 0, 
                         ifelse(g1.tree2[,2] > g1.tree2[,1] & g1.tree2[,2] > g1.tree2[,3], 1,
                                ifelse(g1.tree2[,3] > g1.tree2[,1] & g1.tree2[,3] > g1.tree2[,2], 2,NA)))  
    
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
    R1.tree<-exp(1.5+0.3*X0[,id1]-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
    EYs[i]<-mean(R1.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
    
} else if(Method == "Q-NP" | Method == "Q-NP-REST"){
    
    newdata <- data.frame() 
    A1.vec <- rep(0,N2)
    A1.vec1 <- rep(1,N2)
    A1.vec2 <- rep(2,N2)
    newdata1=data.frame(X0,R1, g1.opt,A1=A1.vec)
    newdata2 <- data.frame(X0,R1, g1.opt,A1=A1.vec1)
    newdata3 <- data.frame(X0,R1, g1.opt,A1=A1.vec2)
    newdata <- rbind(newdata1,newdata2,newdata3)
    g1.tree2<-matrix(predict(QRF1,newdata=newdata), nrow=N2, ncol=3)
    max.Y1.pred <- E.R1.tree <- apply(g1.tree2, 1, max)
    g1.tree <- ifelse(g1.tree2[,1] > g1.tree2[,2] & g1.tree2[,1] > g1.tree2[,3], 0, 
                      ifelse(g1.tree2[,2] > g1.tree2[,1] & g1.tree2[,2] > g1.tree2[,3], 1,
                             ifelse(g1.tree2[,3] > g1.tree2[,1] & g1.tree2[,3] > g1.tree2[,2], 2,NA)))  
    
    select1[i]<-mean(g1.tree==g1.opt, na.rm=TRUE)
    R1.tree<-exp(1.5+0.3*X0[,id1]-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)
    EYs[i]<-mean(R1.tree, na.rm=TRUE)
    
    if(i%%50==0) print(i)
    
  } else(return("error-2method not specified"))
  
  
}
return(c(No.Stages = 1, No.RXperStage = 3, rho=rho, No.Var.H=No.Var.H, No.Var.Hsub=No.Var.Hsub, n.Sample = N, n.TestSet = N2, propensity.model = prop.model, true.DTR = "nontree", B = iter, Method = Method,
         MedianOpt = round(median(select1, na.rm=TRUE),3), IQROpt = round(IQR(select1, na.rm=TRUE),3), MedianEYs = round(median(EYs, na.rm=TRUE),3), IQREYs = round(IQR(EYs, na.rm=TRUE),3)))
} ## end function


