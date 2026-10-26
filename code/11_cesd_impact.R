source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages({library(mice);library(nnet);library(dplyr);library(haven)})
m <- readRDS("out/jamda_full.rds")
# 正确计分：dc013、dc016 反向（5 - x）
CESD <- paste0("dc01", 1:8)
h <- read_dta("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", col_select = all_of(c("ID", CESD)))
Cm <- sapply(h[CESD], function(s) as.numeric(zap_labels(s)))
Cm[, c("dc013","dc016")] <- 5 - Cm[, c("dc013","dc016")]
cs <- data.frame(ID = as.character(h$ID), cesd_fix = rowSums(Cm, na.rm=TRUE))
cs$cesd_fix[rowSums(!is.na(Cm)) < 6] <- NA
cs <- cs %>% distinct(ID, .keep_all=TRUE)
m <- m %>% left_join(cs, by="ID")
cat("cesd 原:", sprintf("%.2f (%.2f)", mean(m$cesd,na.rm=TRUE), sd(m$cesd,na.rm=TRUE)),
    " 修正后:", sprintf("%.2f (%.2f)", mean(m$cesd_fix,na.rm=TRUE), sd(m$cesd_fix,na.rm=TRUE)),
    " 相关 r=", sprintf("%.3f", cor(m$cesd, m$cesd_fix, use="complete.obs")), "\n")
saveRDS(m, "out/jamda_full_fix.rds")

COV2 <- c("age10","male","bmi","nchronic","cesd_fix","smoke","edu1","edu2","edu3","edu4")
b0 <- get_base(m, "grp19")
d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13>0) %>%
      mutate(edu_f=factor(edu), wsc=w13/mean(w13,na.rm=TRUE))
imp <- mice(d0[, c("age10","male","bmi","nchronic","cesd_fix","smoke","edu_f","grp","y")],
            m=20, seed=20261023, printFlag=FALSE)
mn_vcov <- function(beta, X, y, w, cluster){
  K<-2; p<-ncol(X); N<-nrow(X); eta<-X%*%t(beta); P<-exp(cbind(0,eta)); P<-P/rowSums(P)
  Y<-sapply(1:2, function(k) as.numeric(y==k)); A<-matrix(0,K*p,K*p)
  Wm<-list(P[,2]*(1-P[,2]),P[,3]*(1-P[,3]),-P[,2]*P[,3])
  bl<-matrix(c(1,3,3,2),2,2,byrow=TRUE)
  for(k in 1:K) for(l in 1:K){ ww<-w*Wm[[bl[k,l]]]; A[((k-1)*p+1):(k*p),((l-1)*p+1):(l*p)]<-crossprod(X,X*ww) }
  S<-matrix(0,N,K*p); for(k in 1:K) S[,((k-1)*p+1):(k*p)]<-(w*(Y[,k]-P[,k+1]))*X
  Sc<-rowsum(S,group=cluster,reorder=FALSE); G<-nrow(Sc)
  list(V=(G/(G-1))*((N-1)/(N-K*p))*solve(A)%*%crossprod(Sc)%*%solve(A))
}
bs<-list(); vs<-list(); p<-14
for(i in 1:20){
  di<-complete(imp,i); di$y<-factor(di$y,levels=c("alive_free","disabled","dead"))
  di$grp<-relevel(factor(di$grp,levels=GRP),ref="never"); di$edu<-as.numeric(as.character(di$edu_f))
  for(k in 1:4) di[[paste0("edu",k)]]<-as.numeric(di$edu==k)
  di$communityID<-d0$communityID; di$wsc<-d0$wsc
  X<-cbind(1,as.matrix(di[,COV2]),as.numeric(di$grp=="remitted"),as.numeric(di$grp=="incident"),as.numeric(di$grp=="persistent"))
  colnames(X)<-c("const",COV2,paste0("grp",c("remitted","incident","persistent")))
  f<-as.formula(paste("y ~",paste(c(COV2,"grp"),collapse="+")))
  fit<-multinom(f,data=di,weights=di$wsc,trace=FALSE,maxit=500)
  bb<-matrix(coef(fit),2,byrow=FALSE); colnames(bb)<-colnames(X)
  vc<-mn_vcov(bb,X,as.integer(di$y)-1,di$wsc,di$communityID)
  bs[[i]]<-as.numeric(t(bb)); vs[[i]]<-as.numeric(diag(vc$V))
}
b<-Reduce(`+`,bs)/20; W<-Reduce(`+`,vs)/20
Bv<-Reduce(`+`,lapply(bs,function(z)(z-b)^2))/19; T<-W+(1+1/20)*Bv
cat("\n=== 修正 CES-D 后的完全校正 RRR（加权 MI，2019 切点）===\n")
for(k in 1:2) for(j in 12:14){
  est<-b[(k-1)*p+j]; se<-sqrt(T[(k-1)*p+j])
  cat(sprintf("%-9s %-11s %.2f (%.2f-%.2f)\n", c("disabled","dead")[k], c("remitted","incident","persistent")[j-11],
      exp(est), exp(est-1.96*se), exp(est+1.96*se)))
}
