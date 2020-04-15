library(KernSmooth)
library(akima)
library(glmnet)

denPre2D <- function(edge,classLabel,DataFit,newData,method=c("integers","bandwidth")[2]){
  
  Index0 <- which(classLabel==0)
  Index1 <- which(classLabel==1)
  x <- DataFit[c(Index0,Index1),edge[1]]
  y <- DataFit[c(Index0,Index1),edge[2]]
  x0 <- newData[,edge[1]]
  y0 <- newData[,edge[2]]
  N <- max(c(x,y))
  coI <- 1:length(Index0)
  caI <- length(Index0)+(1:length(Index1))
  
  caWidth <- c(bw.nrd0(x[caI]),bw.nrd0(y[caI]))
  coWidth <- c(bw.nrd0(x[coI]),bw.nrd0(y[coI]))
  gridSize <- switch(method,
                     integers  = c(N, N),
                     bandwidth = ceiling(N / c(min(caWidth[1], coWidth[1]), 
                                               min(caWidth[2], coWidth[2]))))
  gridSize <- pmax(gridSize,10) # make sure there are at least 100 points in total

  caSmooth <- bkde2D(x=cbind(x[caI], y[caI]), bandwidth=caWidth,  gridsize=gridSize)
  caP <- caSmooth$fhat
  coSmooth <- bkde2D(x=cbind(x[coI], y[coI]), bandwidth=coWidth, gridsize=gridSize)
  coP <- coSmooth$fhat

  # make sure there are no zeros in the smooth function (since we will take a log of that)
  # caP[caP==0] <- min(caP[caP>0])/100
  caP <- pmax(caP, 1e-10)
  # coP[coP==0] <- min(coP[coP>0])/100
  coP <- pmax(coP, 1e-10)

  cafit <- bicubic(x=caSmooth$x1, y=caSmooth$x2, z=caP, x0=x0,y0=y0)
  cofit <- bicubic(x=coSmooth$x1, y=coSmooth$x2, z=coP, x0=x0,y0=y0)
  cafit$z <- pmax(cafit$z, 1e-10)
  cofit$z <- pmax(cofit$z, 1e-10)
  denPre <- log(cafit$z/cofit$z)
  denPre
}



jdinac <- function(EDGE,classLabel,DataFit,DataPre,zFit=NULL,zPre=NULL,nsplit=10,nfolds=5){
  
  if(missing(DataPre)) {
    DataPre <- DataFit
    zPre <- zFit
  }  
  if(!is.null(zFit) & !is.null(zPre)){
    Pre.Z <- jdinac.z(EDGE=EDGE,classLabel=classLabel,DataFit=DataFit,DataPre=DataPre,
                        zFit=zFit,zPre=zPre,nsplit=nsplit,nfolds=nfolds)
    yPre <- Pre.Z$yPre
    Eset <- Pre.Z$Eset
  } else if(nsplit<=0){       
    stop("nsplit must be positive integer")    
  } else {
    preY <- NULL
    vset <- NULL
    for(i in 1:nsplit){
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      splitid <- c(sample(1:size0,sn0),sample((1:size1)+size0,sn1))
    
      cLabel <- classLabel[splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[splitid,],
                    newData=rbind(DataFit[-splitid,],DataPre))
      y <- classLabel[-splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      vset <- c(vset,coefs)
    
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
      y <- classLabel[splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      vset <- c(vset,coefs)     
    } 
    yPre <- rowMeans(preY) 
    numb <- table(vset)
    Vid <- as.numeric(rownames(numb))  
    Eset <- cbind(EDGE[Vid,],numb)
    Eset <- Eset[order(Eset[,3],decreasing=T),]
    colnames(Eset) <- c("row","col","numb")
  } 
  list(yPre=yPre,Eset=Eset) 
}



jdinac.z <- function(EDGE,classLabel,DataFit,DataPre,zFit,zPre,nsplit=10,nfolds=5){
  
  if(class(zFit)!="matrix" | class(zPre)!="matrix")
    stop("zFit and zPre must be 'matrix'")
  if(nsplit<=0){       
    stop("nsplit must be positive integer")  
  } else {
    preY <- NULL
    vset <- NULL
    for(i in 1:nsplit){
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      splitid <- c(sample(1:size0,sn0),sample((1:size1)+size0,sn1))
    
      cLabel <- classLabel[splitid]
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[splitid,],
                    newData=rbind(DataFit[-splitid,],DataPre))
      y <- classLabel[-splitid]
      zfit <- zFit[-splitid, ,drop=F]
      cv.fit <- cv.glmnet(x=cbind(zfit,denX[1:length(y),]), y=y, family = "binomial", nfolds=nfolds) 

      yp <- predict(cv.fit,newx=cbind(zPre,denX[-(1:length(y)),]), s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-(1:(1+ncol(zFit)))] !=0)
      vset <- c(vset,coefs)
    
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
      y <- classLabel[splitid]
      zfit <- zFit[splitid, ,drop=F]
      cv.fit <- cv.glmnet(x=cbind(zfit,denX[1:length(y),]), y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=cbind(zPre,denX[-(1:length(y)),]), s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-(1:(1+ncol(zFit)))] !=0)
      vset <- c(vset,coefs)     
    } 
    yPre <- rowMeans(preY) 
    numb <- table(vset)
    Vid <- as.numeric(rownames(numb))  
    Eset <- cbind(EDGE[Vid,],numb)
    Eset <- Eset[order(Eset[,3],decreasing=T),]
    colnames(Eset) <- c("row","col","numb")
  } 
  list(yPre=yPre,Eset=Eset) 
}

