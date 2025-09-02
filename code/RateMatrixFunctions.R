### Ancillary functions ####


### make start values function to get regime specific starting values for rate matrix ###

make.start.values.regimes <- function(phy, d, regimes) {
  
  reg <- sort(unique(regimes))
  r.list<-list()
  for(i in 1:length(reg)) {
    phy.tmp <- drop.tip(phy, names(regimes)[which(!regimes %in% reg[i])])
    d.tmp <- treedata(phy.tmp, d)$data
    r.list[[i]]<- mvBM(phy.tmp, d.tmp, model="BM1", method="pic")$sigma
    #r.list[[i]]<-ratematrix(phy.tmp, d.tmp)
  }
  names(r.list) <- reg
  
  list(colMeans(d), lapply(r.list, cov2cor), lapply(lapply(r.list, diag), sqrt))
  
}


### rate matrix acceptance rates specific to individual regimes - check for good mixing ###

acceptance.rates <- function(handle, burn = 0) {
  yy <- read.table(paste(handle$outname, handle$ID, "log",sep="."), header=T, sep=";")
  if(!burn ==0) {
    gens <- nrow(yy)
    yy <- yy[-(1:(gens*burn)),]
  }
  n.regimes <- unique(yy$matrix.corr)[-1]
  for(i in 1:length(n.regimes)) {
  cat(paste("regime " , n.regimes[i], " correlation acc. = ", round(length(which(yy$accepted==1 & yy$matrix.corr==n.regimes[i])) / length(which(yy$matrix.corr==n.regimes[i])), 2),"\n", sep=""))
  cat(paste("sd acc." , n.regimes[i], " = ", round(length(which(yy$accepted==1 & yy$matrix.sd==n.regimes[i])) / length(which(yy$matrix.sd==n.regimes[i])),2),"\n", sep=""))
  }
  cat(paste("root acc. = ", round(length(which(yy$accepted==1 & yy$root==1)) / length(which(yy$root==1)),2),"\n", sep=""))
}


make.2d.mat <- function(x, index) {
  
  tmp <- matrix(data=NA, nrow=2, ncol=2)
  diag(tmp) <- diag(x)[index]
  tmp[upper.tri(tmp)] <- tmp[lower.tri(tmp)] <- x[index[1], index[2]]
  tmp
}

### plot rate matrix posteriors in 2D ###

plot.2d.ellipse <-function(tree, data, posterior, traits, level=0.5, xlab="trait1", ylab ="trait1", main =NULL, plot.mle= F, cex.lab=1, cex.axis=1, xlim = NULL, ylim=NULL) {
  
  xx <- lapply(posterior, make.2d.mat, index= traits)
  
  ellipse.dat <-lapply(xx, ellipse, level= level)
  
  plot(range(unlist(ellipse.dat)), range(unlist(ellipse.dat)), 
       type="n", xlab=xlab, ylab=ylab, main=main, cex.lab=cex.lab, cex.axis=cex.axis, xlim=xlim, ylim=ylim, xaxt='n', ann=FALSE, yaxt='n')
  abline(h=0, v=0, col="gray")
  x<-sample(seq(1, length(ellipse.dat)), size = 50)
  for(i in 1:length(x)){
    lines(ellipse.dat[[x[i]]], col="lightblue", lwd=0.5)
  }
  if(plot.mle==TRUE) {
    mle <- mvMORPH::mvBM(tree, data, model = "BM1", method = "pic" )
    lines(ellipse(make.2d.mat(mle$sigma, index= traits), level=level), col="black")
  }
}




### approximate the integral in Ovaskainen's D by repeated sampling from f(x) ###

Ovaskainens_D <- function(m1, m2, n=1000) {
  xi <- mvtnorm::rmvnorm(n, rep(0, dim(m1)[1]), m1)
  fx <- mvtnorm::dmvnorm(xi, rep(0, dim(m1)[1]), m1)
  gx <- mvtnorm::dmvnorm(xi, rep(0, dim(m1)[1]), m2)
  un.transformed.dist <- 1-(2*mean(gx/(fx+gx)))
  if(un.transformed.dist<0) un.transformed.dist<-0 ## catch instances where distance is negative (?)
  sqrt(un.transformed.dist)
  
}

## compute the distribution of Ovaskainen's D using two posteiror distributions of rate matrices ##

Ovaskainen_distance_test <- function(distribution1, distribution2, n.reps.integral=1000, n.comps=1000) {
  
 joint.post.distances <-  mapply(function(X,Y) {
   Ovaskainens_D(X, Y, n=n.reps.integral)}, X=distribution1, Y=distribution2)
  
  
  dist.2.matrices <- numeric()
  
  upper <- min(c(length(distribution1), length(distribution2)))
  for(i in 1:n.comps) {
    ab <- sample(x=1:upper,size = 2,replace = F)
    dist.2.matrices[i] <- (Ovaskainens_D(distribution1[[ab[1]]], distribution1[[ab[2]]], n=n.reps.integral) + 
                             Ovaskainens_D(distribution2[[ab[1]]], distribution2[[ab[2]]], n=n.reps.integral)) -
      (joint.post.distances[ab[1]] + joint.post.distances[ab[2]]) 
    
  }
  return(list(distances = joint.post.distances, 
              psi.AB = dist.2.matrices, 
              P.same = round((length(which(dist.2.matrices>=0)) / length(dist.2.matrices)),2)
         ))
}

## utility for computing Trait Change Index ##

get.var.expl <- function(matrix, eigenvectors) {
  apply(eigenvectors, 2, function(x) t(x) %*% matrix %*% x)
}
#######

### Robinson and Beckerman's trait change index, examining the difference in the proportion of variance explained ###
### in matrix 2 by the eigenvectors of matrix 1 and the same vectors in matrix 1. If the 95% credible interval    ###
### includes 0, then the hypothesis that the same combinations of traits explain variance cannot be rejected      ###

TCI <- function(matrix1, matrix2, n.vectors) {
  tci <- numeric()
  z<-min(c(length(matrix1), length(matrix2)))
  for(i in 1:z) {
    e1 <- eigen(matrix1[[i]])
     e1.explained <- 100 * (sum(e1$values[1:n.vectors]) / sum(e1$values))
    m2.in.e1 <- get.var.expl(matrix = matrix2[[i]], eigenvectors = e1$vectors)
    
    e2.explained<-100 * (sum(m2.in.e1[1:n.vectors]) /sum(m2.in.e1))
    tci[i] <- round(e1.explained,3) - round(e2.explained,3)
  }
  tci
}


####################################################################
## function to perform Kirkpatrick's (2009) test of the number of ##
## effective dimensions and Robinson and Beckerman's 2013 test of ##  
## the number of significant eigenvectors, both using a Bayesian  ##
## posterior distribution of rate matrices                        ##
####################################################################
n.dims <- function(posterior) {

  n.effective.dims <- function(x) {
    lambda <- eigen(x)$value
    sum(lambda) / lambda[1]
  }
  
  kirkpatric.dim <- unlist(lapply(posterior, n.effective.dims))
  list(av.dim = mean(kirkpatric.dim), quantiles = quantile(kirkpatric.dim, c(0.025, 0.975)))  

}

### Robinson and beckerman version 
  
significant.eigenvectors <- function(posterior, nreps=1000) {
    
    no.traits <- dim(posterior[[1]])[1]
    Peig<-matrix(0,nreps,(no.traits-1))
    
    for(i in 1:nreps) {
      xx <- sample(seq(1, length(posterior),2))
      d1 <- posterior[[xx[1]]]
      d2 <- posterior[[xx[2]]]
      
      env1_1<-eigen(d1)
      env1_2<-eigen(d2)
      
      # No Sig Eigenvectors for each matrix
      
      total1<-sum(env1_1$values)
      total2<-sum(env1_2$values)
      
      for(j in 1:(no.traits-1)){
        idx<-1:j
        tmp1<-env1_1$values[idx]
        tmp2<-env1_2$values[idx]
        Peig[i,][j]<-(((sum(tmp1)-sum(tmp2)) + (total1-total2)) - ((sum(tmp1)-total1) + (sum(tmp2) - total2)))
      }
    }
    

    pp <- apply(Peig, 2, function(x) 1-length(which(x<0))/ length(x))
    return(list(test= Peig, posterior.probability = pp))
  }
  
  
  hdr.extract <- function(x) {
    hpd <- hdr(x)
    c(hpd$mode, hpd$hdr["95%", ])
  }
  
#### 
 
  evolvability <- function(x) {
    sqrt(eigen(x)$values[1])
  }
  
  eccentricity <- function(x) {
    lambda <- eigen(x)$values
    lambda[1] / lambda[2]
  }
  

  rSDE <- function(x) {
    p <- nrow(x)
    xx <- eigen(x)
    sqrt((p*sum((xx$values - mean(xx$values))^2)) / 
      ((p-1)*((sum(xx$values) * sum(xx$values)))))
  }

  add.zeros.cat <- function(x) {
    x <- cbind(rep(0,6), rbind(rep(0, 5), x))
    x <- rbind(x[1:2, ], rep(0,6), x[3:6, ])
    x<- cbind(x[,1:2], rep(0,7), x[,3:6])
    x
  }
  
  
rate.of.dispersion <- function(x)  {
  sum(diag(x))
} 



angle.bw.principal.vectors <- function(matrix1, matrix2, nreps.sig) {

  xx <- lapply(matrix1, eigen)
  xx <- lapply(xx, function(x) abs(x$vectors[,1]))
  
  yy <- lapply(matrix2, eigen)
  yy <- lapply(yy, function(x) abs(x$vectors[,1]))
  
  dot.prod <- numeric()
  for(i in 1:length(xx)) {
    dot.prod[i] <- sum(xx[[i]]*yy[[i]])
  }
  angles <- acos(dot.prod) * (180/pi)
  
  sig.test <- numeric()
  upper <- min(c(length(matrix1), length(matrix2)))
  for(i in 1:(nreps.sig)) {
    ab <- sample(x=1:upper,size = 2,replace = F)
    xx1 <- abs(eigen(matrix1[[ab[1]]])$vectors[,1])
    xx2 <- abs(eigen(matrix1[[ab[2]]])$vectors[,1])
    yy1 <- abs(eigen(matrix2[[ab[1]]])$vectors[,1])
    yy2 <- abs(eigen(matrix2[[ab[2]]])$vectors[,1])
    
    find.angle <- function(x,y)     acos(sum(x*y)) * (180/pi)
    
    sig.test[i] <- ((find.angle(xx1, xx2) + find.angle(yy1, yy2)) - (angles[ab[1]] + angles[ab[2]]))
    
  }
    
  return(list(angles=angles, P.same = (length(which(sig.test>=0)) / length(sig.test)) , sig.test=sig.test))
  
}

### perform pairwise comparison of rate matrix posteriors  ###

pairwise.compare.matrix <- function(X) {
  if(is.list(X) == FALSE) stop("the posterior sample needs to be in the form of a list")

  comp <- combn(seq(1, length(X)), m=2)
  
  for(i in 1: ncol(comp)) {
    xx<- table(apply(cbind(X[[comp[1,i]]], X[[comp[2,i]]]), 1, which.max))
    if(length(xx)==1) {
      yy <- setNames(c(0,0), c(1,2))
      yy[c(1,2) %in% names(xx) ] <- xx
      xx<-yy
    }
    print((paste("posterior probability", names(X)[comp[which.max(xx),i]], "is greater than", 
          names(X)[comp[which.min(xx),i]], "=", round(xx[which.max(xx)]/sum(xx),2))))
   }
    
}


pairwise.TCI <- function(X) {
  if(is.list(X) == FALSE) stop("the posterior sample needs to be in the form of a list")
  
  comp <- combn(seq(1, length(X)), m=2)
  res <- list()
  for(i in 1: (2*ncol(comp))) {
    if(i <= ncol(comp)) {
    res[[i]] <- TCI(X[[comp[1,i]]],X[[comp[2, i]]] , n.vectors = 1)
    } else{
      res[[i]] <- TCI(X[[comp[2,(i-ncol(comp))]]],X[[comp[1, (i-ncol(comp))]]] , n.vectors = 1)
    }
  }
  names(res) <- c(apply(comp, 2, function(x) paste("TCI", names(X)[x[2]], "by",names(X)[x[1]], sep=" ")), apply(comp, 2, function(x) paste("TCI", x[1], "by" ,x[2], sep=" ")))
  res
}

pairwise.angles <- function(X) {
  if(is.list(X) == FALSE) stop("the posterior sample needs to be in the form of a list")
  
  comp <- combn(seq(1, length(X)), m=2)
  res <- list()
  for(i in 1: ncol(comp)) {
      tmp <- angle.bw.principal.vectors(X[[comp[1,i]]], X[[comp[2,i]]],nreps.sig = 1000)
      print(paste("probability that angles bw", names(X)[comp[1,i]], "and", 
                   names(X)[comp[2,i]], "are same =", tmp$P.same))
      res[[i]] <- tmp$angles
  }   
  names(res) <- paste("angle_between_", names(X)[comp[1,]], "_and_", names(X)[comp[2,]], sep="")
  res
    }


pairwise.Ovaskainen <- function(X) {
  comp <- combn(seq(1, length(X)), m=2)
  res <- list()
  for(i in 1: ncol(comp)) {
  ovd <- Ovaskainen_distance_test(X[[comp[1,i]]], X[[comp[2,i]]],n.comps = 1000)
  print(paste("probability that matrices", names(X)[comp[1,i]], "and", 
              names(X)[comp[2,i]], "differ =", ovd$P.same))
  res[[i]] <- ovd$psi.AB
  }
  names(res) <- paste("distance_between_", names(X)[comp[1,]], "_and_", names(X)[comp[2,]], sep="")
  res
}


