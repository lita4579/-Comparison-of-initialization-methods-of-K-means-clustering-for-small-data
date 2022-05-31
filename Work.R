#Bachelors thesis: Comparison of initialization methods of the K-Means clustering 
#algorithm for small data:


#Function to calculate the complexity
Complexity<-function(x,K){
  require(mclust)
  TrueCl<-rep(1:K,each=nrow(x)/K) 
  RealCenters<-matrix(0,K,ncol(x))
  for (i in 1:K){ 
    a<-((i-1)*(nrow(x)/K))+1
    b<-(i*(nrow(x)/K))
    RealCenters[i,]=apply(x[a:b,], 2, function(x){mean(x)})
  }
  RealCentersKMeans<-kmeans(x,RealCenters)
  adjRand<-adjustedRandIndex(TrueCl,RealCentersKMeans$cluster)
  return(adjRand)
}

#Initialization method based on the variable with highest variance
variancem<-function(x,K){
  Cl<-matrix(NA,K,ncol(x))
  varvec<-diag(var(x)) #Compute variance of each variable
  
  x<-x[order(x[,which.max(varvec)]),] #Sort X based on this variable
  
  for (j in 0:(K-1)){
    Cl[j+1,]=x[(nrow(x)*(2*j+1)/(2*K)),]
  }
  return(Cl)
}

#Initialization method based on highest distance to previously selected S_i centers
ppmethod<-function(x,K){
  Cl=matrix(NA,K,ncol(x))
  DistanceToCl<-matrix(0,K+2,nrow(x))
  Cl[1,]=x[round(runif(1,1,nrow(x))),]
  for (j in 2:K){
    DistanceToCl[j-1,]<-apply(x, 1, function(x){dist1(x,Cl[j-1,])})
    for (i in 1:nrow(x)){
      DistanceToCl[K+1,i]=min(DistanceToCl[1:(j-1),i])
    }
    DistanceToCl[K+2,]=DistanceToCl[K+1,]/sum(DistanceToCl[K+1,])
    Cl[j,]=x[match(sample(DistanceToCl[K+1,],1,prob = DistanceToCl[K+2,]),DistanceToCl[K+1,]),]
  }
  return(Cl)
}

#Euclidean squared distance
dist1<-function(x1,x2){
  d<-0
  for (i in 1:length(x1)){
    d=d+(x1[i]-x2[i])^2
  }
  return(d)
}

#Initialization method, where the first center is the mean, subsequent centers
Meanm<-function(x,K){
  Cl<-matrix(NA,K,ncol(x)) #Let # Clusters= # Variables
  mean=apply(x, 2, function(x){mean(x)})
  distancetomean<-apply(x, 1, function(x){dist1(x,mean)})
  x<-x[order(rank(distancetomean)),] #order
  for (i in 1:nrow(Cl)){
    Cl[i,]=x[1+(i-1)*nrow(x)/K,]
  }
  return(Cl)
}

#SSE given obs m is the j'th center.
SSE<-function(x,Cl){ #Given data, centers matrix, return sse
  SSE=0
  for (i in 1:nrow(x)){ 
    dobs2center<-rep(NA,nrow(Cl)-sum(is.na(Cl[,1]))) #Create vector of distance with # elements=#  centers
    for (j in 1:length(dobs2center)){
      dobs2center[j]=dist1(x[i,],Cl[j,]) #Distance of obs i to center j
    }
    SSE=SSE+min(dobs2center) 
  }
  return (SSE)
}

#Function that runs a Hierarchical clustering and cuts the tree at K clusters,
# then assigns the centers of those clusters as the centers 
HierMethod<-function(x,K){ 
  Cl<-matrix(NA,K,ncol(x))
  d<-hclust(dist(x)) #Input pairwise distance of observations, output hier.
  for (i in 1:K){ 
    Cl[i,]=colMeans(x[cutree(d,K) == i, , drop = FALSE])
  }
  return(Cl)
}

#Initialization method to assign each center randomly to an observation
RMethod<-function(x,K){
  Cl<-matrix(NA,K,ncol(x)) #Let # Clusters= # Variables
  for (i in 1:10000){ #Loop to make sure no duplicates from "runif()"
    randomvec<-round (runif(nrow(Cl),1,nrow(x)))
    if (any(duplicated(randomvec))==FALSE){ #If no duplicates, cut loop
      break
    }
  }
  for (i in 1:nrow(Cl)){ #Allocate each center as an observation
    Cl[i,]=x[randomvec[i],]
  }
  return(Cl)
}

#Function to plot the initial partition of centers
points<-function(N,V,c,K){
  if (V!=2){
    print("V has to be 2, as this is a 2d plot")
  }
  else{require(ggplot2)
    x<-si(N,V,o,K)
    m4<-as.data.frame(variancem(x,K))
    m1<-as.data.frame(RMethod(x,K))
    m5<-as.data.frame(ppmethod(x,K))
    m3<-as.data.frame(Meanm(x,K))
    m2<-as.data.frame(HierMethod(x,K))
    x<-as.data.frame(x)
    m<-rbind(m1,m2,m3,m4,m5)
    Method<-c(rep("Random",4),rep("Hierarchical",4),rep("Mean",4),rep("Variance",4),rep("PlusPlus",4))
    m<-cbind(m,Method)
    ggplot(x,aes(x=Var1,y=Var2))+
      geom_point()+geom_point(data = m, 
                              mapping = aes(x = V1, y = V2,color=Method), size=5)+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'))
  }
}

#Function to plot clusters in 2d
twodplot<-function(N,V,c,K){
  if (V!=2){
    print("V has to be 2, as this is a 2d plot")
  }
  else{require(ggplot2)
    x<-si(N,V,c,K)
    
    x<-as.data.frame(x)
    Cluster<-rep(1:K,each=N/K) 
    Cluster<-cut(Cluster,breaks=c(0.3,1.3,2.3,3.3,4.3,5.3,6.3,7.3,8.3),
                 labels=c("1","2","3","4","5","6","7","8"))
    x<-(cbind(x,Cluster))
    
    ggplot(x,aes(x=Var1,y=Var2,color=Cluster))+
      geom_point()+
      theme(legend.position = "none")
  }
}

#Function to plot clusters in 3d
threedplot<-function(N,V,c,K){
  if (V!=3){
    print("V has to be 3, as this is a 3d plot")
  }
  else{require(plotly)
    x<-si(N,V,c,K)
    x<-as.data.frame(x)
    Cluster<-rep(1:K,each=N/K) 
    Cluster<-cut(Cluster,breaks=c(0.3,1.3,2.3,3.3,4.3,5.3,6.3,7.3,8.3),
                 labels=c("1","2","3","4","5","6","7","8"))
    x<-(cbind(x,Cluster))
    plot_ly(x=x$Var1, y=x$Var2, z=x$Var3, type="scatter3d", mode="markers", showlegend = F,color=x$Cluster,size = 0.05)
  }
}

# Function to generate data, given parameters.
# Note: For specifying the number of variables, change the number of "a" in
# "expand grid" to the number of variables desired.
si<-function(N,V,c,K){
  if (K>2^V || N%%K!=0){ #If clusters>orthants or unbalanced cluster sizes
    print("Too many cl or N not divisible by K")
  }
  else{
    require(MASS)
    TrueCl<-rep(1:K,each=N/K) 
    x<-matrix(0,N,V) 
    a=c*c(1,-1)
    TrueCenters<-as.matrix(expand.grid(a,a,a)) #V=3
    TrueCenters<-TrueCenters[1:K,]
    for (i in 1:K){ 
      x<-x+mvrnorm(N,TrueCenters[i,],1*diag(V))*(TrueCl==i)
    }
    return(x)
  }
}

#specify parameters for sampling of data
#Sample m times, calculate the measures of performance for each IM m times,
# Output results: for 5 initialization methods and 4 criterion, 20 columns
an<-function(m,N,V,c,K){ #m=iter, N=x size, V=variables,o=overlap,K=clusters
  TrueCl<-rep(1:K,each=N/K) 
  MeanmMat<-matrix(NA,m,4)
  HierMethodMat<-matrix(NA,m,4)
  RMethodmat<-matrix(NA,m,4)
  ppMethodmat<-matrix(NA,m,4)
  variancemmat<-matrix(NA,m,4)
  require(mclust)
  set.seed(10)
  for (i in 1:m){ 
    x<-si(N,V,c,K)
    CentersMeanm<-Meanm(x,K)
    MeanmMat[i,1]=SSE(x,CentersMeanm) # SSE before
    kmeansmethod<-kmeans(x,CentersMeanm,iter.max=30)
    MeanmMat[i,2]=sum(kmeansmethod$withinss)
    MeanmMat[i,3]=adjustedRandIndex(TrueCl,kmeansmethod$cluster)
    MeanmMat[i,4]=kmeansmethod$iter
    CentersHierMethod<-HierMethod(x,K)
    HierMethodMat[i,1]=SSE(x,CentersHierMethod) # SSE before
    kmeansmethod<-kmeans(x,CentersHierMethod,iter.max=30)
    HierMethodMat[i,2]=sum(kmeansmethod$withinss)
    HierMethodMat[i,3]=adjustedRandIndex(TrueCl,kmeansmethod$cluster)
    HierMethodMat[i,4]=kmeansmethod$iter
    CentersRmethod<-RMethod(x,K)
    RMethodmat[i,1]=SSE(x,CentersRmethod) # SSE before
    kmeansmethod<-kmeans(x,CentersRmethod,iter.max=30)
    RMethodmat[i,2]=sum(kmeansmethod$withinss)
    RMethodmat[i,3]=adjustedRandIndex(TrueCl,kmeansmethod$cluster)
    RMethodmat[i,4]=kmeansmethod$iter
    Centersppmethod<-ppmethod(x,K)
    ppMethodmat[i,1]=SSE(x,Centersppmethod) # SSE before
    kmeansmethod<-kmeans(x,Centersppmethod,iter.max=30)
    ppMethodmat[i,2]=sum(kmeansmethod$withinss)
    ppMethodmat[i,3]=adjustedRandIndex(TrueCl,kmeansmethod$cluster)
    ppMethodmat[i,4]=kmeansmethod$iter
    CentersVariancem<-variancem(x,K)
    variancemmat[i,1]=SSE(x,CentersVariancem) # SSE before
    kmeansmethod<-kmeans(x,CentersVariancem,iter.max=30)
    variancemmat[i,2]=sum(kmeansmethod$withinss)
    variancemmat[i,3]=adjustedRandIndex(TrueCl,kmeansmethod$cluster)
    variancemmat[i,4]=kmeansmethod$iter
  }
  return(cbind(MeanmMat,HierMethodMat,RMethodmat,variancemmat,ppMethodmat))
}

# function to execute simulation, return mean if option=0, std.dev. if option=1 
results<-function(N,V,c,K,option){
  options(scipen=999) #Output decimal numbers
  m=1000
  if (option==0){
    r<-an(m,N,V,c,K)
    Mean<-cbind(colMeans(r[,1:4]),colMeans(r[,5:8]),colMeans(r[,9:12]),
                colMeans(r[,13:16]),colMeans(r[,17:20]))
    Mean<-as.data.frame(Mean)
    colnames(Mean)<-c("M","H","R","V","P")
    return(Mean)
  }
  if (option==1){
    r<-an(m,N,V,c,K)
    stddev<-cbind(sqrt(diag(var(r[,1:4]))),sqrt(diag(var(r[,5:8]))),sqrt(diag(var(r[,9:12]))),
                  sqrt(diag(var(r[,13:16]))),sqrt(diag(var((r[,17:20])))))
    stddev<-as.data.frame(stddev)
    colnames(stddev)<-c("MeanMat","HierMethod","R","variance","ppmethod")
    return(stddev)
  }
}
