numbparents<-rpois(1,lambdaParent*dim)

xxParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
yyParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)

numbdaughter<-rpois(numbparents,(lambdaDaughter))
sumdaughter<-sum(numbdaughter)

phi<-2*pi*runif(sumdaughter)
rho<-radiusCluster*sqrt(runif(sumdaughter))

xx0=rho*cos(phi)
yy0=rho*sin(phi)

xx<-rep(xxParent,numbdaughter)
yy<-rep(yyParent,numbdaughter)

xx<-xx+xx0

yy<-yy+yy0
cds<-data.frame(xx,yy)
is_outlier<-function(x){
  x > dim| x < 0
}
cds<-cds[!(is_outlier(cds$xx)|is_outlier(cds$yy)),]
sampleselect<-sample(1:nrow(cds),hosts,replace=F)
cds<-cds%>%slice(sampleselect)

randfunction<-function(x){
  x<-runif(length(x),0,dim)
}
randselect<-sample(1:nrow(cds),floor(hosts*randmod),replace=F)
cds[randselect,]<-apply(cds[randselect,],1,randfunction)
landscape<-ppp(x=cds$xx,y=cds$yy,window=owin(xrange=c(0,dim),yrange=c(0,dim)))
ggplot(data.frame(landscape))+geom_point(aes(x=x,y=y))+coord_equal()+theme_minimal()
