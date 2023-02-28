fitcurve<-function(h, d, model, start){
  if ( any( is.na(h) ) )  stop("NAs values are not allowed for height data.")
  if ( any( is.na(d) ) )  stop("NAs values are not allowed for dbh data.")
  if(model != "weibull" & model != "probola" & model != "chapman-richards" & model != "logistic" & model != "prodan" & 
	model != "gompertz" & model != "korf" & model != "sibbesen" & model != "katkowsky" & model != "hossfeldiv" )
	stop("model's name is not implemented or misspelled. Please check the manual for guidelines.")
  if(length(start) != 3) stop("length of vector of starting values must be three.")
  x <- seq( min(d), max(d), 0.01 )
  a <- min(d) + 11
  b <- max(h) - 1
  k <- 1.25
  if(model=="weibull"){
    relation<-as.formula(h~1.3+beta1*(1-exp(-beta2*d^beta3)))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1*(1-exp(-b2*x^b3))
    }
    g<-expression(H==paste(1.3+beta[1](1-e^{-beta[2]*D^beta[3]})))
  }
  if(model=="probola"){
    relation<-as.formula(h~1.3+beta1*(1-exp(-beta2*d))^beta3)
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1*(1-exp(-b2*x))^b3
    }
    g<-expression(H==paste(1.3+beta[1](1-e^{beta[2]*D})^beta[3]))
  }
  if(model=="chapman-richards") {
    relation<-as.formula(h~1.3+beta1+beta2/(d+beta3))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1+b2/(x+b3)
    }
    g<-expression(H==1.3+beta[1]+paste(frac(beta[2],D+beta[3])))
  }
  if (model=="logistic"){
    relation<-as.formula(h~1.3+beta1/(1+beta2*exp(-beta3*d)))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1/(1+b2*exp(-b3*x))
    }
    g<-expression(H==paste(1.3+frac(beta[1],1+beta[2]*e^{-beta[3]*D})))
  }
  if (model=="prodan"){
    relation<-as.formula(h~1.3+d^2/(beta1*d^2+beta2*d+beta3))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+x^2/(b1*x^2+b2*x+b3)
    }
    g<-expression(H==paste(1.3+frac(D^2,beta[1]*D^2+beta[2]*D+beta[3])))
  }
  if (model=="gompertz"){
    relation<-as.formula(h~1.3+beta1*exp(-beta2*exp(-beta3*d)))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1*exp(-b2*exp(-b3*x))
    }
    g<-expression(H==paste(1.3+beta[1]*e^{-beta[2]*e^{-beta[3]*D}}))
  }
  if (model=="korf"){
    relation<-as.formula(h~1.3+beta1*exp(-beta2/d^beta3))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1*exp(-b2/x^b3)
    }
    g<-expression(H==paste(1.3+beta[1]*e^{-beta[2]*D^{-beta[3]}}))
  }
  if (model=="sibbesen"){
    relation<-as.formula(h~1.3+beta1*d^(beta2/d^beta3))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1*x^(b2/x^b3)
    }
    g<-expression(H==paste(1.3+beta[1]*D^{beta[2]*D^{-beta[3]}}))
  }
  if (model=="katkowsky"){
    relation<-as.formula(h~1.3+beta1*e^(-beta2/(d+beta3)))
    f<-function(x,par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1*exp(-b2/(x+b3))
    }
    g<-expression(H==paste(1.3+beta[1]*e^{-frac(beta[2],D+beta[3])}))
  }
  if (model=="hossfeldiv"){
    relation<-as.formula(h~1.3+beta1/(1+1/(beta2*d^beta3)))
    f<-function(x, par){
      b1<-par[1]
      b2<-par[2]
      b3<-par[3]
      1.3+b1/(1+1/(b2*x^b3))
    }
    g<-expression(H==paste(1.3+frac(beta[1],1+frac(1,beta[2]*D^beta[3]))))
  }

it <- 1
 i <- 0
while(it <= 1){ 
random1 <- runif(1, min(0, start[1] - 5), start[1] + 5)
random2 <- runif(1, min(0, start[2] - 5), start[2] + 5)
random3 <- runif(1, min(0, start[3] - 5), start[3] + 5)
out1 <- tryCatch( summary( nls( relation, start = list( beta1 = random1, beta2 = random2, beta3 = random3 ) ) ), 
	error=function(e)( "fail" )  )
if( out1[1] == "fail" ){
	it <- 1
	}else{
	it <- 2
	}
	i <- i + 1
}
  out2 <- plot(d, h, main = "Height Vs. Diameter", xlab = "Diameter (cm)", ylab = "Height (m)", cex = k, cex.lab = k, cex.axis = k, col = 'black', lwd = k )
  lines(x, f(x, out1$parameters[, 1]), col = 'blue', cex = 0.5)
  text(a,b,g)
  return( list( estimate = out1$parameters, residuals = out1[[2]], 
"var-cov" = out1[[5]], "residual Std. Error" = out1[[3]], iteration = i) )
  out2
}
fitmixturegrouped<-function(family,r,f,K,initial=FALSE,starts){
  m<-length(r)
  N<-1000
  cri<-0.005
  h<-2
  eps<-1
  R<-r[m]-r[1]
  prob<-E<-anderson<-cramer<-F.hat<-pgks<-Fplus.hat<-Fminus.hat<-rep(0,m-1)
  pg<-rep(0,m)
  n<-sum(f)
  edf<-cumsum(f)/n
  mid<-(r[-m]+r[-1])/2
  xmid<-rep(mid,f)
  omega<-alpha<-beta<-lambda<-rep(0,K)
  alpha.matrix<-beta.matrix<-lambda.matrix<-omega.matrix<-matrix(1,nrow=N,ncol=K)
  clust<-kmeans(xmid,K)
  if(family=="log-normal"){
    F<-function(x,a,b,omega){
      d<-matrix(0,nrow=length(x)-1,ncol=length(a))
      for(k in 1:length(a))d[,k]<-diff(plnorm(x,meanlog=a[k],sdlog=b[k]))
      return(omega%*%t(d))
    }
    Fi<-function(x,a,b) diff(plnorm(x,meanlog=a,sdlog=b))
    int.beta<-Vectorize(function(w,a,b)integrate(function(x)(log(x)-a)^2*dlnorm(x,meanlog=a,sdlog=b),lower=0,upper=w)$value,"w")
    if (initial==FALSE){
      for(k in 1:K){
        omega.matrix[1,k]<-sum(clust$cluster==k)/n
        group.x<-xmid[clust$cluster==k]
        alpha.matrix[1,k]<-log(median(group.x))
        beta.matrix[1,k]<-sqrt(2*abs(log(mean(group.x)/median(group.x))))
      }
    }
    if(initial==TRUE){
      omega.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
    }
    while (eps>cri && h<=N){
      for (k in 1:K){
        ss<-Fi(r,alpha.matrix[h-1,k],beta.matrix[h-1,k])/F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,])
        omega.matrix[h,k]<-1/n*omega.matrix[h-1,k]*sum(f*ss,na.rm=TRUE)
        rr<-(log(r)-alpha.matrix[h-1,k])/beta.matrix[h-1,k]
        num.alpha<-(-beta.matrix[h-1,k]/sqrt(2*pi)*diff(exp(-rr^2/2))+alpha.matrix[h-1,k]*diff(pnorm(rr)))/
          F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,])
        alpha.matrix[h,k]<-omega.matrix[h-1,k]*sum(f*num.alpha,na.rm=TRUE)/(n*omega.matrix[h,k])
        num.beta<-diff(int.beta(r,alpha.matrix[h-1,k],beta.matrix[h-1,k]))/F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,])
        beta.matrix[h,k]<-sqrt(omega.matrix[h-1,k]*sum(f*num.beta,na.rm=TRUE)/(n*omega.matrix[h,k]))
      }
      eps<-sum(c(abs(omega.matrix[h-1,]-omega.matrix[h,])),c(abs(alpha.matrix[h-1,]-alpha.matrix[h,])),c(abs(beta.matrix[h-1,]-beta.matrix[h,])))
      h<-h+1
    }
    omega<-round(omega.matrix[h-1,],digits=7)
    omega[K]<-1-(sum(omega)-omega[K])
    alpha<-alpha.matrix[h-1,]
    beta<-beta.matrix[h-1,]
  }
  if(family=="weibull"){
    F<-function(x,a,b,omega){
      d<-matrix(0,nrow=length(x)-1,ncol=length(a))
      for(k in 1:length(a)){
        d[,k]<-diff(pweibull(x,shape=a[k],scale=b[k]))
      }
      return(omega%*%t(d))
    }
    Fi<-function(x,a,b) diff(pweibull(x,shape=a,scale=b))
    if (initial==FALSE){
      for(k in 1:K){
        omega.matrix[1,k]<-sum(clust$cluster==k)/n
        group.x<-xmid[clust$cluster==k]
        h.grouped<-function(u){(gamma(1+2/u)-(gamma(1+1/u))^2)/(gamma(1+1/u))^2-var(group.x)/(mean(group.x))^2}
        alpha.matrix[1,k]<-uniroot(h.grouped,lower=0.02,upper=500000)$root
        beta.matrix[1,k]<-mean(group.x)/gamma(1+1/alpha.matrix[1,k])
      }
    }
    if(initial==TRUE){
      omega.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
    }
    while (eps>cri && h<=N){
      for (k in 1:K){
        ss<-omega.matrix[h-1,k]*Fi(r,alpha.matrix[h-1,k],beta.matrix[h-1,k])/F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,])
        omega.matrix[h,k]<-1/n*sum(f*ss,na.rm=TRUE)
        num<-sum(f*Fi(r,alpha.matrix[h-1,k],beta.matrix[h-1,k]))
        rr<-(r/beta.matrix[h-1,k])^alpha.matrix[h-1,k]
        denom<--omega.matrix[h-1,k]*sum(f*diff((1+rr*log(rr))*exp(-rr))/F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,]))
        alpha.matrix[h,k]<-n*omega.matrix[h,k]*alpha.matrix[h-1,k]/denom
        rr<-(r/beta.matrix[h-1,k])^alpha.matrix[h,k]
        numb<-beta.matrix[h-1,k]^alpha.matrix[h,k]*omega.matrix[h-1,k]*sum(f*diff(-exp(-rr))/
                                                                             F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,]))
        beta.matrix[h,k]<-(numb/(n*omega.matrix[h,k]))^(1/alpha.matrix[h,k])
      }
      eps<-sum(c(abs(omega.matrix[h-1,]-omega.matrix[h,])),c(abs(alpha.matrix[h-1,]-alpha.matrix[h,])),c(abs(beta.matrix[h-1,]-beta.matrix[h,])))
      h<-h+1
    }
    omega<-round(omega.matrix[h-1,],digits=7)
    omega[K]<-1-(sum(omega)-omega[K])
    alpha<-alpha.matrix[h-1,]
    beta<-beta.matrix[h-1,]
  }
  if(family=="gamma"){
    intalpha<- function(x, a, b) {
      y <- rep(0, length(f))
      g1 <- function(z) {
        integrand<-log(z/b)*dgamma(z,shape=a,scale=b)
        return(ifelse(z==0 | z==Inf,0,integrand))
      }
      for (i in 2:(length(x))) {
        y[i - 1] <-suppressWarnings(integrate(g1, lower=x[i-1],
                                              upper=x[i])$value)
      }
      return(y)
    }
    F<-function(x,a,b,omega){
      d<-matrix(0,nrow=length(x)-1,ncol=length(a))
      for(k in 1:length(a)){
        d[,k]<-diff(pgamma(x,shape=a[k],scale=b[k]))
      }
      return(omega%*%t(d))
    }
    Fi<-function(x,a,b) diff(pgamma(x,shape=a,scale=b))
    if (initial==FALSE){
      for(k in 1:K){
        omega.matrix[1,k]<-sum(clust$cluster==k)/sum(f)
        group.x<-sort(xmid[clust$cluster==k])
        sk<-mean((group.x-mean(group.x))^3)/(sqrt(var(group.x)))^3
        alpha.matrix[1,k]<-uniroot(function(u) trigamma(u)-var(log(group.x)[is.finite(log(group.x))]),c(0,1000000))$root
        beta.matrix[1,k]<-mean(group.x)/alpha.matrix[1,k]
      }
    }
    if(initial==TRUE){
      omega.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
    }
    while (eps>cri && h<=N){
      for (k in 1:K){
        ss<-omega.matrix[h-1,k]*Fi(r,alpha.matrix[h-1,k],beta.matrix[h-1,k])/
          F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,])
        omega.matrix[h,k]<-1/n*sum(f*ss,na.rm=TRUE)
        alpha.root<-function(w) {-digamma(w)*
            sum(f*Fi(r,alpha.matrix[h-1,k],beta.matrix[h-1,k])/F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,]),na.rm=TRUE)+
            sum(f*intalpha(r,alpha.matrix[h-1,k],beta.matrix[h-1,k])/F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,]))
        }
        alpha.matrix[h,k]<-uniroot(alpha.root,lower=.01,upper=30000)$root
        rrb<-r/beta.matrix[h-1,k]
        numb<-beta.matrix[h-1,k]*omega.matrix[h-1,k]*sum(f*Fi(rrb,alpha.matrix[h,k]+1,1)/
                                                           F(r,alpha.matrix[h-1,],beta.matrix[h-1,],omega.matrix[h-1,]))
        beta.matrix[h,k]<-numb/(n*omega.matrix[h,k])
      }
      eps<-sum(c(abs(omega.matrix[h-1,]-omega.matrix[h,])),c(abs(alpha.matrix[h-1,]-alpha.matrix[h,])),c(abs(beta.matrix[h-1,]-beta.matrix[h,])))
      h<-h+1
    }
    omega<-round(omega.matrix[h-1,],digits=7)
    omega[K]<-1-(sum(omega)-omega[K])
    alpha<-alpha.matrix[h-1,]
    beta<-beta.matrix[h-1,]
  }
  if(family=="skew-normal"){
    Ix<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      g1<-function(u){
        aa=u;bb=-1/2*(u-a)^2/b^2;cc=(pnorm(l*(u-a)/(b)))
        return(2/(sqrt(2*pi)*b)*aa*exp(bb)*(cc))
      }
      for (i in 2:length(x)){
        out[i-1]<-integrate(g1,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
      }
      return(out)
    }
    Imx2<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      g2<-function(u){
        aa=(u-a)^2;bb=-1/2*(u-a)^2/b^2;cc=(pnorm(l*(u-a)/(b)))
        return(2/(sqrt(2*pi)*b)*aa*exp(bb)*(cc))
      }
      for (i in 2:length(x)){
        out[i-1]<-integrate(g2,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
      }
      return(out)
    }
    Iux<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      g31<-function(u){
        aa=u;bb=-1/2*(u-a)^2/b^2*(1+l^2)
        return(aa*exp(bb))
      }
      g32<-function(u){
        aa=u*(u-a);bb=-1/2*(u-a)^2/b^2;cc=pnorm(l*(u-a)/b)
        return(aa*exp(bb)*cc)
      }

      for (i in 2:length(x)){
        out1<-1/(b*pi*sqrt((1+l^2)))*integrate(g31,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        out2<-2*l/(b^2*sqrt(2*pi*(1+l^2)))*integrate(g32,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        out[i-1]<-out1+out2
      }
      return(out)
    }
    Iumx<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      g81<-function(u){
        aa=u-a;bb=-1/2*(u-a)^2/b^2*(1+l^2)
        return(aa*exp(bb))
      }
      g82<-function(u){
        aa=(u-a)^2;bb=-1/2*(u-a)^2/b^2;cc=pnorm(l*(u-a)/b)
        return(aa*exp(bb)*cc)
      }

      for (i in 2:length(x)){
        out1<-1/(b*pi*sqrt((1+l^2)))*integrate(g81,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        out2<-2*l/(b^2*sqrt(2*pi*(1+l^2)))*integrate(g82,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        out[i-1]<-out1+out2
      }
      return(out)
    }
    Iu<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      for (i in 2:length(x)){
        g41<-function(u){
          bb=-1/2*((u-a)/b)^2*(1+l^2);
          return(1/(b*pi*sqrt(1+l^2))*exp(bb))
        }
        g42<-function(u){
          aa=(u-a);bb=-1/2*((u-a)/b)^2;cc=pnorm(l*(u-a)/(b))
          return(2*l/(b^2*sqrt(2*pi*(1+l^2)))*aa*exp(bb)*cc)
        }
        out1<-integrate(g41,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        out2<-integrate(g42,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        out[i-1]<-out1+out2
      }
      return(out)
    }
    Iu2<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      for (i in 2:length(x)){
        g51<-function(u){
          aa=u^2;bb=-1/2*u^2;cc=pnorm(sqrt(1+l^2)*(x[i]-a)/b-l*u)
          return(aa*exp(bb)*cc)
        }
        g52<-function(u){
          aa=u^2;bb=-1/2*u^2;cc=pnorm(sqrt(1+l^2)*(x[i-1]-a)/b-l*u)
          return(aa*exp(bb)*cc)
        }
        out1<-integrate(g51,lower=0,upper=Inf,rel.tol=10e-12)$value
        out2<-integrate(g52,lower=0,upper=Inf,rel.tol=10e-12)$value
        out[i-1]<-2/sqrt(2*pi)*(out1-out2)
      }
      return(out)
    }
    FiM<-function(x,a,b,l){
      out<-rep(0,(length(x)-1))
      for (i in 2:length(x)){
        g6<-function(u){
          bb=-1/2*(u-a)^2/b^2;cc=pnorm(l/b*(u-a))
          return(exp(bb)*cc)
        }
        out[i-1]<-2/(sqrt(2*pi*b^2))*integrate(g6,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
      }
      return(out)
    }
    FM<-function(x,a1,b1,l1,o1){
      d<-matrix(0,nrow=(length(x)-1),ncol=length(a1))
      for(k in 1:length(a1)){
        a<-a1[k];b<-b1[k];l<-l1[k];
        for (i in 2:length(x)){
          g7<-function(u){
            bb<--1/2*(u-a)^2/b^2;cc=pnorm(l/b*(u-a))
            return(exp(bb)*cc)
          }
          d[i-1,k]<-o1[k]*2/(sqrt(2*pi*b^2))*integrate(g7,lower=x[i-1],upper=x[i],rel.tol=10e-12)$value
        }
      }
      return(apply(d,1,sum))
    }
    if (initial==FALSE){
      for(k in 1:K){
        omega.matrix[1,k]<-sum(clust$cluster==k)/n
        y<-sort(xmid[clust$cluster==k])
        sk1<-mean((y-mean(y))^3)/(sqrt(var(y)))^3
        sk<-ifelse(abs(sk1)<0.99,sk1,0.99*sign(sk1))
        lambda.matrix[1,k]<-suppressWarnings(uniroot(function(u)sk-sqrt(2)*(4-pi)*u^3/(pi+(pi-2)*u^2)^(3/2),c(-100000000,100000000))$root)
        beta.matrix[1,k]<-sqrt(var(y)/(1-2*lambda.matrix[1,k]^2/(pi*(1+lambda.matrix[1,k]^2))))
        alpha.matrix[1,k]<-mean(y)-beta.matrix[1,k]*sqrt(2/pi)*lambda.matrix[1,k]/sqrt(1+lambda.matrix[1,k]^2)
      }
    }
    if(initial==TRUE){
      omega.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
    }
    while (eps>cri && h<=N){
      for (k in 1:K){
        IX<-Ix(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        IUX<-Iux(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        IU<-Iu(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        IMX2<-Imx2(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        IU2<-Iu2(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        IUMX<-Iumx(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        FIM<-FiM(r,alpha.matrix[h-1,k],beta.matrix[h-1,k],lambda.matrix[h-1,k])
        FT<-FM(r,alpha.matrix[h-1,],beta.matrix[h-1,],lambda.matrix[h-1,],omega.matrix[h-1,])
        ss<-omega.matrix[h-1,k]*FIM/FT
        omega.matrix[h,k]<-sum(f*ss,na.rm=TRUE)/n
        numa<-sum(f*IX/FT,na.rm=TRUE)-lambda.matrix[h-1,k]*beta.matrix[h-1,k]/sqrt(1+lambda.matrix[h-1,k]^2)*sum(f*IU/FT,na.rm=TRUE)
        alpha.matrix[h,k]<-omega.matrix[h-1,k]*numa/sum(f*ss,na.rm=TRUE)
        fbeta<-function(w){
          omega.matrix[h-1,k]*n*w^2+
            omega.matrix[h-1,k]*w*lambda.matrix[h-1,k]*sqrt(1+lambda.matrix[h-1,k]^2)*sum(f*IUX/FT,na.rm=TRUE)-
            omega.matrix[h-1,k]*w*lambda.matrix[h-1,k]*sqrt(1+lambda.matrix[h-1,k]^2)*alpha.matrix[h,k]*sum(f*IU/FT,na.rm=TRUE)-
            omega.matrix[h-1,k]*(1+lambda.matrix[h-1,k]^2)*sum(f*IMX2/FT,na.rm=TRUE)
        }
        beta.matrix[h,k]<-uniroot(fbeta,c(0.01,1000))$root
        flambda<-function(w){
          omega.matrix[h-1,k]*n*w/(1+w^2)-
            omega.matrix[h-1,k]*w*sum(f*IU2/FT,na.rm=TRUE)-
            omega.matrix[h-1,k]*w*sum(f*IMX2/FT,na.rm=TRUE)/beta.matrix[h-1,k]^2+
            omega.matrix[h-1,k]*(sqrt(1+w^2)/beta.matrix[h-1,k]+w^2/(sqrt(1+w^2)*beta.matrix[h-1,k]))*sum(f*IUMX/FT,na.rm=TRUE)
        }
        lambda.matrix[h,k]<-uniroot(flambda,c(-200000,200000))$root
      }
      eps<-sum(c(abs(omega.matrix[h-1,]-omega.matrix[h,])),c(abs(alpha.matrix[h-1,]-alpha.matrix[h,])),c(abs(beta.matrix[h-1,]-beta.matrix[h,])))
      h<-h+1
    }
    omega<-round(omega.matrix[h-1,],digits=7)
    omega[K]<-1-(sum(omega)-omega[K])
    alpha<-alpha.matrix[h-1,]
    beta<-beta.matrix[h-1,]
    lambda<-lambda.matrix[h-1,]
  }
  if(family=="skew-normal"){
    n.p<-4*K-1
    E<-n*FT
    prob<-E/sum(E)
    F1<-FM(c(-Inf,r[1]),alpha,beta,lambda,omega)[1]
    CS<-cumsum(FM(r,alpha,beta,lambda,omega))
    pg<-c(F1,CS+F1)
    pgks<-CS+F1
    out1<-cbind(omega,alpha,beta,lambda)
    colnames(out1)<-c("weight","alpha","beta","lambda")
  }
  else
  {
    n.p<-3*K-1
    E<-n*F(r,alpha,beta,omega)
    prob<-E/sum(E)
    F1<-F(c(-Inf,r[1]),alpha,beta,omega)[1]
    CS<-cumsum(F(r,alpha,beta,omega))
    pg<-c(F1,CS+F1)
    pgks<-CS+F1
    out1<-cbind(omega,alpha,beta)
    colnames(out1)<-c("weight","alpha","beta")
  }
  F.hat<-c(0,cumsum(f[-(m-1)])/(n+1))
  for(i in 2:(m-1)){
    anderson[i]<-F.hat[i]^2*log(pg[i+1]/pg[i]*(1-pg[i])/(1-pg[i+1]))+2*F.hat[i]*log((1-pg[i+1])/(1-pg[i]))
    cramer[i]<-F.hat[i]^2*(pg[i+1]-pg[i])-F.hat[i]*(pg[i+1]^2-pg[i]^2)+1/3*(pg[i+1]^3-pg[i]^3)
  }
  frn<-1
  frd<-0
  Fminus.hat<-c((cumsum(f)-frn)/(n+frd))
  Fplus.hat<-c(cumsum(f)/n)
  KS<-max(pgks-Fminus.hat,Fplus.hat-pgks)
  AD<-n*(sum(anderson)-(pg[m]-pg[1])-log((1-pg[m])/(1-pg[1])))
  CVM<-n*(sum(cramer)+1/3*(pg[2]^3-pg[1]^3))
  Chi.stat<-sum((f-E)^2/E)
  log.likelihood<-dmultinom(f,n,prob,log=TRUE)
  CAIC<--2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC<--2*log.likelihood + 2*n.p
  BIC<--2*log.likelihood + n.p*log(n)
  HQIC<--2*log.likelihood + 2*log(log(n))*n.p
  out2<-cbind(AIC, CAIC, BIC, HQIC, AD, Chi.stat, CVM, KS, log.likelihood)
  colnames(out2)<-c("AIC","CAIC","BIC","HQIC","AD","Chi-square","CVM","KS","log.likelihood")
  return(list("estimate"=out1,"measures"=out2))
}
fitWeibull<-function(data, location, method, starts){
  if(location==FALSE){
    if(method!="mlm" & method!="ml" &
       method!="wmle" & method!="rank" &
       method!="greg2" & method!="lm" &
       method!="greg1" &  method!="wreg" &
       method!="moment" & method!="ustat" &
       method!="pm" & method!="reg")
    {stop ("method class misspelled. Please check the manual for guidelines.") }}
  if(location==TRUE){
    if(method!="moment" & method!="ml" &
       method!="mm1" & method!="mm2" & method!="mm3"&
       method!="mml1" & method!="mml2" &
       method!="mml3" &  method!="mml4" & method!="wml"&
       method!="tlm" & method!="mps")
    {stop ("Baseline method misspelled. Please check the manual for guidelines.") }}
  n<-length(data)

  if(location==TRUE){
    if(method=="mps"){
      sort.x<-sort(data)
      pdf<-function(x,par){dweibull(x-par[3],par[1],par[2])}
      cdf<-function(x,par){pweibull(x-par[3],par[1],par[2])}
      f1<-function(x,par){
        z<-c()
        z<-diff(c(0,cdf(x,par),1))
        for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf(x[j],par)}}
        z[z<1e-16]<-1e-16
        -sum(log(z))}
      out<-suppressWarnings(optim(c(starts[1],starts[2],starts[3]),fn=f1,x=sort.x,method="Nelder-Mead")$par)
    }

    if(method=="wml"){
      n<-length(data)
      f2<-function(x,par){sum(-dweibull(x-par[3],par[1],scale=par[2],log=TRUE))}
      out<-suppressWarnings(optim(c(starts[1],starts[2],starts[3]),fn=f2,x=data,method="Nelder-Mead")$par)
      if (n>100){
        out<-c(out[1],out[2],out[3])
      }
      else
      {
        J1<-c(
          0.69340,0.83913,0.89139,0.91786,0.93423,0.94503	,
          0.95278,0.95865,0.96328,0.96683,0.96984,0.97232	,
          0.97443,0.97634,0.97786,0.97921,0.98046,0.98156	,
          0.98247,0.98336,0.98415,0.98487,0.98551,0.98611	,
          0.98672,0.98718,0.98770,0.98811,0.98853,0.98890	,
          0.98931,0.98962,0.98994,0.99021,0.99045,0.99076	,
          0.99099,0.99129,0.99148,0.99167,0.99191,0.99212	,
          0.99225,0.99245,0.99261,0.99277,0.99293,0.99309	,
          0.99319,0.99335,0.99345,0.99360,0.99370,0.99385	,
          0.99394,0.99407,0.99417,0.99429,0.99435,0.99447	,
          0.99455,0.99465,0.99472,0.99479,0.99489,0.99496	,
          0.99503,0.99510,0.99515,0.99525,0.99533,0.99536	,
          0.99545,0.99553,0.99556,0.99562,0.99569,0.99575	,
          0.99578,0.99583,0.99587,0.99595,0.99601,0.99605	,
          0.99612,0.99612,0.99617,0.99623,0.99627,0.99630	,
          0.99634,0.99638,0.99640,0.99645,0.99649,0.99655	,
          0.99655,0.99659,0.99666,0.996668)
        J2<-c(
          0.00000,0.27466,0.51747,0.63821,0.71001,0.75757	,
          0.79182,0.81712,0.83714,0.85306,0.86623,0.87744	,
          0.88628,0.89445,0.90121,0.90750,0.91274,0.91748	,
          0.92194,0.92570,0.92929,0.93244,0.93527,0.93794	,
          0.94037,0.94275,0.94468,0.94677,0.94847,0.95028	,
          0.95186,0.95340,0.95468,0.95618,0.95734,0.95856	,
          0.95942,0.96057,0.96166,0.96270,0.96350,0.96441	,
          0.96514,0.96589,0.96674,0.96734,0.96811,0.96873	,
          0.96944,0.96999,0.97061,0.97107,0.97168,0.97228	,
          0.97273,0.97321,0.97361,0.97419,0.97463,0.97510	,
          0.97537,0.97575,0.97609,0.97652,0.97678,0.97716	,
          0.97760,0.97787,0.97813,0.97848,0.97881,0.97911	,
          0.97941,0.97980,0.97999,0.98026,0.98048,0.98070	,
          0.98107,0.98111,0.98144,0.98159,0.98189,0.98204	,
          0.98230,0.98250,0.98274,0.98290,0.98315,0.98327	,
          0.98343,0.98363,0.98380,0.98401,0.98413,0.98424	,
          0.98447,0.98463,0.98480,0.985002)
        b1<-c(0.999419	,2.40369,3.83844,5.25078,6.70944,8.14346,9.52698,11.02088,12.47313,13.87016,
              15.25009,16.73677,18.15034,19.54663,21.09361,22.53010,23.82014,25.29579,26.70240,28.11381,
              29.52738,31.02171,32.47831,33.85212,35.36131,36.64201,38.17227,39.62372,41.05824,42.50542,
              43.87207,45.40919,46.78731,48.23230,49.78447,50.94642,52.43192,54.04720,55.32917,56.86716,
              58.27290,59.92009,61.20656,62.47552,64.04051,65.26954,66.93428,68.14113,69.67673,71.05012,
              72.55710,73.90928,75.29166,76.78865,78.26389,79.60301,80.93184,82.80045,84.22411,85.27033,
              87.11295,88.35019,89.55595,91.23598,92.77666,94.33027,95.34514,97.14380,98.32372,99.78025,
              101.13203,102.51998,103.72630,105.34980,107.13926,108.48973,109.70111,111.05117,112.15830,
              114.01928,116.01311,117.51730,118.10279,120.22575,121.04122,122.67439,124.08054,126.01121,
              126.57793,128.66007,129.70088,131.44117,132.91747,134.11285,135.60132,137.42157,138.36162,
              139.50138,141.14033,142.738630)
        b2<-c(1.000564	,2.38668,3.76955,5.14519,6.52922,7.91658,9.29997,10.70300,12.05786,13.47720,
              14.85522,16.19320,17.63708,18.99150,20.32042,21.75881,23.17091,24.49732,25.86710,27.24707,
              28.51663,29.99103,31.31774,32.79055,34.13763,35.46776,36.86314,38.33804,39.74500,41.00351,
              42.36087,43.89888,45.03703,46.52529,48.17510,49.35581,50.83383,52.27589,53.66510,54.95379,
              56.13727,57.46110,58.91282,60.48241,61.45341,62.95862,64.57392,66.21517,67.14057,68.77071,
              70.00901,71.30213,72.71345,73.97832,75.52757,77.31261,78.44624,80.07955,81.11361,82.43356,
              83.70449,85.34948,86.80341,87.79154,89.46271,90.42891,92.49190,93.99860,95.02369,96.02400,
              97.59045,98.95257,100.66270,101.84892,103.20334,104.31768,106.11673,107.25136,109.02587,
              110.22441,111.37984,112.93382,114.30203,115.74197,117.23615,118.55862,120.18218,121.10918,
              122.75864,123.83220,125.45253,126.73449,127.77614,129.58122,130.76601,132.21417,133.36799,
              134.70337,136.58863,137.424293)
        b3<-c(1.000195	,2.30542,3.60065,4.86645,6.11014,7.41961,8.64374,9.91321,11.13558,12.38813,
              13.66403,14.89621,16.11206,17.43501,18.63122,19.92408,21.15430,22.32046,23.74534,24.90473,
              26.17291,27.43755,28.60067,29.80842,31.30604,32.37293,33.57610,34.84576,36.07813,37.21921,
              38.56902,39.66478,40.93755,42.38144,43.50831,44.77444,46.02000,47.18365,48.31325,49.79495,
              50.88252,52.17011,53.48956,54.61035,55.81674,56.97552,58.34210,59.69581,60.88914,62.21924,
              63.34458,64.74456,65.96320,67.14392,68.22133,69.69377,70.90223,72.04439,73.47661,74.56693,
              75.62081,76.64199,78.18101,79.63297,80.76720,81.92504,83.40543,84.75744,86.36410,87.09991,
              88.25168,89.28729,90.73631,91.55075,93.13311,94.06731,95.56108,97.17017,98.18341,99.46701,
              100.59728,101.39081,103.29274,104.43916,105.17867,106.71837,108.18658,109.00555,110.69544,
              112.08485,113.41549,114.25492,115.64680,116.76379,118.07531,119.60246,120.51164,121.68645,
              122.16554,123.88978)
        b4<-c(0.999971	,2.21977,3.34996,4.47383,5.51848,6.57966,7.62637,8.61047,9.65870,10.68593,
              11.73894,12.61740,13.69749,14.63440,15.65282,16.58458,17.59279,18.62360,19.59686,20.51252,
              21.43462,22.48547,23.36335,24.38346,25.33662,26.28502,27.28960,28.17728,29.17154,30.04091,
              31.17642,32.08354,32.95523,34.02376,34.88262,35.77040,36.84907,37.87764,38.65829,39.63878,
              40.49853,41.40171,42.39290,43.19569,44.19390,45.18901,46.31572,47.24058,48.12679,49.12404,
              49.82614,51.03764,51.80719,52.62621,53.63460,54.63599,55.50899,56.29713,57.35872,58.30508,
              59.02388,60.22414,60.81050,62.07104,62.98124,63.94690,64.66725,65.69104,66.46643,67.69863,
              68.61025,69.42285,70.38871,71.61254,72.21700,73.09851,73.75254,75.05505,75.92815,76.74237,
              77.31276,78.48023,79.58264,80.90273,81.15164,82.18862,83.35053,84.64462,85.60832,86.11300,
              87.04442,88.11497,88.92334,89.77340,90.63447,91.88174,92.48561,93.37355,94.23893,95.463713)
        b5<-c(0.996472	,2.10278,3.06924,3.96978,4.81562,5.62821,6.43151,7.16495,7.89587,8.60808,
              9.36662,10.05582,10.74438,11.43316,12.03157,12.73752,13.43477,14.08647,14.61323,15.34230,
              15.97495,16.53045,17.13992,17.81857,18.40829,19.02063,19.59926,20.22399,20.78133,21.34597,
              21.95380,22.50925,23.09105,23.81801,24.38717,24.87731,25.49792,26.04568,26.54291,27.31364,
              27.82623,28.20748,28.94449,29.41787,29.99535,30.59824,31.05490,31.62618,32.15814,32.67728,
              33.18654,33.83011,34.40781,34.87585,35.55083,35.99780,36.47847,36.99387,37.53880,37.81647,
              38.76525,39.03579,39.69058,40.28051,40.67256,41.41254,41.78470,42.14900,42.82274,43.51694,
              43.94128,44.30513,44.89677,45.40447,46.02973,46.27556,46.81203,47.41376,47.93676,48.42347,
              48.95905,49.45992,49.91452,50.71609,51.13429,51.37643,52.03785,52.46086,52.97077,53.85307,
              54.29694,54.71287,54.85817,55.59169,55.90096,56.43494,56.77928,57.75950,58.12053,58.549002)
        b6<-c(1.000196	,1.99030,2.81205,3.51868,4.16587,4.74759,5.32893,5.86577,6.38462,6.83923,
              7.34450,7.76121,8.23533,8.61970,9.06532,9.49624,9.88269,10.29285,10.63846,11.02887,11.36246,
              11.78019,12.10889,12.47425,12.77896,13.07602,13.43002,13.78611,14.09811,14.44406,14.79079,
              15.01301,15.33809,15.70616,15.99263,16.30830,16.61772,16.90987,17.20681,17.44492,17.71759,
              18.05578,18.29872,18.58609,18.93141,19.22764,19.40701,19.65032,19.94359,20.24321,20.52529,
              20.77194,21.01384,21.27935,21.52844,21.76766,22.02346,22.31122,22.56081,22.79343,23.18409,
              23.35571,23.54428,23.87096,24.07365,24.22840,24.61390,24.74033,25.06939,25.27360,25.47626,
              25.65355,26.02591,26.16731,26.40300,26.53790,26.89221,27.20662,27.25878,27.56878,27.67879,
              27.95406,28.34934,28.64513,28.63284,28.95600,29.04384,29.33999,29.55270,29.86158,29.94223,
              30.28504,30.42563,30.60974,30.96228,31.00309,31.32467,31.52457,31.60065,31.849782)
        b7<-c(0.999314	,1.88719,2.57599,3.13205,3.63167,4.06707,4.47511,4.84220,5.18637,5.50847,
              5.80998,6.10146,6.37533,6.61851,6.91297,7.13995,7.38415,7.63689,7.84061,8.07189,8.27679,
              8.48186,8.68562,8.86874,9.06159,9.23659,9.40549,9.58358,9.74598,9.95594,10.11137,10.25876,
              10.46182,10.62677,10.73090,10.88898,11.08312,11.17984,11.37647,11.50089,11.63824,11.82247,
              11.94986,12.07388,12.17139,12.32177,12.46272,12.56700,12.71522,12.87449,12.96957,13.09991,
              13.23694,13.34446,13.45254,13.57921,13.70137,13.83660,13.95533,14.04745,14.20660,14.28288,
              14.38101,14.48047,14.64452,14.72823,14.84201,14.94054,15.04008,15.19630,15.27327,15.32566,
              15.44946,15.57988,15.68466,15.77853,15.87207,15.98169,16.04324,16.16120,16.31488,16.37536,
              16.46770,16.51017,16.65114,16.81005,16.88910,16.94576,17.04232,17.10823,17.21106,17.27973,
              17.41461,17.45630,17.55699,17.64537,17.73035,17.77916,17.88631,17.960293)
        b8<-c(0.999313	,1.79956,2.38519,2.83185,3.20546,3.54772,3.82961,4.08761,4.33165,4.54456,
              4.74419,4.93058,5.11529,5.28766,5.45189,5.59469,5.75664,5.88953,6.02516,6.14841,6.27954,
              6.38477,6.52237,6.63764,6.73031,6.85435,6.94532,7.05680,7.14425,7.24954,7.34991,7.42338,
              7.51603,7.59433,7.67782,7.77718,7.80458,7.91923,8.01339,8.07692,8.16227,8.20089,8.28489,
              8.35928,8.44692,8.51275,8.58706,8.61122,8.69214,8.77283,8.83065,8.90949,8.96210,9.03180,
              9.09699,9.12670,9.20488,9.25399,9.31376,9.37610,9.42642,9.52052,9.54833,9.58419,9.64671,
              9.69201,9.76021,9.79894,9.85650,9.91076,9.97786,10.01790,10.05014,10.09650,10.15798,
              10.16661,10.24708,10.26003,10.32319,10.38774,10.43855,10.47729,10.51818,10.53054,
              10.58725,10.66083,10.65822,10.70976,10.78579,10.79213,10.86947,10.89583,10.96527,
              10.96435,11.01067,11.06136,11.09744,11.14734,11.18324,11.195746)
        b9<-c(1.001712	,1.73573,2.21682,2.57481,2.88535,3.12157,3.33993,3.53106,3.71570,3.84315,
              3.99344,4.13971,4.24437,4.36185,4.46240,4.56669,4.66580,4.75844,4.83760,4.90944,4.99703,
              5.07223,5.14526,5.21290,5.27959,5.34777,5.39952,5.45864,5.51329,5.57709,5.63713,5.68813,
              5.73984,5.79700,5.81976,5.89586,5.94014,5.97187,6.03161,6.06019,6.11181,6.15849,6.18791,
              6.23379,6.26967,6.29588,6.34017,6.37678,6.40855,6.45794,6.48963,6.52739,6.55929,6.58798,
              6.62154,6.65562,6.67756,6.71431,6.73369,6.77514,6.80959,6.82230,6.87419,6.88817,6.91283,
              6.94904,6.96675,7.01370,7.01562,7.04580,7.08196,7.10396,7.13121,7.16878,7.17286,7.22021,
              7.22656,7.25879,7.25202,7.29646,7.31816,7.33465,7.34777,7.36925,7.41397,7.43417,7.45066,
              7.46395,7.50323,7.51835,7.55227,7.54456,7.59264,7.59888,7.61631,7.63189,7.64291,7.68615,
              7.69150,7.7125825)
        b10<-c(1.00265	,1.66817,2.08180,2.39358,2.62218,2.81960,2.98409,3.13332,3.24384,3.35977,
               3.45256,3.55606,3.63295,3.72190,3.78654,3.86331,3.91876,3.98394,4.04344,4.08871,4.13740,
               4.19805,4.23774,4.28808,4.32639,4.36624,4.40593,4.45164,4.46739,4.51017,4.55308,4.58706,
               4.61805,4.64381,4.67016,4.70130,4.72340,4.76899,4.80166,4.81243,4.84487,4.86756,4.88985,
               4.92066,4.94104,4.96407,4.99087,5.00618,5.02965,5.05341,5.06226,5.09132,5.10967,5.12648,
               5.15031,5.16878,5.18130,5.20751,5.22234,5.22948,5.26046,5.26621,5.28325,5.30394,5.31904,
               5.34608,5.35033,5.38039,5.38433,5.39325,5.41439,5.42920,5.44360,5.45734,5.47042,5.48325,
               5.49855,5.50126,5.52403,5.53312,5.54743,5.55859,5.57451,5.58296,5.59481,5.61407,5.62522,
               5.62962,5.63921,5.66400,5.66604,5.69557,5.69954,5.70874,5.71769,5.72711,5.73455,5.74685,
               5.74939,5.7714421)
        b11<-c(1.000197	,1.61344,1.97956,2.22541,2.42127,2.57163,2.70417,2.81655,2.91219,2.98803,
               3.06902,3.12649,3.19813,3.24815,3.30346,3.35616,3.39832,3.44758,3.48122,3.52209,3.54911,
               3.59082,3.62186,3.65593,3.67808,3.70352,3.73089,3.75838,3.78415,3.81034,3.83272,3.84684,
               3.87534,3.89669,3.91666,3.93976,3.95092,3.97091,3.98843,4.00718,4.01924,4.03588,4.05096,
               4.06685,4.08504,4.10408,4.10870,4.11633,4.13955,4.15089,4.16445,4.17888,4.18957,4.19971,
               4.20908,4.22368,4.23001,4.24837,4.25543,4.26828,4.28752,4.29263,4.30364,4.31207,4.31957,
               4.32598,4.34218,4.34858,4.36277,4.36729,4.37841,4.38078,4.39804,4.39978,4.40632,4.41474,
               4.42481,4.43908,4.44022,4.44976,4.45468,4.46237,4.47772,4.49004,4.48731,4.49771,4.49807,
               4.50733,4.51298,4.52723,4.52835,4.54027,4.54135,4.54683,4.56313,4.55861,4.57599,4.58016,
               4.57904,4.5869005)
        b12<-c(0.999312	,1.56029,1.88643,2.09415,2.25851,2.38369,2.49030,2.57168,2.64861,2.71148,
               2.76953,2.82061,2.86526,2.90268,2.94945,2.98525,3.01844,3.05361,3.07669,3.10772,3.13223,
               3.15606,3.17907,3.19922,3.22080,3.23969,3.25673,3.27343,3.29096,3.31207,3.32343,3.33990,
               3.35558,3.37200,3.38284,3.39527,3.41094,3.41954,3.43840,3.44257,3.45597,3.47008,3.48010,
               3.48848,3.49464,3.50503,3.51766,3.52506,3.53239,3.54654,3.55018,3.56004,3.57125,3.58092,
               3.58356,3.58993,3.60068,3.60830,3.61601,3.61997,3.63163,3.63543,3.63987,3.64541,3.65741,
               3.66243,3.66464,3.67068,3.67964,3.68723,3.69269,3.69409,3.70061,3.70701,3.71348,3.71737,
               3.72431,3.73085,3.73467,3.73879,3.74794,3.75085,3.75365,3.75515,3.76438,3.77149,3.77580,
               3.78001,3.77959,3.78500,3.78890,3.79458,3.80104,3.80354,3.80798,3.80914,3.81465,3.81594,
               3.82123,3.8250084)
        b13<-c(0.999313	,1.51778,1.80697,1.98980,2.12464,2.23317,2.31459,2.38403,2.44554,2.49448,
               2.53894,2.57812,2.61388,2.65054,2.68036,2.70190,2.73166,2.75485,2.77701,2.79617,2.81779,
               2.83043,2.85028,2.86811,2.88154,2.89779,2.91026,2.92427,2.93638,2.95140,2.96162,2.97118,
               2.97886,2.99050,2.99884,3.01249,3.01331,3.02944,3.03723,3.04413,3.05449,3.05844,3.06625,
               3.07059,3.08258,3.08952,3.09653,3.09816,3.10812,3.11155,3.12058,3.12733,3.13229,3.13701,
               3.14467,3.14406,3.15587,3.15833,3.16239,3.16828,3.17208,3.18068,3.18321,3.18667,3.19226,
               3.19473,3.20284,3.20315,3.20828,3.21315,3.21967,3.22017,3.22407,3.22656,3.23277,3.23099,
               3.23990,3.23688,3.24389,3.25053,3.25159,3.25718,3.26005,3.25720,3.26262,3.26786,3.26781,
               3.27247,3.27891,3.27668,3.28344,3.28500,3.28903,3.28838,3.29134,3.29815,3.29731,3.30107,
               3.30497,3.3055204)
        b14<-c(1.001713	,1.48848,1.73878,1.89472,2.01567,2.09900,2.17416,2.22847,2.28359,2.31701,
               2.35526,2.39503,2.41740,2.44484,2.46773,2.49109,2.51238,2.52868,2.54512,2.55909,2.57575,
               2.58831,2.60119,2.61427,2.62591,2.63724,2.64718,2.65355,2.66428,2.67596,2.68437,2.69213,
               2.69887,2.70995,2.71004,2.72381,2.72932,2.73472,2.74410,2.74597,2.75467,2.75957,2.76513,
               2.77036,2.77391,2.77616,2.78365,2.78965,2.79172,2.79880,2.80271,2.80699,2.81181,2.81431,
               2.81835,2.82040,2.82490,2.83092,2.83196,2.83588,2.84131,2.84187,2.84925,2.85001,2.85243,
               2.85475,2.85810,2.86311,2.86236,2.86747,2.86936,2.87457,2.87602,2.87821,2.88080,2.88789,
               2.88611,2.88867,2.88735,2.89287,2.89563,2.89558,2.89723,2.89824,2.90536,2.90638,2.90685,
               2.91055,2.91416,2.91501,2.91938,2.91795,2.92186,2.92236,2.92386,2.92359,2.92684,2.93199,
               2.93152,2.9332782)
        b15<-c(1.002684	,1.45286,1.67738,1.82583,1.92148,1.99847,2.05667,2.11003,2.14431,2.18155,
               2.20805,2.23815,2.25859,2.28551,2.30106,2.31711,2.33240,2.34922,2.36312,2.37130,2.38412,
               2.39723,2.40504,2.41618,2.42605,2.43434,2.44173,2.45149,2.45367,2.46150,2.47036,2.47739,
               2.48187,2.48746,2.49076,2.49788,2.50198,2.51130,2.51656,2.51677,2.52372,2.52675,2.53054,
               2.53795,2.53931,2.54372,2.54742,2.55069,2.55353,2.55871,2.55873,2.56371,2.56634,2.56828,
               2.57224,2.57637,2.57703,2.58180,2.58449,2.58482,2.58977,2.59014,2.59305,2.59503,2.59822,
               2.60274,2.60230,2.60639,2.60682,2.60701,2.61230,2.61341,2.61473,2.61618,2.61954,2.62167,
               2.62264,2.62271,2.62668,2.62719,2.62826,2.63088,2.63361,2.63420,2.63535,2.63868,2.63971,
               2.63986,2.64099,2.64436,2.64535,2.64736,2.64790,2.65040,2.65030,2.65144,2.65332,2.65572,
               2.65496,2.6588466)
        b16<-c(1.000195	,1.42565,1.63334,1.75403,1.84364,1.90792,1.96137,2.00369,2.03703,2.06426,
               2.09221,2.10932,2.13325,2.14861,2.16442,2.18045,2.19090,2.20592,2.21425,2.22736,2.23128,
               2.24454,2.25107,2.26147,2.26511,2.27020,2.27728,2.28358,2.29152,2.29763,2.30239,2.30478,
               2.31350,2.31749,2.31957,2.32532,2.32854,2.33431,2.33623,2.34149,2.34280,2.34667,2.35043,
               2.35319,2.35757,2.36006,2.36005,2.36183,2.36806,2.36975,2.37219,2.37347,2.37672,2.37782,
               2.37963,2.38343,2.38381,2.38870,2.38989,2.39152,2.39564,2.39561,2.39826,2.39965,2.39955,
               2.40124,2.40367,2.40641,2.40830,2.40829,2.41113,2.41026,2.41417,2.41366,2.41435,2.41671,
               2.41786,2.42100,2.42094,2.42187,2.42224,2.42440,2.42776,2.42973,2.42933,2.43058,2.43065,
               2.43264,2.43236,2.43565,2.43515,2.43828,2.43729,2.43796,2.44182,2.43971,2.44327,2.44423,
               2.44256,2.4448841)
        b17<-c(0.999316	,1.39520,1.58884,1.69756,1.77777,1.83419,1.87919,1.91211,1.94609,1.97020,
               1.98954,2.00865,2.02395,2.03676,2.05275,2.06614,2.07549,2.08814,2.09386,2.10375,2.10912,
               2.11789,2.12536,2.13029,2.13653,2.14184,2.14643,2.15128,2.15540,2.16301,2.16603,2.16813,
               2.17367,2.17757,2.17935,2.18271,2.18775,2.18939,2.19378,2.19542,2.19886,2.20310,2.20482,
               2.20601,2.20734,2.20954,2.21185,2.21545,2.21559,2.21959,2.22000,2.22303,2.22504,2.22951,
               2.22878,2.22916,2.23141,2.23429,2.23526,2.23662,2.23905,2.23931,2.24013,2.24183,2.24416,
               2.24574,2.24496,2.24572,2.24833,2.25030,2.25097,2.25159,2.25278,2.25394,2.25626,2.25687,
               2.25847,2.25944,2.25977,2.26056,2.26321,2.26374,2.26299,2.26408,2.26652,2.26704,2.26898,
               2.26871,2.26807,2.26977,2.27133,2.27163,2.27301,2.27288,2.27428,2.27441,2.27678,2.27579,
               2.27661,2.2776023)
        b18<-c(0.999313	,1.37223,1.54877,1.65001,1.71875,1.76990,1.80974,1.84091,1.86818,1.88823,
               1.90611,1.92170,1.93395,1.94842,1.96028,1.96703,1.97776,1.98605,1.99511,2.00150,2.00772,
               2.01230,2.01906,2.02443,2.02847,2.03382,2.03816,2.04314,2.04666,2.05122,2.05367,2.05624,
               2.05779,2.06270,2.06651,2.06952,2.06786,2.07488,2.07584,2.07784,2.08153,2.08258,2.08305,
               2.08480,2.08820,2.09072,2.09281,2.09251,2.09576,2.09627,2.09927,2.10052,2.10224,2.10248,
               2.10548,2.10462,2.10857,2.10878,2.10924,2.11142,2.11167,2.11398,2.11510,2.11564,2.11740,
               2.11776,2.11931,2.11987,2.12082,2.12211,2.12355,2.12398,2.12486,2.12568,2.12719,2.12588,
               2.12905,2.12711,2.13001,2.13198,2.13039,2.13307,2.13318,2.13170,2.13343,2.13610,2.13496,
               2.13529,2.13722,2.13587,2.13938,2.13903,2.14031,2.13912,2.13927,2.14206,2.14175,2.14234,
               2.14309,2.1427779)
        b19<-c(1.001715	,1.35564,1.51254,1.60281,1.66904,1.71175,1.75004,1.77546,1.79953,1.81512,
               1.83215,1.85034,1.85928,1.87013,1.88005,1.88851,1.89782,1.90343,1.91045,1.91605,1.92196,
               1.92624,1.93039,1.93545,1.93940,1.94283,1.94688,1.94942,1.95242,1.95765,1.95944,1.96259,
               1.96484,1.96891,1.96787,1.97418,1.97503,1.97681,1.97912,1.97968,1.98397,1.98371,1.98663,
               1.98883,1.98912,1.98953,1.99155,1.99317,1.99420,1.99595,1.99731,1.99875,2.00085,2.00093,
               2.00115,2.00163,2.00326,2.00599,2.00614,2.00778,2.00944,2.00922,2.01070,2.01068,2.01186,
               2.01139,2.01313,2.01403,2.01499,2.01653,2.01547,2.01876,2.01822,2.01837,2.01898,2.02143,
               2.02133,2.02170,2.02102,2.02211,2.02308,2.02348,2.02350,2.02338,2.02639,2.02741,2.02637,
               2.02754,2.02865,2.02834,2.02980,2.02928,2.03059,2.03005,2.03030,2.03020,2.03109,2.03263,
               2.03288,2.0328423)
        b20<-c(1.002683	,1.33465,1.48104,1.56949,1.62283,1.66575,1.69624,1.72449,1.74068,1.75860,
               1.77025,1.78439,1.79278,1.80531,1.81265,1.81847,1.82521,1.83271,1.83883,1.84013,1.84743,
               1.85190,1.85403,1.86076,1.86358,1.86647,1.86970,1.87262,1.87384,1.87651,1.88001,1.88291,
               1.88429,1.88639,1.88699,1.89041,1.89096,1.89490,1.89730,1.89769,1.89951,1.89936,1.90104,
               1.90504,1.90414,1.90589,1.90773,1.90910,1.90973,1.91228,1.91136,1.91246,1.91348,1.91413,
               1.91522,1.91706,1.91697,1.91974,1.91888,1.92015,1.92198,1.92155,1.92292,1.92266,1.92377,
               1.92590,1.92494,1.92682,1.92690,1.92651,1.92886,1.92928,1.92947,1.92954,1.93030,1.93104,
               1.93215,1.93171,1.93345,1.93265,1.93228,1.93466,1.93539,1.93505,1.93569,1.93611,1.93637,
               1.93608,1.93648,1.93804,1.93867,1.93872,1.93949,1.94015,1.93916,1.94022,1.94035,1.94135,
               1.94054,1.9423633)
        b21<-c(0.998716	,1.31806,1.45690,1.53384,1.58647,1.62302,1.64862,1.67193,1.68934,1.70560,
               1.71645,1.72978,1.73634,1.74741,1.75164,1.75894,1.76475,1.77155,1.77680,1.78019,1.78501,
               1.78802,1.79014,1.79417,1.79842,1.80060,1.80332,1.80449,1.80894,1.81101,1.81343,1.81529,
               1.81717,1.81711,1.82043,1.82190,1.82385,1.82509,1.82575,1.82732,1.82905,1.83108,1.83214,
               1.83206,1.83338,1.83439,1.83594,1.83782,1.83722,1.83772,1.83925,1.84020,1.84076,1.84151,
               1.84145,1.84425,1.84411,1.84428,1.84592,1.84594,1.84777,1.84785,1.84790,1.84805,1.84842,
               1.85034,1.85049,1.85033,1.85189,1.85167,1.85161,1.85393,1.85364,1.85319,1.85439,1.85528,
               1.85588,1.85687,1.85642,1.85649,1.85679,1.85874,1.85896,1.85901,1.85885,1.85905,1.86045,
               1.85944,1.86174,1.86123,1.86172,1.86193,1.86170,1.86235,1.86159,1.86381,1.86386,1.86257,
               1.86412,1.8628100)
        b22<-c(1.002518	,1.30563,1.43036,1.49915,1.54779,1.58266,1.60921,1.63084,1.64584,1.65889,
               1.67076,1.67899,1.68800,1.69629,1.70114,1.70776,1.71256,1.71771,1.72184,1.72596,1.72947,
               1.73250,1.73487,1.73732,1.74113,1.74359,1.74653,1.74625,1.75004,1.75123,1.75451,1.75456,
               1.75740,1.75948,1.75938,1.76272,1.76378,1.76394,1.76590,1.76688,1.76826,1.76912,1.76872,
               1.77031,1.77182,1.77318,1.77367,1.77481,1.77567,1.77566,1.77773,1.77889,1.77852,1.77897,
               1.77885,1.78116,1.78136,1.78217,1.78333,1.78290,1.78277,1.78304,1.78508,1.78446,1.78526,
               1.78511,1.78602,1.78743,1.78737,1.78767,1.78864,1.78899,1.78875,1.79012,1.78997,1.79076,
               1.79024,1.79191,1.79166,1.79338,1.79185,1.79293,1.79285,1.79348,1.79383,1.79442,1.79487,
               1.79412,1.79482,1.79477,1.79565,1.79586,1.79610,1.79687,1.79681,1.79688,1.79694,1.79725,
               1.79799,1.7978944)
        b23<-c(1.000092	,1.28622,1.41019,1.47509,1.52097,1.54847,1.57159,1.59185,1.60609,1.61693,
               1.62803,1.63747,1.64578,1.65050,1.65677,1.66166,1.66605,1.67042,1.67489,1.67806,1.68051,
               1.68342,1.68732,1.68999,1.69208,1.69470,1.69635,1.69702,1.70033,1.69990,1.70279,1.70529,
               1.70533,1.70688,1.70852,1.70998,1.71051,1.71354,1.71284,1.71370,1.71464,1.71591,1.71657,
               1.71816,1.71954,1.71829,1.71950,1.72073,1.72234,1.72189,1.72314,1.72356,1.72542,1.72455,
               1.72569,1.72636,1.72682,1.72761,1.72748,1.72917,1.72876,1.72961,1.72993,1.72999,1.73079,
               1.73129,1.73177,1.73126,1.73211,1.73282,1.73234,1.73358,1.73360,1.73375,1.73391,1.73542,
               1.73364,1.73552,1.73602,1.73559,1.73645,1.73645,1.73578,1.73671,1.73845,1.73727,1.73755,
               1.73808,1.73754,1.73907,1.73906,1.73945,1.73891,1.73875,1.74051,1.73944,1.74005,1.74087,
               1.74011,1.7404808)
        b24<-c(1.005307	,1.27293,1.38824,1.45132,1.49050,1.51744,1.54128,1.55701,1.57211,1.58046,
               1.59316,1.59888,1.60569,1.61262,1.61814,1.62002,1.62696,1.62943,1.63290,1.63455,1.63818,
               1.64212,1.64418,1.64549,1.64831,1.65041,1.65254,1.65348,1.65542,1.65623,1.65920,1.65945,
               1.66181,1.66167,1.66286,1.66456,1.66532,1.66677,1.66798,1.66882,1.66980,1.67003,1.66989,
               1.67132,1.67253,1.67152,1.67241,1.67373,1.67520,1.67494,1.67569,1.67596,1.67627,1.67739,
               1.67798,1.67768,1.67902,1.68002,1.68049,1.67973,1.68098,1.68043,1.68240,1.68166,1.68201,
               1.68232,1.68261,1.68337,1.68299,1.68337,1.68465,1.68483,1.68466,1.68496,1.68592,1.68583,
               1.68616,1.68597,1.68612,1.68670,1.68766,1.68583,1.68688,1.68851,1.68835,1.68824,1.68890,
               1.68844,1.68853,1.68968,1.68951,1.68953,1.68953,1.68984,1.69040,1.69015,1.69092,1.69031,
               1.69077,1.6914009)
        b25<-c(0.998827	,1.26181,1.37058,1.42564,1.46532,1.49006,1.51291,1.52554,1.54017,1.54796,
               1.55784,1.56540,1.56990,1.57473,1.58066,1.58619,1.58932,1.59039,1.59569,1.59938,1.60101,
               1.60227,1.60597,1.60670,1.60860,1.60949,1.61231,1.61347,1.61571,1.61715,1.61755,1.61837,
               1.61991,1.62041,1.62259,1.62517,1.62419,1.62461,1.62626,1.62646,1.62802,1.62946,1.62885,
               1.62977,1.63102,1.63148,1.63147,1.63338,1.63239,1.63361,1.63359,1.63452,1.63335,1.63584,
               1.63544,1.63699,1.63633,1.63625,1.63783,1.63830,1.63832,1.63821,1.63869,1.63960,1.63951,
               1.64039,1.64058,1.64038,1.64116,1.64162,1.64132,1.64204,1.64224,1.64328,1.64267,1.64256,
               1.64314,1.64322,1.64314,1.64358,1.64424,1.64404,1.64407,1.64529,1.64438,1.64483,1.64497,
               1.64550,1.64574,1.64540,1.64643,1.64684,1.64635,1.64650,1.64669,1.64704,1.64722,1.64659,
               1.64774,1.6478516)
        b26<-c(1.000190	,1.25061,1.35309,1.40481,1.44133,1.46650,1.48551,1.50030,1.51070,1.51947,
               1.52951,1.53300,1.53958,1.54511,1.54941,1.55314,1.55573,1.55996,1.56171,1.56651,1.56610,
               1.56969,1.57143,1.57472,1.57455,1.57607,1.57698,1.57798,1.58096,1.58174,1.58213,1.58374,
               1.58522,1.58615,1.58666,1.58759,1.58794,1.58990,1.58936,1.59148,1.59081,1.59211,1.59283,
               1.59342,1.59516,1.59465,1.59417,1.59490,1.59643,1.59657,1.59662,1.59686,1.59709,1.59769,
               1.59831,1.59882,1.59879,1.60023,1.59998,1.60034,1.60082,1.60106,1.60166,1.60134,1.60171,
               1.60138,1.60201,1.60319,1.60295,1.60275,1.60360,1.60331,1.60420,1.60395,1.60338,1.60418,
               1.60450,1.60502,1.60546,1.60500,1.60551,1.60541,1.60634,1.60668,1.60666,1.60653,1.60639,
               1.60734,1.60668,1.60722,1.60782,1.60812,1.60788,1.60760,1.60843,1.60799,1.60900,1.60891,
               1.60825,1.6092222)
        b27<-c(0.999319	,1.23764,1.33766,1.38650,1.42239,1.44447,1.46130,1.47302,1.48588,1.49443,
               1.49903,1.50561,1.51054,1.51440,1.52058,1.52447,1.52776,1.53077,1.53278,1.53567,1.53607,
               1.53859,1.54034,1.54255,1.54392,1.54514,1.54615,1.54745,1.54906,1.55031,1.55157,1.55141,
               1.55297,1.55425,1.55374,1.55517,1.55637,1.55692,1.55771,1.55774,1.55940,1.55952,1.56009,
               1.56000,1.56094,1.56108,1.56103,1.56243,1.56252,1.56341,1.56290,1.56433,1.56442,1.56575,
               1.56534,1.56520,1.56556,1.56633,1.56683,1.56709,1.56747,1.56710,1.56724,1.56730,1.56849,
               1.56860,1.56802,1.56795,1.56902,1.56952,1.56941,1.56971,1.57019,1.57006,1.57025,1.57049,
               1.57110,1.57109,1.57126,1.57143,1.57202,1.57201,1.57167,1.57161,1.57280,1.57178,1.57270,
               1.57258,1.57253,1.57290,1.57333,1.57320,1.57385,1.57339,1.57364,1.57364,1.57428,1.57398,
               1.57406,1.5746339)
        b28<-c(0.999313	,1.22819,1.32080,1.37242,1.40114,1.42311,1.43942,1.45138,1.46219,1.46963,
               1.47686,1.48194,1.48590,1.49055,1.49496,1.49705,1.49962,1.50306,1.50581,1.50802,1.50991,
               1.51007,1.51270,1.51461,1.51540,1.51749,1.51868,1.52037,1.52048,1.52235,1.52217,1.52274,
               1.52369,1.52492,1.52696,1.52726,1.52656,1.52860,1.52880,1.52835,1.53024,1.53013,1.52984,
               1.53040,1.53099,1.53264,1.53313,1.53222,1.53393,1.53373,1.53427,1.53485,1.53445,1.53488,
               1.53582,1.53527,1.53627,1.53640,1.53673,1.53681,1.53661,1.53826,1.53789,1.53842,1.53806,
               1.53831,1.53933,1.53906,1.53902,1.53912,1.53938,1.53918,1.53984,1.53984,1.54019,1.53981,
               1.54091,1.54015,1.54072,1.54161,1.54050,1.54183,1.54141,1.54088,1.54141,1.54241,1.54194,
               1.54171,1.54240,1.54156,1.54303,1.54226,1.54328,1.54281,1.54240,1.54326,1.54292,1.54325,
               1.54371,1.5432129)
        b29<-c(1.001710	,1.22285,1.30683,1.35250,1.38446,1.40342,1.42081,1.43006,1.44042,1.44621,
               1.45306,1.46070,1.46321,1.46782,1.47119,1.47371,1.47734,1.47925,1.48162,1.48316,1.48564,
               1.48605,1.48752,1.48961,1.48964,1.49173,1.49239,1.49364,1.49493,1.49631,1.49685,1.49771,
               1.49806,1.50021,1.49906,1.50118,1.50127,1.50154,1.50250,1.50249,1.50406,1.50380,1.50486,
               1.50597,1.50551,1.50513,1.50581,1.50585,1.50620,1.50699,1.50729,1.50796,1.50791,1.50838,
               1.50766,1.50816,1.50844,1.50903,1.50957,1.51007,1.51054,1.51016,1.51064,1.51058,1.51110,
               1.51040,1.51114,1.51155,1.51142,1.51214,1.51126,1.51260,1.51268,1.51213,1.51247,1.51320,
               1.51374,1.51351,1.51254,1.51327,1.51314,1.51362,1.51330,1.51318,1.51405,1.51481,1.51472,
               1.51482,1.51506,1.51422,1.51531,1.51502,1.51526,1.51520,1.51479,1.51507,1.51539,1.51599,
               1.51561,1.5160695)
        b30<-c(1.002688	,1.21259,1.29450,1.34180,1.36680,1.38772,1.40188,1.41506,1.41982,1.42895,
               1.43225,1.43862,1.44095,1.44757,1.44901,1.45174,1.45424,1.45627,1.45925,1.45945,1.46152,
               1.46337,1.46345,1.46710,1.46747,1.46877,1.46924,1.47076,1.47115,1.47159,1.47311,1.47383,
               1.47418,1.47543,1.47515,1.47591,1.47594,1.47746,1.47881,1.47844,1.47894,1.47858,1.47940,
               1.48079,1.48041,1.48108,1.48148,1.48203,1.48174,1.48326,1.48222,1.48294,1.48294,1.48352,
               1.48327,1.48405,1.48360,1.48559,1.48434,1.48492,1.48609,1.48535,1.48622,1.48548,1.48627,
               1.48665,1.48617,1.48726,1.48673,1.48672,1.48779,1.48749,1.48768,1.48791,1.48780,1.48790,
               1.48846,1.48837,1.48901,1.48851,1.48829,1.48907,1.48966,1.48934,1.48927,1.48914,1.48968,
               1.48936,1.48917,1.49019,1.48996,1.49020,1.49061,1.49032,1.48967,1.49044,1.49037,1.49084,
               1.49042,1.4907879)
        J3<-matrix(rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,
                         b21,b22,b23,b24,b25,b26,b27,b28,b29,b30),ncol=100,nrow=30)
        hat<-matrix(0,ncol=2,nrow=11)
        n<-length(data)
        hat[1,]<-c(out[1],out[3])
        jj<-splinefun(seq(.1,3,.1),J3[,n])(hat[1,1])
        for (r in 1:10){
          f3<-function(x,par){(-sum(J2[n]/(n*par[1]))-sum(log(abs(x-par[2]))/n)+sum(log(abs(x-par[2]))*
                                                                                      (abs(x-par[2]))^par[1])/sum((abs(x-par[2]))^par[1]))^2-(sum(1/(n*(abs(x-par[2]))))*sum((abs(x-par[2]))^par[1])/
                                                                                                                                                sum((abs(x-par[2]))^(par[1]-1))-jj)^2}
          hat[r+1,]<-suppressWarnings(optim(c(abs(hat[1,1]),hat[1,2]),fn=f3,x=data,method="Nelder-Mead")$par)
          jj<-splinefun(seq(.1,3,.1),J3[,n])(hat[r+1,1])
        }
        beta<-(sum((abs(data-hat[11,2]))^(hat[11,1]))/(n*J1[n]))^(1/hat[11,1])
        alpha<-hat[11,1]
        mu<-hat[11,2]
        out<-c(alpha,beta,mu)
      }
    }
    if(method=="ml"){
      f4<-function(x,par){sum(-dweibull(x-par[3],par[1],par[2],log=TRUE))}
      out<-suppressWarnings(optim(c(starts[1],starts[2],starts[3]),fn=f4,x=data,method="Nelder-Mead")$par)
    }
    if(method=="tlm"){
      s.x<-sort(data)
      i2<-seq(2,(n-2))
      i1<-seq(2,(n-1))
      i3<-seq(3,(n-1))
      s.x2<-s.x[2:(n-2)]
      s.x1<-s.x[2:(n-1)]
      s.x3<-s.x[3:(n-1)]
      m1<-mean(data)
      m11<-6/(n*(n-1)*(n-2))*sum((i1-1)*(n-i1)*s.x1)
      m21<-12/(n*(n-1)*(n-2)*(n-3))*(sum((n-i3)*(i3-1)*(i3-2)/2*s.x3)-sum((n-i2)*(n-i2-1)/2*(i2-1)*s.x2))
      f5<-function(p){gamma(1/p+1)*(3*2^(-1/p)-2*3^(-1/p)-1)/(6*gamma(1/p+1)*(2^(-1/p-1)-2*
                                                                                3^(-1/p-1)+4^(-1/p-1))-6*gamma(1/p+1)*(3^(-1/p-1)-4^(-1/p-1)))-(m11-m1)/m21}
      if (f5(0.02)*f5(500)>0){
        out<-list("alpha"=0,"beta"=0,"mu"=0)
      }else{
        alpha<-uniroot(f5,lower=0.02,upper=500)$root
        beta<-m21/(6*gamma(1/alpha+1)*(2^(-1/alpha-1)-2*3^(-1/alpha-1)+4^(-1/alpha-1)-3^(-1/alpha-1)+4^(-1/alpha-1)))
        mu<-m11-beta*gamma(1/alpha+1)*(3*2^(-1/alpha)-2*3^(-1/alpha))
        out<-c(alpha,beta,mu)
      }
    }

    if(method=="moment"){
      n<-length(data)
      f6<-function(u){(gamma(1+3/u)-3*gamma(1+2/u)*gamma(1+1/u)+(gamma(1+1/u))^3)/
          (gamma(1+2/u)-(gamma(1+1/u))^2)^1.5-mean((data-mean(data))^3)/(var(data))^1.5}
      alpha<-uniroot(f6,lower=0.02,upper=500)$root
      beta<-sqrt(var(data)/(2*gamma(1+2/alpha)-(gamma(1+1/alpha))^2))
      mu<-mean(data)-beta*gamma(1+1/alpha)
      out<-c(alpha,beta,mu)
    }

    if(method=="mm1"){
      h1<-function(u){(gamma(1+2/u)-(gamma(1+1/u))^2)/(gamma(1+1/u)-(-log(n/(n+1)))^(1/u))^2-var(data)/(mean(data)-min(data))^2}
      if (h1(0.02)*h1(500)>0){
        out<-list("alpha"=0,"beta"=0,"mu"=0)
      }else{
        alpha<-uniroot(h1,lower=.02,upper=500)$root
        beta<-sd(data)/sqrt(gamma(1+2/alpha)-(gamma(1+1/alpha))^2)
        mu<-mean(data)-beta*gamma(1+1/alpha)
        out<-c(alpha,beta,mu)
      }
    }

    if(method=="mm2"){
      h2<-function(u){(gamma(1+2/u)-(gamma(1+1/u))^2)/(gamma(1+1/u)*(1-n^(-1/u)))^2-var(data)/(mean(data)-min(data))^2}
      if (h2(0.02)*h2(500)>0){
        out<-list("alpha"=0,"beta"=0,"mu"=0)
      }else{
        alpha<-uniroot(h2,lower=.02,upper=500)$root
        beta<-sd(data)/sqrt(gamma(1+2/alpha)-(gamma(1+1/alpha))^2)
        mu<-mean(data)-beta*gamma(1+1/alpha)
        out<-c(alpha,beta,mu)
      }
    }

    if(method=="mm3"){
      h3<-function(u){(gamma(1+2/u)-(gamma(1+1/u))^2)/(gamma(1+1/u)-(log(2))^(1/u))-var(data)/(mean(data)-median(data))^2}
      if (h3(0.02)*h3(500)>0){
        out<-list("alpha"=0,"beta"=0,"mu"=0)
      }else{
        alpha<-uniroot(h3,lower=.02,upper=500)$root
        beta<-sd(data)/sqrt(gamma(1+2/alpha)-(gamma(1+1/alpha))^2)
        mu<-mean(data)-beta*gamma(1+1/alpha)
        out<-c(alpha,beta,mu)
      }
    }

    if(method=="mml1"){
      h4<-function(x,par){-n*log(par[1])+n*par[1]*log(par[2])-(par[1]-1)*log(-log(.5))-(par[1]-1)*sum(log(x-par[3]))+
          sum(((x-par[3])/par[2])^par[1])-log(.5)}
      s.x2<-sort(data)[-1]
      out<-suppressWarnings(optim(c(starts[1],starts[2],starts[3]),fn=h4,x=s.x2,method="Nelder-Mead")$par)
    }

    if(method=="mml2"){
      h5<-function(x,par){-n*log(par[1])+n*par[1]*log(par[2])-(par[1]-1)*sum(log(x-min(x)+par[2]*gamma(1+1/par[1])*n^(-1/par[1])))+
          sum(((x-min(x)+par[2]*gamma(1+1/par[1])*n^(-1/par[1]))/par[2])^par[1])}
      out<-suppressWarnings(optim(c(starts[1],starts[2],starts[3]),fn=h5,x=data,method="Nelder-Mead")$par)
    }

    if(method=="mml3"){
      h6<-function(x,par){-n*log(par[1])+n*par[1]*log(par[2])-(par[1]-1)*sum(log(abs(x-mean(x)+par[2]*gamma(1+1/par[1]))))+
          sum((abs((x-mean(x)+par[2]*gamma(1+1/par[1])))/par[2])^par[1])}
      out<-suppressWarnings(optim(c(starts[1],starts[2]),fn=h6,x=data,method="Nelder-Mead")$par)
      mu<-mean(data)-out[2]*gamma(1+1/out[1])
      out<-c(out[1],out[2],mu)
    }

    if(method=="mml4"){
      h7<-function(x,par){-n*log(par[1])+n*par[1]*log(par[2])-(par[1]-1)*sum(log(abs(x-median(x)+par[2]*(log(2))^(1/par[1]))))+
          sum((x-median(x)+par[2]*(log(2))^(1/par[1]))^(par[1]))}
      out<-suppressWarnings(optim(c(starts[1],starts[2]),fn=h7,x=data,method="Nelder-Mead")$par)
      mu<-median(data)-out[2]*(log(2))^(1/out[1])
      out<-c(out[1],out[2],mu)
    }
  }

  if(location==FALSE){

    if(method=="wml"){
      if (starts[1]>3){
        h8<-function(x,par){sum(-dweibull(x,par[1],scale=par[2],log=TRUE))}
        out<-suppressWarnings(optim(c(starts[1],starts[2]),fn=h8,x=data,method="Nelder-Mead")$par)
      }
      if (n>=100 & starts[1]<3){
        h9<-function(x,par){sum(-dweibull(x,par[1],scale=par[2],log=TRUE))}
        out<-suppressWarnings(optim(c(starts[1],starts[2]),fn=h9,x=data,method="Nelder-Mead")$par)
      }
      if (n<=100 & starts[1]<3){
        J1<-c(
          0.693399510865922	,0.839125940117070,0.891392777598858,0.917864269058571,
          0.934225124007816,0.945033812270250,
          0.952777477374756	,0.958649499627512,0.963284287519060,0.966829917795996,
          0.969837377579613,0.972318590280105,
          0.974430644815776	,0.976343073832322,0.977860635417786,0.979211775254001,
          0.980463131367962,0.981558956270962,
          0.982469493725515	,0.983357555900490,0.984153061789624,0.984872194861913,
          0.985511273135968,0.986114872460018,
          0.986715968519795	,0.987181375441105,0.987701578360708,0.988109191812995,
          0.988528392866783,0.988897231222285,
          0.989311980068170	,0.989617690382056,0.989939910592766,0.990208268465778,
          0.990448612206555,0.990757801777452,
          0.990991853716937	,0.991288148841666,0.991483298379355,0.991670917088605,
          0.991913942276601,0.992115489420633,
          0.992254208167335	,0.992447902562333,0.992607973035660,0.992774850846298,
          0.992933225790896,0.993094360726243,
          0.993191487808104	,0.993349892249525,0.993446763742623,0.993596572362274,
          0.993698424459367,0.993852110559752,
          0.993935526305090	,0.994069584343061,0.994174974190556,0.994288023930874,
          0.994354723058738,0.994469140023294,
          0.994547077302967	,0.994650869484763,0.994718753882411,0.994794515328013,
          0.994890780337426,0.994961617728987,
          0.995034634354075	,0.995098757176877,0.995148620665368,0.995251277668251,
          0.995328721196434,0.995358206531323,
          0.995445169920828	,0.995532541997083,0.995561067480534,0.995616547211185,
          0.995690753099331,0.995749055388893,
          0.995777614171302	,0.995825678321107,0.995865176734199,0.995946387066735,
          0.996006140281170,0.996052316892683,
          0.996119497507878	,0.996122842630375,0.996174812424118,0.996229257476399,
          0.996271273371078,0.996301700864851,
          0.996335198179122	,0.996378436492582,0.996402998079485,0.996453125266745,
          0.996486128783805,0.996550403208885,
          0.996554365199593	,0.996591110435042,0.996655503821045,0.996668235534587)
        J2<-c(
          0.0000000, 0.2746639, 0.5174700, 0.6382124, 0.7100117, 0.7575679,
          0.7918226 ,0.8171166 ,0.8371389 ,0.8530596 ,0.8662250 ,0.8774398,
          0.8862835 ,0.8944473 ,0.9012117 ,0.9075023 ,0.9127426 ,0.9174826,
          0.9219437 ,0.9257018 ,0.9292919 ,0.9324352 ,0.9352715 ,0.9379400,
          0.9403699 ,0.9427507 ,0.9446836 ,0.9467658 ,0.9484713 ,0.9502837,
          0.9518578 ,0.9533995 ,0.9546821 ,0.9561776 ,0.9573394 ,0.9585641,
          0.9594171 ,0.9605684 ,0.9616638 ,0.9626955 ,0.9635044 ,0.9644070,
          0.9651397 ,0.9658877 ,0.9667386 ,0.9673406 ,0.9681123 ,0.9687300,
          0.9694359 ,0.9699877 ,0.9706071 ,0.9710743 ,0.9716842 ,0.9722846,
          0.9727341 ,0.9732096 ,0.9736057 ,0.9741854 ,0.9746275 ,0.9750970,
          0.9753708 ,0.9757499 ,0.9760938 ,0.9765220 ,0.9767770 ,0.9771591,
          0.9775977 ,0.9778742 ,0.9781287 ,0.9784835 ,0.9788126 ,0.9791067,
          0.9794128 ,0.9797964 ,0.9799889 ,0.9802592 ,0.9804821 ,0.9806958,
          0.9810675 ,0.9811123 ,0.9814423 ,0.9815888 ,0.9818914 ,0.9820367,
          0.9823018 ,0.9825038 ,0.9827401 ,0.9829047 ,0.9831506 ,0.9832705,
          0.9834256 ,0.9836325 ,0.9838002 ,0.9840143 ,0.9841296 ,0.9842398,
          0.9844669 ,0.9846324 ,0.9847958 ,0.9850023)
        h10<-function(x,par){(J2[length(x)]/par[1]+mean(log(x))-sum(log(x)*x^par[1])/sum(x^par[1]))^2}
        alpha<-suppressWarnings(optim(starts[[1]],fn=h10,x=data,method="Nelder-Mead")$par)
        beta<-(sum(data^alpha)/(n*J1[n]))^(1/alpha)
        out<-c(alpha,beta)
      }
    }

    if(method=="rank"){
      d<-cor(data,rank(data))*sqrt(var(data))/(sqrt(3)*mean(data))*sqrt((n+1)/(n-1))
      alpha<--log(2)/log(1-d)
      beta<-(sum(data^(alpha))/n)^(1/alpha)
      out<-c(alpha,beta)
    }

    if(method=="reg"){
      k<-seq(1,n)
      fit<-lm(log(-log(1-(k-.3)/(n+.4)))~log(sort(data)))
      alpha<-fit[[1]][[2]]
      beta<-exp(-fit[[1]][[1]]/alpha)
      out<-c(alpha,beta)
    }

    if(method=="moment"){
      h11<-function(u){(gamma(1+2/u)-(gamma(1+1/u))^2)/(gamma(1+1/u))^2-var(data)/(mean(data))^2}
      alpha<-uniroot(h11,lower=0.02,upper=500000)$root
      beta<-mean(data)/gamma(1+1/alpha)
      out<-c(alpha,beta)
    }

    if(method=="wreg"){
      y<-sort(data)
      k<-seq(1,n)
      X<-cbind(rep(1,n),log(-log(1-(k)/(n+1))))
      V<-matrix(0,n,n)
      W<-diag(k/(n+1-k)*1/(log((n+1-k)/(n+1)))*1/(log((n+1-k)/(n+1))))
      U<-solve(t(X)%*%solve(W)%*%X)%*%t(X)%*%solve(W)%*%log(y)
      beta<-exp(U[1])
      alpha<-1/U[2]
      out<-c(alpha,beta)
    }

    if(method=="greg1"){
      y<-sort(data)
      k<-seq(1,n)
      X<-cbind(rep(1,n),log(-log(1-k/(n+1))))
      V<-matrix(0,n,n)
      for(i in 1:n){for(j in i:n){V[i,j]<-i/(n+1-i)*(log((n+1-i)/(n+1))*log((n+1-j)/(n+1)))^(-1)}}
      U<-solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%log(y)
      beta<-exp(U[1])
      alpha<-1/U[2]
      out<-c(alpha,beta)
    }

    if(method=="greg2"){
      y<-sort(data)
      k<-seq(1,n)
      X<-cbind(rep(1,n),log(-log(1-k/(n+1)))+(k*(n-k+1))/((n+1)^2*(n+2))*(log(1-k/(n+1))+1)/((1-k/(n+1))*log(1-k/(n+1)))^2)
      V<-matrix(0,n,n)
      for(i in 1:n){for(j in i:n){V[i,j]<-i/(n+1-i)*(log((n+1-i)/(n+1))*log((n+1-j)/(n+1)))^(-1)}}
      U<-solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%log(y)
      beta<-exp(U[1])
      alpha<-1/U[2]
      out<-c(alpha,beta)
    }

    if(method=="ml"){
      h12<-function(x,par){sum(-dweibull(x,shape=par[1],scale=par[2],log=TRUE))}
      out<-suppressWarnings(optim(c(starts[1],starts[2]),fn=h12,x=data,method="Nelder-Mead")$par)
    }

    if(method=="lm"){
      y<-sort(data)
      k<-seq(1,n)
      m1<-mean(data)
      m2<-2/(n*(n-1))*sum((k-1)*y)-m1
      alpha<--log(2)/log(1-m2/m1)
      beta<-m1/gamma(1+1/alpha)
      out<-c(alpha,beta)
    }

    if(method=="mlm"){
      alpha<-sqrt(pi^2/(6*var(log(data))))
      out<-c(alpha,exp((mean(log(data))+0.5772156)/alpha))
    }

    if(method=="pm"){
      p<-0.31
      qu<-quantile(data,0.6321206)[[1]]
      alpha<-log(-log(1-p))/(log(quantile(data,p)[[1]]/qu))
      out<-c(alpha,qu)
    }

    if(method=="ustat"){
      sa<-sb<-0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          sa<-sa+(log(min(data[i],data[j]))-(log(data[i])+log(data[j]))/2)/log(2)
          sb<-sb-0.57721566*log(min(data[i],data[j]))/log(2)+(log(data[i])+log(data[j]))/
            2*(1+0.57721566/log(2))
        }}
      alpha<--(n*(n-1))/(2*sa)
      beta<-exp(2*sb/(n*(n-1)))
      out<-c(alpha,beta)
    }
  }
  n.p<-ifelse(location==TRUE,3,2)
  pdf0<-function(par,x){shape=par[1]; scale=par[2]; dweibull(x,shape,scale)}
  cdf0<-function(par,x){shape=par[1]; scale=par[2]; pweibull(x,shape,scale)}
  if (location==TRUE){
    pdf0<-function(par,x){shape=par[1]; scale=par[2]; dweibull(x-par[3],shape,scale)}
    cdf0<-function(par,x){shape=par[1]; scale=par[2]; pweibull(x-par[3],shape,scale)}
  }
  log.likelihood=suppressWarnings(sum(log(round(pdf0(out,data), digits=20))))
  von<-u<-c()
  anderson<-c()
  sx<-sort(data)
  for(i in 1:n){
    u[i]<-ifelse(cdf0(out,sx[i])==1,0.99999999,cdf0(out,sx[i]))
    von[i]<-(cdf0(out,sx[i])-(2*i-1)/(2*n))^2
    anderson[i]<-(2*i-1)*log(cdf0(out,sx[i]))+(2*n+1-2*i)*log(1-cdf0(out,sx[i]))
  }
  anderson.stat<-suppressWarnings(-n-mean(anderson))
  von.stat<-suppressWarnings(sum(von)+1/(12*n))
  CAIC<--2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC<--2*log.likelihood + 2*n.p
  BIC<--2*log.likelihood + n.p*log(n)
  HQIC<--2*log.likelihood + 2*log(log(n))*n.p
  ks.stat<-suppressWarnings(ks.test(data, "cdf0",par=out))[[1]][[1]]
  out1<-cbind(AIC, CAIC, BIC, HQIC, anderson.stat, von.stat, ks.stat, log.likelihood)
  oux<-cbind(out[1],out[2])
  colnames(oux)<-c("shape","scale")
  if (location==TRUE){
    oux<-cbind(out[1],out[2],out[3])
    colnames(oux)<-c("alpha","beta","mu")
  }
  colnames(out1)<-c("AIC","CAIC","BIC","HQIC","AD","CVM","KS","log.likelihood")
  list("estimate"=oux,"measures"=out1)
}
dmixture<-function(data,g,K,param){
  if(g!="birnbaum-saunders" &  g!="burr"& g!="chen" &  g!="f" & g!="frechet" & g!="gamma"&
     g!="gompertz" & g!="log-logistic" & g!="log-normal" & g!="lomax" & g!="skew-normal" & g!="weibull"){
    stop ("Baseline distribution is not implemented or misspelled.")}
  # if (sum(param[1:K])!= 1){stop ("The weight vector must be added to one.")}
  y<-rep(NA,length(data));
  if(g=="birnbaum-saunders"){
    den=function(x,par){a=par[1]; b=par[2]; ifelse(x==0,0,(sqrt(b/x)+(b/x)^(3/2))/(2*a*b)*dnorm((sqrt(x/b)-sqrt(b/x))/a))}
  }
  if(g=="burr"){
    den<-function(x,par){a=par[1];b=par[2];a*b*x^(a-1)*(1+x^a)^(-b-1)}
  }
  if(g=="chen"){
    den=function(x,par){a=par[1]; b=par[2]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}
  }
  if(g=="f"){
    den=function(x,par){df1=par[1]; df2=par[2]; df(x,df1,df2)}
  }
  if(g=="frechet"){
    den=function(x,par){a=par[1]; b=par[2]; ifelse(x==0,0,a*exp(-(x/b)^(-a))*(x/b)^(-a-1)/b)}
  }
  if(g=="gamma"){
    den=function(x,par){a=par[1]; b=par[2]; dgamma(x,a,scale=b)}
  }
  if(g=="gompertz"){
    den=function(x,par){a=par[1]; b=par[2]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
  }
  if(g=="log-logistic"){
    den=function(x,par){a=par[1]; b=par[2]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
  }
  if(g=="log-normal"){
    den=function(x,par){a=par[1]; b=par[2]; dlnorm(x,a,sdlog=b)}
  }
  if(g=="lomax"){
    den=function(x,par){a=par[1]; b=par[2]; (a*b)/((1+a*x)^(b+1))}
  }
  if(g=="skew-normal"){
    den=function(x,par){a=par[1]; b=par[2]; l=par[3]; 2/b*dnorm((x-a)/b)*pnorm(l*(x-a)/b)}
  }
  if(g=="weibull"){
    den=function(x,par){a=par[1]; b=par[2]; dweibull(x,a,scale=b)}
  }
  omega<-param[1:K]
  alpha<-param[(K+1):(2*K)]
  beta<-param[(2*K+1):(3*K)]
  lambda<-param[(3*K+1):(4*K)]
  for (i in 1:length(data)){
    pdf<-0
    for (j in 1:K){
      pdf<-pdf+omega[j]*den(data[i],c(alpha[j],beta[j],lambda[j]))
    }
    y[i]<-pdf
  }
  return(y)
}

pmixture<-function(data,g,K,param){
  if(g!="birnbaum-saunders" &  g!="burr"& g!="chen" &  g!="f" & g!="frechet" & g!="gamma"&
     g!="gompertz" & g!="log-logistic" & g!="log-normal" & g!="lomax" & g!="skew-normal" & g!="weibull"){
    stop ("Baseline distribution is not implemented or misspelled.")}
  # if (sum(param[1:K])!= 1){stop ("The weight vector must be added to one.")}
  y<-rep(NA,length(data));
  if(g=="birnbaum-saunders"){
    cum=function(x,par){a=par[1]; b=par[2]; pnorm((sqrt(x/b)-sqrt(b/x))/a)}
  }
  if(g=="burr"){
    cum<-function(x,par){a=par[1];b=par[2];1-(1+x^a)^(-b)}
  }
  if(g=="chen"){
    cum=function(x,par){a=par[1]; b=par[2]; 1-exp(-b*(exp(x^a)-1))}
  }
  if(g=="f"){
    cum=function(x,par){df1=par[1]; df2=par[2]; pf(x,df1,df2)}
  }
  if(g=="frechet"){
    cum=function(x,par){a=par[1]; b=par[2]; exp(-(x/b)^(-a))}
  }
  if(g=="gamma"){
    cum=function(x,par){a=par[1]; b=par[2]; pgamma(x,a,scale=b)}
  }
  if(g=="gompertz"){
    cum=function(x,par){a=par[1]; b=par[2]; 1-exp(-(exp(a*x)-1)*b/a)}
  }
  if(g=="log-logistic"){
    cum=function(x,par){a=par[1]; b=par[2]; 1/((x/b)^(-a)+1)}
  }
  if(g=="log-normal"){
    cum=function(x,par){a=par[1]; b=par[2]; plnorm(x,a,sdlog=b)}
  }
  if(g=="lomax"){
    cum=function(x,par){a=par[1]; b=par[2]; 1-(1+a*x)^(-b)}
  }
  if(g=="skew-normal"){
    cum=function(x,par){a=par[1]; b=par[2]; l=par[3];
    gi<-function(u){
      bb=-1/2*(u-a)^2/b^2;cc=pnorm(l/b*(u-a))
      return(exp(bb)*cc)
    }
    ifelse(x<=-8*b+a,.Machine$double.xmin,2/(b*sqrt(2*pi))*integrate(gi,lower=-8*b+a,upper=x,rel.tol=1e-8)$value)
    }
  }
  if(g=="weibull"){
    cum=function(x,par){a=par[1]; b=par[2]; pweibull(x,a,scale=b)}
  }
  omega<-param[1:K];alpha<-param[(K+1):(2*K)];beta<-param[(2*K+1):(3*K)];lambda<-param[(3*K+1):(4*K)]
  for (i in 1:length(data)){
    cdf<-0
    for (j in 1:K){
      cdf<-cdf+omega[j]*cum(data[i],c(alpha[j],beta[j],lambda[j]))
    }
    y[i]<-cdf
  }
  return(y)
}

rmixture<-function(n, g, K, param)
{
  if(g!="birnbaum-saunders" &  g!="burr"& g!="chen" &  g!="f" & g!="frechet" & g!="gamma"&
     g!="gompertz" & g!="log-logistic" & g!="log-normal" & g!="lomax" &  g!="skew-normal"&
     g!="skew-normal" & g!="weibull")
  {
    stop ("Baseline distribution is not implemented or misspelled.")
  }
  omega <- param[1:K]
  alpha <- param[(K+1):(2*K)]
  beta  <- param[(2*K+1):(3*K)]
  lambda<- param[(3*K+1):(4*K)]
  label <- rep(NA, K)
  y <- rep(NA, n)
  label <- apply( rmultinom(n, 1, omega), 2, which.max )
  if(g=="birnbaum-saunders")
  {
    for (i in 1:n)
    {
      zq <- alpha[label[i]]*qnorm(runif(1))
      y[i] <- beta[label[i]]/4*(zq+sqrt(zq^2+4))^2
    }
  }
  if(g=="burr")
  {
    for (i in 1:n)
    {
      y[i] <-	((1-runif(1))^(-1/beta[label[i]])-1)^(1/alpha[label[i]])
    }
  }
  if(g=="chen")
  {
    for (i in 1:n)
    {
      y[i] <- (log(1-log(1-runif(1))/beta[label[i]]))^(1/alpha[label[i]])
    }
  }
  if(g=="f")
  {
    for (i in 1:n)
    {
      y[i] <-	rf(1,df1=alpha[label[i]],df2=beta[label[i]])
    }
  }
  if(g=="frechet")
  {
    qf <- function(par,u){a=par[1]; b=par[2]; b*(-log(u))^(-1/a)}
    for (i in 1:n)
    {
      y[i] <-	qf(c(alpha[label[i]],beta[label[i]]),runif(1))
    }
  }
  if(g=="gamma")
  {
    for (i in 1:n)
    {
      y[i] <-	rgamma(1,alpha[label[i]],scale=beta[label[i]])
    }
  }
  if(g=="gompertz")
  {
    for (i in 1:n)
    {
      y[i] <-
        log(1-alpha[label[i]]/beta[label[i]]*log(1-runif(1)))/alpha[label[i]]
    }
  }
  if(g=="log-logistic")
  {
    for (i in 1:n)
    {
      y[i] <-	qf(runif(1), alpha[label[i]], beta[label[i]])
    }
  }
  if(g=="log-normal")
  {
    for (i in 1:n)
    {
      y[i] <-	qlnorm(runif(1),meanlog=alpha[label[i]],sdlog=beta[label[i]])
    }
  }
  if(g=="lomax")
  {
    for (i in 1:n)
    {
      y[i] <-	((1-runif(1))^(-1/beta[label[i]])-1)/alpha[label[i]]
    }
  }
  if(g=="skew-normal")
  {
    rskewn<-function(nu,a,b,lamb)
    {
      SN <- c()
      for (i in 1:nu)
      {
        SN[i]<-a+b*(lamb*abs(rnorm(1))+rnorm(1))/sqrt(1+lamb^2)
      }
      return(SN)
    }
    for (i in 1:n)
    {
      y[i] <-	rskewn(1,alpha[label[i]],beta[label[i]],lambda[label[i]])
    }
  }
  if(g=="weibull")
  {
    for (i in 1:n)
    {
      y[i] <-	rweibull(1,alpha[label[i]],scale=beta[label[i]])
    }
  }
  return(y)
}

fitmixture<-function(data,family,K,initial=FALSE,starts){
  n<-length(data)
  N<-50000
  cri<-5e-5
  von<-u<-anderson<-d.single<-d.mix<-cdf0<-pdf0<-clustering<-y<-alpha<-beta<-delta<-gamma<-lambda<-portion<-c()
  alpha.matrix<-beta.matrix<-lambda.matrix<-gamma.matrix<-p.matrix<-matrix(0,ncol=K,nrow=N)
  weight<-alpha<-beta2<-lambda<-delta<-Delta<-Gama<-rep(0,K)
  eu<-eu2<-tau.matrix<-z<-s.pdf<-matrix(0,ncol=K,nrow=n)
  clustering<-rep(0,n)
  clust<-kmeans(data,K,50,1,algorithm="Hartigan-Wong")
  if (family=="birnbaum-saunders"){
    if (initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        inv.data<-1/y[is.finite(1/y)]
        inv.mean<-mean(inv.data)
        inv<-(mean(inv.data))^(-1)
        alpha[i]<-sqrt(2*(sqrt(mean(y)/inv)))
        beta[i]<-sqrt(mean(y)*inv)
      }
      p.matrix[1,]<-portion
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denbir<-function(x,par){(sqrt(par[2]/x)+(par[2]/x)^(3/2))/(2*par[1]*par[2])*dnorm((sqrt(x/par[2])-sqrt(par[2]/x))/par[1])}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denbir(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hbir<-function(x,par){sum(-z[,j]*log(denbir(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hbir,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denbir(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"birnbaum-saunders",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"birnbaum-saunders",K,c(weight,alpha,beta))
  }
  if (family=="burr"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        z1<-c()
        burr.log<-function(x,par){-nn*log(par[1])-nn*log(par[2])-(par[1]-1)*sum(log(x))+(par[2]+1)*sum(log(1+x^par[1]))}
        p.hat<-suppressWarnings(optim(c(1,1),burr.log, x=y, method="BFGS")$par)
        alpha[i]<-p.hat[1]
        beta[i]<-p.hat[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denburr<-function(x,par){a=par[1];b=par[2];a*b*x^(a-1)*(1+x^a)^(-b-1)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denburr(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hburr<-function(x,par){sum(-z[,j]*log(denburr(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hburr,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denburr(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"burr",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"burr",K,c(weight,alpha,beta))
  }
  if (family=="chen"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        y<-(y[2:nn]-min(y))
        sd.y<-(y-mean(y))/sd(y)
        beta[i]<-ifelse(max(sd.y)<=1,-log(1-nn/(nn+.4))/(exp(1)-1),-log(1-length(sd.y[sd.y<=1])/length(sd.y))/(exp(1)-1))
        alpha[i]<-abs(log(log(1-log(1-0.5)/beta[i]))/log(median(y)))
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denchen<-function(x,par){a=par[1]; b=par[2]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denchen(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1;eps<-1;
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hchen<-function(x,par){sum(-z[,j]*log(denchen(data,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hchen,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denchen(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"chen",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"chen",K,c(weight,alpha,beta))
  }
  if (family=="f"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        inv.mean<-mean(1/y)
        mean.y=mean(y)
        alpha[i]<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
        beta[i]<-ifelse(mean.y>1,2*mean.y/(mean.y-1),1)
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denf<-function(x,par){a=par[1]; b=par[2]; df(x,a,b)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denf(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hf<-function(x,par){sum(-z[,j]*log(denf(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hf,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denf(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"f",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"f",K,c(weight,alpha,beta))
  }
  if (family=="frechet"){
    if(initial==FALSE){
      clust<-kmeans(1/data,K)
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(1/data[clust$cluster==i])
        nn<-length(y)
        k<-seq(1,nn)
        X<-cbind(rep(1,nn),log(-log(1-k/(nn+1)))+(k*(nn-k+1))/((nn+1)^2*(nn+2))*(log(1-k/(nn+1))+1)/
                   ((1-k/(nn+1))*log(1-k/(nn+1)))^2)
        V<-matrix(0,nn,nn)
        for(ii in 1:nn){
          for(jj in ii:nn){
            V[ii,jj]<-ii/(nn+1-ii)*(log((nn+1-ii)/(nn+1))*log((nn+1-jj)/(nn+1)))^(-1)
          }
        }
        U<-solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%log(y)
        beta[i]<-exp(U[1])
        alpha[i]<-1/U[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denfre<-function(x,par){a=par[1]; b=par[2]; a/x*(b/x)^a*exp(-(b/x)^a)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denfre(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hfre<-function(x,par){sum(-z[,j]*log(denfre(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hfre,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denfre(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"frechet",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"frechet",K,c(weight,alpha,beta))
  }
  if (family=="gamma"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        alpha[i]<-uniroot(function(u) trigamma(u)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
        beta[i]<-mean(y)/alpha[i]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*dgamma(data,shape=alpha.matrix[1,i],scale=beta.matrix[1,i])
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hgam<-function(x,par){sum(-z[,j]*dgamma(x,shape=par[1],scale=par[2],log=TRUE))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hgam,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*dgamma(data,shape=alpha.matrix[r+1,j],scale=beta.matrix[r+1,j])
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"gamma",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"gamma",K,c(weight,alpha,beta))
  }
  if (family=="gompertz"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        qp3<-y[floor(.75*nn)]
        m.data<-y[floor(.5*nn)]
        a.hat<-log(log(1-.75)/log(1-.5))/(qp3-m.data)
        b.hat<--a.hat*log(0.5)/(exp(a.hat*m.data)-1)
        z1<-c()
        gompertz.log<-function(par){
          z1<--nn*log(par[2])-par[1]*sum(y)+par[2]/par[1]*sum(exp(par[1]*y)-1)
          z1[z1<1e-16]<-1e-16
        }
        p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log)$par)
        beta[i]<-p.hat[1]
        alpha[i]<-p.hat[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    dengom<-function(x,par){a=par[1]; b=par[2]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*dengom(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hgom<-function(x,par){sum(-z[,j]*log(dengom(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hgom,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*dengom(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"gompertz",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"gompertz",K,c(weight,alpha,beta))
  }
  if (family=="log-logistic"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        qp3<-sort(y)[floor(.75*nn)]
        alpha[i]<-log(0.75/(1-0.75))/log(qp3/median(y))
        beta[i]<-median(y)
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denlog<-function(x,par){a=par[1]; b=par[2]; (a*b^(-a)*(x)^(a-1))/(((x/b)^a +1)^2)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denlog(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hlog<-function(x,par){sum(-z[,j]*log(denlog(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hlog,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denlog(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"log-logistic",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"log-logistic",K,c(weight,alpha,beta))
  }
  if (family=="log-normal"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        median.data<-median(y)
        alpha[i]<-log(median.data)
        beta[i]<-sqrt(2*abs(log(mean(y)/median.data)))
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*dlnorm(data,meanlog=alpha.matrix[1,i],sdlog=beta.matrix[1,i])
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hlogn<-function(x,par){sum(-z[,j]*dlnorm(x,meanlog=par[1],sdlog=par[2],log=TRUE))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hlogn,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*dlnorm(data,meanlog=alpha.matrix[r+1,j],sdlog=beta.matrix[r+1,j])
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"log-normal",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"log-normal",K,c(weight,alpha,beta))
  }
  if (family=="lomax"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        a.hat<-1
        b.hat<-(.5^(-1/a.hat)-1)/median(y)
        z1<-c()
        lomax.log<-function(p) {
          z1<--nn*log(p[1])-nn*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
          z1[z1<1e-16]<-1e-16}
        p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log)$par)
        alpha[i]<-p.hat[1]
        beta[i]<-p.hat[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    denlom<-function(x,par){a=par[1]; b=par[2]; (a*b)/((1+a*x)^(b+1))}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*denlom(data,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hlom<-function(x,par){sum(-z[,j]*log(denlom(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hlom,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*denlom(data,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"lomax",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"lomax",K,c(weight,alpha,beta))
  }
  if (family=="weibull"){
    if(initial==FALSE){
      for  (i in 1:K){
        portion[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust$cluster==i])
        nn<-length(y)
        k<-seq(1,nn)
        X<-cbind(rep(1,nn),log(-log(1-k/(nn+1)))+(k*(nn-k+1))/((nn+1)^2*(nn+2))*(log(1-k/(nn+1))+1)/
                   ((1-k/(nn+1))*log(1-k/(nn+1)))^2)
        V<-matrix(0,nn,nn)
        for(ii in 1:nn){
          for(jj in ii:nn){
            V[ii,jj]<-ii/(nn+1-ii)*(log((nn+1-ii)/(nn+1))*log((nn+1-jj)/(nn+1)))^(-1)
          }
        }
        U<-solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%log(y)
        beta[i]<-exp(U[1])
        alpha[i]<-1/U[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-portion
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){z[i,j]<-1}
        }
      }
    }
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) z[,i]<-p.matrix[1,i]*dweibull(data,shape=alpha.matrix[1,i],scale=beta.matrix[1,i])
      d<-cbind(c(1:n),apply(z, 1, which.max))
      z<-matrix(0,nrow=n,ncol=K)
      z[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hweib<-function(x,par){sum(-z[,j]*dweibull(x,shape=par[1],scale=par[2],log=TRUE))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hweib,x=data)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*dweibull(data,shape=alpha.matrix[r+1,j],scale=beta.matrix[r+1,j])
      }
      z<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        maxim<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= maxim){
            maxim<-tau.matrix[i,t]
            count<-t
          }
        }
        z[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmixture(data,"weibull",K,c(weight,alpha,beta))
    cdf0<-pmixture(sort(data),"weibull",K,c(weight,alpha,beta))
  }
  if (family=="skew-normal"){
    S1<-S2<-S3<-z<-tal<- matrix(0, n, K)
    dSN <- function(y, alpha = 0, beta2 = 1, lambda=1){
      dens <- 2*dnorm(y, alpha, sqrt(beta2))*pnorm(lambda*((y - alpha)/sqrt(beta2)))
      return(dens)
    }
    d.mixedSN <- function(x, pi1, alpha, beta2, lambda){
      K <- length(pi1)
      dens <- 0
      for (j in 1:K) dens <- dens + pi1[j]*dSN(x, alpha[j], beta2[j], lambda[j])
      return(dens)
    }
    if(initial==FALSE){
      for(i in 1:K){
        weight[i]<-sum(clust$cluster==i)/n
        y<-sort(data[clust[[1]]==i])
        sk1<-mean((y-mean(y))^3)/(sqrt(var(y)))^3
        sk<-ifelse(abs(sk1)<0.99,sk1,0.99*sign(sk1))
        lambda[i]<-sign(sk)
        beta2[i]<-clust$withinss[i]
        alpha[i]<-mean(y)
      }
    }
    if(initial==TRUE){
      weight<-starts[1:K]
      alpha<-starts[(K+1):(2*K)]
      beta2<-(starts[(2*K+1):(3*K)])^2
      lambda <-starts[(3*K+1):(4*K)]
    }
    for (j in 1:K){
      delta[j] <- lambda[j] / (sqrt(1 + lambda[j]^2))
      Delta[j] <- sqrt(beta2[j])*delta[j]
      Gama[j]  <- beta2[j] - Delta[j]^2
    }
    alpha.old  <- alpha
    Delta.old  <- Delta
    Gama.old   <- Gama
    count      <- 0
    eps        <- 1
    while((eps > cri) && (count <= N)){
      count <- count + 1
      for (j in 1:K){
        Mtij2 <- 1/(1+(Delta[j]^2)*(Gama[j]^(-1)))
        Mtij  <- sqrt(Mtij2)
        alphatij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(data-alpha[j])
        prob  <- pnorm(alphatij/Mtij)
        if(length(which(prob==0))>0) prob[which(prob == 0)]<-.Machine$double.xmin
        E  <- dnorm(alphatij/Mtij)/prob
        u  <- rep(1,n)
        d1 <- dSN(data, alpha[j], beta2[j], lambda[j])
        if(length(which(d1 == 0)) > 0) d1[which(d1==0)]<-.Machine$double.xmin
        d2 <- d.mixedSN(data, weight, alpha, beta2, lambda)
        if(length(which(d2 == 0)) > 0) d2[which(d2==0)]<-.Machine$double.xmin
        tal[,j]    <- d1*weight[j]/d2
        S1[,j]     <- tal[,j]*u
        S2[,j]     <- tal[,j]*(alphatij*u + Mtij*E)
        S3[,j]     <- tal[,j]*(alphatij^2*u + Mtij2 + Mtij*alphatij*E)
        weight[j]  <- (1/n)*sum(tal[,j])
        alpha[j]   <- sum(S1[,j]*data-Delta.old[j]*S2[,j])/sum(S1[,j])
        Delta[j]   <- sum(S2[,j]*(data-alpha[j]))/sum(S3[,j])
        Gama[j]    <- sum(S1[,j]*(data-alpha[j])^2-2*(data-alpha[j])*Delta[j]*S2[,j]+Delta[j]^2*S3[,j])/sum(tal[,j])
        beta2[j]   <- Gama[j] + Delta[j]^2
        lambda[j]  <- ((beta2[j]^(-1/2))*Delta[j])/(sqrt(1-(Delta[j]^2)*(beta2[j]^(-1))))
      }
      weight[K]    <- 1-(sum(weight)-weight[K])
      eps          <- max(sum(c(abs(Delta.old-Delta))),sum(c(abs(alpha.old-alpha))),sum(c(abs(Gama.old-Gama))))
      alpha.old    <- alpha
      Delta.old    <- Delta
      Gama.old     <- Gama
    }
    beta<-sqrt(beta2)
    for (i in 1:n)z[i,max.col(tal)[i]]<-1
    pdf0<-dmixture(data,"skew-normal",K,c(weight,alpha,beta,lambda))
    cdf0<-pmixture(sort(data),"skew-normal",K,c(weight,alpha,beta,lambda))
  }
  for(i in 1:n)clustering[i]<-which(z[i,]==1)[1]
  for(i in 1:n){
    u[i]<-ifelse(cdf0[i]==1,0.99999999,cdf0[i])
    von[i]<-(cdf0[i]-(2*i-1)/(2*n))^2
    anderson[i]<-suppressWarnings((2*i-1)*log(cdf0[i])+(2*n+1-2*i)*log(1-cdf0[i]))
  }
  von.stat<-suppressWarnings(sum(von)+1/(12*n))
  n.p<-K*3-1
  log.likelihood<-suppressWarnings(sum(log(pdf0)))
  I<-seq(1,n)
  ks.stat<-suppressWarnings(max(I/n-cdf0,cdf0-(I-1)/n))
  anderson.stat<-suppressWarnings(-n-mean(anderson))
  CAIC<--2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC<--2*log.likelihood + 2*n.p
  BIC<--2*log.likelihood + n.p*log(n)
  HQIC<--2*log.likelihood + 2*log(log(n))*n.p
  out2<-cbind(AIC, CAIC, BIC, HQIC, anderson.stat, von.stat, ks.stat, log.likelihood)
  colnames(out2)<-c("AIC","CAIC","BIC","HQIC","AD","CVM","KS","log.likelihood")
  if (family=="skew-normal"){
    out1<-cbind(weight,alpha,beta,lambda)
    colnames(out1)<-c("weight","alpha","beta","lambda")
  }else{
    out1<-cbind(weight,alpha,beta)
    colnames(out1)<-c("weight","alpha","beta")
  }
  out3<-clustering
  list("estimate"=out1,"measures"=out2,"cluster"=out3)
}
fitgrouped1<-function(r,f,family,method1,starts,method2){
  if (family=="ge"){
    G<-function(x,par)pexp(x-par[3],rate=par[2])^par[1]
    g<-function(x,par){a=par[1];b=par[2];d=par[3];log(a*dexp(x-d,rate=b)*pexp(x-d,rate=b)^(a-1))}
    EMG<-function(r,f){
      E1i<-E2i<-alpha.hat<-beta.hat<-mu.hat<-con<-vec<-c()
      m<-length(r)
      n<-sum(f)
      mu0<-r[1]-1/n
      x.mid<-(r[-m]+r[-1])/2-mu0
      x.bar<-sum(f*x.mid)/n
      x.bar2<-sum(f*x.mid^2)/n
      f.root<-function(z){6*(z*digamma(z)+z*0.5772156649+1)^2/(6*digamma(z)^2*z^2+12*0.5772156649*
                                                                 digamma(z)*z^2+6*0.5772156649^2*z^2+pi^2*z^2-6*trigamma(z)*z^2+12*z*
                                                                 digamma(z)+12*z*0.5772156649+12)-x.bar^2/x.bar2}
      alpha0<-max(0.001,uniroot(f.root,c(0.001,100))$root)
      beta0<-(alpha0*digamma(alpha0)+alpha0*0.5772156649+1)/(alpha0*x.bar)
      alpha.hat[1]<-alpha0
      beta.hat[1]<-beta0
      mu.hat[1]<-mu0
      j<-2
      tol<-1
      con[1]<-10000
      while (tol>0.5){
        I1<-Vectorize(function(d,par) integrate(function(par,u){
          par[1]*(pexp(u-par[3],rate=par[2]))^(par[1]-1)*dexp(u-par[3],rate=par[2])*
            log(pexp(u-par[3],rate=par[2]))},lower=par[3]+0.00000001,upper=d,par=par)$value,"d")
        I2<-Vectorize(function(d,par) integrate(function(par,u){
          par[1]*(pexp(u-par[3],rate=par[2]))^(par[1]-1)*dexp(u-par[3],rate=par[2])*
            (u-par[3])},lower=par[3]+.0000001,upper=d,par=par)$value,"d")
        hat<-c(alpha.hat[j-1],beta.hat[j-1],mu.hat[j-1])
        E1<-diff(I1(r,hat))/diff(G(r,hat))
        E2<-diff(I2(r,hat))/diff(G(r,hat))
        alpha.hat[j]<--n/sum(f*E1)
        beta.hat[j]<-n/sum(f*E2)
        h<-function(par)-sum(f*log(G(r[-1],c(alpha.hat[j],beta.hat[j],par))-G(r[-m],c(alpha.hat[j],beta.hat[j],par))))
        mu.hat[j]<-suppressWarnings(optimize(h,c(0,mu0))$minimum)
        vec<-c(abs((alpha.hat[j]-alpha.hat[j-1])/(alpha.hat[j-1]+alpha.hat[j])),
               abs((beta.hat[j]-beta.hat[j-1])/(beta.hat[j-1]+beta.hat[j])),
               abs((mu.hat[j]-mu.hat[j-1])/(mu.hat[j-1]+mu.hat[j])))
        con[j]<-max(vec[!is.na(vec)])
        if(abs(con[j]-con[j-1])>0.0001){
          tol<-1
          j<-j+1}
        else{
          tol<-0.25
        }
      }
      return(c((alpha.hat[j-1]+alpha.hat[j])/2,(beta.hat[j-1]+beta.hat[j])/2,(mu.hat[j-1]+mu.hat[j])/2))
    }
  }
  if (family=="weibull"){
    G<-function(x,par)pweibull(x-par[3],shape=par[1],scale=par[2])
    g<-function(x,par){a=par[1];b=par[2];d=par[3];dweibull(x-d,shape=a,scale=b,log=TRUE)}
    EMG<-function(r,f){
      E1i<-E2i<-E3i<-alpha.hat<-beta.hat<-mu.hat<-u<-T1<-T2<-con<-vec<-c()
      m<-length(r)
      R<-r[m]-r[1]
      mu0<-r[1]
      n<-sum(f)
      xmid<-(r[-m]+r[-1])/2-mu0
      x.bar<-sum(f*xmid)/n
      x.bar2<-sum(f*xmid^2)/n
      f.root<-function(z)gamma(1+1/z)^2/gamma(1+2/z)-x.bar^2/x.bar2
      a0<-uniroot(f.root,c(0.05,100))$root
      beta0<-x.bar2/x.bar*gamma(1+1/a0)/gamma(1+2/a0)
      alpha0<-log(-log(1-.995))/log((r[m]-mu0)/beta0)
      alpha.hat[1]<-alpha0
      beta.hat[1]<-beta0
      mu.hat[1]<-mu0
      I1<-Vectorize(function(w) integrate(function(b)exp(-b*w)/b,lower=1,upper=max(2,10/w))$value,"w")
      j<-2
      tol<-1
      con[1]<-10000
      while (tol>0.5){
        u<-((r-mu.hat[j-1])/beta.hat[j-1])^alpha.hat[j-1];d<-diff(G(r,c(alpha.hat[j-1],beta.hat[j-1],mu.hat[j-1])))
        dd<-diff(G(r[1:2],c(alpha.hat[j-1],beta.hat[j-1],mu.hat[j-1])))
        T1<-c((-exp(-u[2])*log(u[2])-I1(u[2]))/dd,diff((-exp(-u[-1])*log(u[-1])-I1(u[-1])))/d[-1])
        T2<-c((-(1+(1+u[2])*log(u[2]))*exp(-u[2])-I1(u[2])+1)/dd,diff((-(1+(1+u[-1])*log(u[-1]))*exp(-u[-1])-I1(u[-1])))/d[-1])
        E1i<-beta.hat[j-1]^alpha.hat[j-1]*diff(-exp(-u)*(1+u))/d
        if(abs(u[1])<1e-29){
          E2i<-T1}
        else{
          E2i<-diff(-exp(-u)*log(u)-I1(u))/d
        }
        if(abs(u[1])<1e-29){
          E3i<-T2}
        else{
          E3i<-diff(-(1+(1+u)*log(u))*exp(-u)-I1(u))/d
        }
        alpha.hat[j]<-alpha.hat[j-1]*n/(-sum(f*E2i)+sum(f*E3i))
        beta.hat[j]<-(sum(f*E1i)/n)^(1/alpha.hat[j-1])
        h<-function(par)-sum(f*log(diff(G(r,c(alpha.hat[j],beta.hat[j],par[1])))))
        mu.hat[j]<-ifelse(alpha.hat[j]<2,mu0,suppressWarnings(optimize(h,c(-1000000,mu0))$minimum))
        vec<-c(abs((alpha.hat[j]-alpha.hat[j-1])/(alpha.hat[j-1]+alpha.hat[j])),
               abs((beta.hat[j]-beta.hat[j-1])/(beta.hat[j-1]+beta.hat[j])),abs((mu.hat[j]-
                                                                                   mu.hat[j-1])/(mu.hat[j-1]+mu.hat[j])))
        con[j]<-max(vec[!is.na(vec)])
        if(abs(con[j]-con[j-1])>0.0001){
          tol<-1
          j<-j+1}
        else{
          tol<-0.25
        }
      }
      return(c((alpha.hat[j-1]+alpha.hat[j])/2,(beta.hat[j-1]+beta.hat[j])/2,(mu.hat[j-1]+mu.hat[j])/2))
    }
  }
  if (family=="birnbaum-saunders"){
    G<-function(x,par)pnorm(1/par[1]*(sqrt((x-par[3])/par[2])-sqrt(par[2]/(x-par[3]))))
    g<-function(x,par){a=par[1];b=par[2];d=par[3];log(1/(2*sqrt(2*pi)*a*(x-d))*(sqrt(b/(x-d))+sqrt((x-d)/b))*exp(-1/(2*a^2)*((x-d)/b+b/(x-d)-2)))}
    EMG<-function(r,f){
      Xij<-E1i<-E2i<-E3i<-alpha.hat<-beta.hat<-mu.hat<-u<-T1<-T2<-con<-vec<-c()
      m<-length(r)
      n<-sum(f)
      x.mid<-(r[-m]+r[-1])/2
      x.bar<-sum(f*x.mid)/n
      x.bar2<-sum(f*x.mid^2)/n
      f.root<-function(z)(5*z^4+4*z^2)/(z^2+2)-sqrt(var(x.mid))/(x.bar)
      mu0<-r[1]-1/n
      alpha0<-max(0,uniroot(f.root,c(0,100))$root)
      beta0<-2*(x.bar-mu0)/(alpha0^2+2);
      alpha0<-sqrt((r[m]-mu0)/beta0)-sqrt(beta0/(r[m]-mu0))/qnorm(.97)
      beta0<-2*(x.bar-mu0)/(alpha0^2+2);
      N<-2000;M<-100;alpha.hat[1]<-alpha0;beta.hat[1]<-beta0;mu.hat[1]<-mu0
      for (j in 1:M){
        for (i in 2:m){
          lo<-G(r[i-1],c(alpha.hat[j],beta.hat[j],mu.hat[j]))
          up<-G(r[i],c(alpha.hat[j],beta.hat[j],mu.hat[j]))
          qz<-alpha.hat[j]*qnorm(runif(N)*(up-lo)+lo)
          Xij<-beta.hat[j]/4*(qz+sqrt(qz^2+4))^2+mu.hat[j]
          E1i[i-1]<-mean(Xij-mu.hat[j])
          E2i[i-1]<-mean(1/(Xij-mu.hat[j]))
        }
        beta.hat[j+1]<-sqrt(sum(f*E1i)/sum(f*E2i))
        alpha.hat[j+1]<-sqrt(sum(f*E1i)/(n*beta.hat[j+1])+beta.hat[j+1]*sum(f*E2i)/n-2)
        h<-function(par)-sum(f*log(G(r[-1],c(alpha.hat[j+1],beta.hat[j+1],par))-
                                     G(r[-m],c(alpha.hat[j+1],beta.hat[j+1],par))))
        mu.hat[j+1]<-optimize(h,c(0,mu0))$minimum
      }
      return(c(alpha.hat[M],beta.hat[M],mu.hat[M]))
    }
  }
  MLG<-function(r,f,starts){
    m<-length(r)-1
    LLog=function(par){L=c();a=par[1];b=par[2];d=par[3];
    for(i in 2:(m+1)){L[i-1]=f[i-1]*log(G(r[i],c(a,b,d))-G(r[i-1],c(a,b,d)))}
    -sum(L)}
    out<-suppressWarnings(optim(starts,fn=LLog,method=method2)$par)
    return(out)
  }
  AMLG<-function(r,f,starts){
    m<-length(r)
    mid<-(r[-1]+r[-m])/2
    ALLog<-function(x,par)-sum(f*g(x,par))
    out<-suppressWarnings(optim(starts,fn=ALLog,x=mid,method=method2)$par)
    return(out)
  }
  Stat<-function(r,f,starts){
    m<-length(r)
    pg<-rep(0,m)
    prob<-anderson<-cramer<-E<-F.hat<-pgks<-Fminus.hat<-Fplus.hat<-rep(0,m-1)
    n<-sum(f)
    F.hat<-c(0,cumsum(f[-(m-1)])/(n+1))
    pg<-G(r,starts)
    for(i in 2:(m-1)){
      anderson[i]<-F.hat[i]^2*log(pg[i+1]/pg[i]*(1-pg[i])/(1-pg[i+1]))+2*F.hat[i]*log((1-pg[i+1])/(1-pg[i]))
      cramer[i]<-F.hat[i]^2*(pg[i+1]-pg[i])-F.hat[i]*(pg[i+1]^2-pg[i]^2)+1/3*(pg[i+1]^3-pg[i]^3)
    }
    frn<-1
    frd<-0
    Fminus.hat<-c((cumsum(f)-frn)/(n+frd))
    Fplus.hat<-c(cumsum(f)/n)
    pgks<-G(r,starts)[-1]
    KS<-max(pgks-Fminus.hat,Fplus.hat-pgks)
    AD<-n*(sum(anderson)-(pg[m]-pg[1])-log((1-pg[m])/(1-pg[1])))
    CVM<-n*(sum(cramer)+1/3*(pg[2]^3-pg[1]^3))
    E<-n*diff(pg)
    Chi.stat<-sum((f-E)^2/E)
    prob<-E/sum(E)
    log.likelihood<-dmultinom(f,n,prob,log=TRUE)
    n.p<-2
    CAIC<--2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
    AIC<--2*log.likelihood + 2*n.p
    BIC<--2*log.likelihood + n.p*log(n)
    HQIC<--2*log.likelihood + 2*log(log(n))*n.p
    return(c(AIC,CAIC,BIC,HQIC,AD,CVM,Chi.stat,KS,log.likelihood))
  }
  if(method1=="aml"){out<-AMLG(r,f,starts)}
  if(method1=="em"){out<-EMG(r,f)}
  if(method1=="ml"){out<-MLG(r,f,starts)}
  out1<-cbind(out[1],out[2],out[3])
  out2<-suppressWarnings(Stat(r,f,out))
  out3<-cbind(out2[1],out2[2],out2[3],out2[4],out2[5],out2[6],out2[7],out2[8],out2[9])
  colnames(out1)<-c("alpha","beta","mu")
  colnames(out3)<-c("AIC","CAIC","BIC","HQIC","AD","Chi-square","CVM","KS","log.likelihood")
  return(list("estimate"=out1,"measures"=out3))
}

fitgrouped2 <- function(r, f, param, start, cdf, pdf, method = "Nelder-Mead", lb = 0, ub = Inf, level = 0.05)
{
  d <- length(start)
  m <- length(r)
  n <- sum(f)
  K <- rep(NA, m)
  T0 <- c(lb, r)
  PDF <- CDF <- function(x){}
  body(CDF) <- bquote(.(cdf))
  body(PDF) <- bquote(.(pdf))
  stat_A <- stat_C <- rep(NA, m-1)
  Alpha <- U <- rep(NA, m)
  Alpha[1] <- f[1]/n
  if( length(param) != d ) stop("The length of parameter vector and initial values must be the same.")
  D2 <- matrix (NA, nrow = d, ncol = d)
  d1_PDF <- d1_PDFc <- matrix(NA, nrow = d, ncol = m)
  d2_PDF<- d2_PDFc <- array(NA, dim = c(d, m, d) )
  integrand0 <- integrand1 <- integrand2 <- integrand3 <- function(x){}
  body(integrand0) <- bquote(.(pdf))
  loglike  <- function(x, par)
  {
    for(i in 1:d) assign(param[i], par[i])
    -sum( f*log( diff(eval(cdf)) )  )
  }
  out <- suppressWarnings( optim(start, fn = loglike, x = T0, method = method)$par )
  mle <- out
  for(k in 1:d) assign(param[k], mle[k])
  if (is.nan( integrand0(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
  I0 <- Vectorize( function(w) integrate(integrand0, lower = lb, upper = w)$value, "w" )
  C0 <- I0(r)
  d0_PDF <- diff( c( 0, C0 ) )
  d0_PDFc <- 1-C0
  E_X <- d0_PDF
  first <- sapply(1:d, function(i) D(pdf, param[i]))
  for (i in 1:d)
  {
    body(integrand1) <- bquote(.(first[[i]]))
    I1 <- Vectorize( function(w) quadinf(integrand1, lb, w)$Q, "w" )
    C1 <- I1(r)
    d1_PDF[i,] <- diff( c( 0, C1 ) )
    d1_PDFc[i,] <- -C1
  }
  for (i in 1:d)
  {
    for(j in i:d)
    {
      second <- D(first[[i]], param[j])
      body(integrand2) <- bquote(.(second))
      I2 <- Vectorize( function(w) quadinf(integrand2, lb, w)$Q, "w" )
      C2 <- I2(r)
      d2_PDF[i,,j] <- diff( c( 0, C2 ) )
      d2_PDF[j,,i] <- d2_PDF[i,,j]
      d2_PDFc[i,,j] <- -C2
      d2_PDFc[j,,i] <- d2_PDFc[i,,j]
      D2[i,j] <- -n*sum( ( E_X*( d2_PDF[i,,j]/d0_PDF - d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2 ) )  )
      D2[j,i] <- D2[i,j]
    }
  }
  for (i in 2:m)
  {
    Prod <- 1
    for (j in 2:i)
    {
      K <- n-sum(f[1:(j-1)])
      if (f[j]==0 && K==0)
      {
        Prod <- Prod
      }else{
        Prod <- Prod*( 1 - f[j]/K )
      }
    }
    Alpha[i] <- 1 - Prod*( 1-f[1]/n )
  }
  U <- CDF(r)
  for (i in 1:(m-1))
  {
    stat_A[i] <- Alpha[i]^2*log( U[i+1]*( 1-U[i] )/( U[i]*( 1-U[i+1] ) ) ) +
      2*Alpha[i]*log( ( 1-U[i+1] )/( 1-U[i] ) )
    stat_C[i] <- Alpha[i]*( U[i+1]-U[i] )*( Alpha[i]-U[i+1]-U[i] )
  }
  Anderson <- n*sum(stat_A) - n*log( 1 - U[m] ) - n*U[m]
  Cramer <- n*sum(stat_C) + n/3*U[m]^3
  KS <- max( abs( Alpha - U ) )
  out1 <- cbind( out, sqrt( diag(solve(D2)) ), out + sqrt(diag(solve(D2)))*qnorm(level/2), out + sqrt(diag(solve(D2)))*qnorm(1 - level/2) )
  colnames(out1) <- c("estimate", "std. error", "lower bound", "upper bound")
  rownames(out1) <- param
  out2 <- rbind( c(Anderson, Cramer, KS) )
  colnames(out2) <- c("AD", "CVM", "KS")
  rownames(out2) <- c("measure")
  return( list( "estimate" = out1,  "measures" = out2 ) )
}


fitbayesJSB<-function(data,n.burn=8000,n.simul=10000){
  lambda.sampler<-function(gamma,delta,xi,x,nrep=1000){
    Max<-max(x)
    n<-length(x)
    target.ratio<-function(lam.new,lam.old,gamma,delta,xi,x){
      n<-length(x)
      n*log(lam.new/lam.old)+(1-gamma*delta)*sum(log((lam.old+xi-x)/(lam.new+xi-x)))-
        delta^2/2*sum((log((x-xi)/(lam.new+xi-x)))^2-(log((x-xi)/(lam.old+xi-x)))^2)
    }
    lambda<-c()
    lambda[1]<-Max-xi+1/n
    for (i in 2:nrep){
      lambda.new<-Max-xi-log(1-runif(1))
      alpha<-target.ratio(lambda.new,lambda[i-1],gamma,delta,xi,x)+(lambda.new-lambda[i-1])
      if(runif(1)<exp(alpha)){
        lambda[i]<-lambda.new
      }
      else{
        lambda[i]<-lambda[i-1]
      }
    }
    return(lambda[nrep])
  }
  xi.sampler<-function(gamma,delta,lambda,x,nrep=1000){
    Min<-min(x)
    Max<-max(x)
    n<-length(x)
    target.ratio<-function(xi.n,xi.o,gam,del,lam,x){sum(log((x-xi.o)/(x-xi.n)))+sum(log((lam+xi.o-x)/(lam+xi.n-x)))-1/2*
        sum((gam+del*log((x-xi.n)/(lam+xi.n-x)))^2)+1/2*sum((gam+del*log((x-xi.o)/(lam+xi.o-x)))^2)
    }
    xi<-c()
    xi[1]<-runif(1,(Max-lambda),Min)
    for (i in 2:nrep){
      xi.new<-runif(1,(Max-lambda),Min)
      alpha<-target.ratio(xi.new,xi[i-1],gamma,delta,lambda,x)
      if(runif(1)<exp(alpha)){
        xi[i]<-xi.new
      }
      else{
        xi[i]<-xi[i-1]
      }
    }
    return(xi[nrep])
  }
  stat.JSB<-function(y,hat){
    n<-length(y)
    cramer<-u<-anderson<-log.like<-ks<-c()
    sx<-sort(y)
    pdf0<-function(par,x){del=par[1];gam=par[2];lam=par[3];xi=par[4];
    del*lam/(sqrt(2*pi)*(x-xi)*(lam+xi-x))*exp(-1/2*(gam+del*log((x-xi)/(lam+xi-x)))^2)
    }
    cdf0<-function(par,z){del=par[1];gam=par[2];lam=par[3];xi=par[4];
    f<-function(x){del*lam/(sqrt(2*pi)*(x-xi)*(lam+xi-x))*
        exp(-1/2*(gam+del*log((x-xi)/(lam+xi-x)))^2)
    }
    integrate(f,lower=xi,upper=z)$value
    }
    for(i in 1:n){
      u[i]<-ifelse(cdf0(hat,sx[i])==1,0.99999999,cdf0(hat,sx[i]))
      cramer[i]<-(cdf0(hat,sx[i])-(2*i-1)/(2*n))^2
      anderson[i]<-(2*i-1)*log(cdf0(hat,sx[i]))+(2*n+1-2*i)*log(1-cdf0(hat,sx[i]))
      log.like[i]<-log(round(pdf0(hat,sx[i]),digits=20))
      ks[i]<-max(i/n-cdf0(hat,sx[i]),cdf0(hat,sx[i])-(i-1)/n)
    }
    anderson.stat<-suppressWarnings(-n-mean(anderson))
    cramer.stat<-suppressWarnings(sum(cramer)+1/(12*n))
    ks.stat<-suppressWarnings(max(ks))
    loglike.stat<-suppressWarnings(sum(log.like))
    out<-c(anderson.stat,cramer.stat,ks.stat,loglike.stat)
    return(out)
  }
  if(n.burn>n.simul){stop ("n.burn must be less than n.simul")}
  gamma.hat<-delta.hat<-lambda.hat<-xi.hat<-c()
  n<-length(data)
  delta.hat[1]<-1
  lambda.hat[1]<-max(data)-min(data)+2/n
  xi.hat[1]<-min(data)-1/n
  gamma.hat[1]<-log(1/median((data-xi.hat[1])/lambda.hat[1])-1);
  for (r in 2:n.simul){
    dd<-log((data-xi.hat[r-1])/(lambda.hat[r-1]+xi.hat[r-1]-data))
    mu.gamma<--delta.hat[r-1]*sum(dd)/n
    sigma.gamma<-sqrt(1/n)
    gamma.hat[r]<-rnorm(1,mu.gamma,sigma.gamma)
    f<-function(x,gamma,lambda,xi,y){
      t1<-gamma*sum(log((y-xi)/(lambda+xi-y)))
      t2<-sum((log((y-xi)/(lambda+xi-y)))^2)
      length(y)*log(x)-t2/2*(x+t1/t2)^2
    }
    fprima<-function(x,gamma,lambda,xi,y){
      t1<-gamma*sum(log((y-xi)/(lambda+xi-y)))
      t2<-sum((log((y-xi)/(lambda+xi-y)))^2)
      length(y)/x-t2*(x+t1/t2)
    }
    delta.hat[r]<-ars(1,f,fprima,x=c(0.001,5,33),m=3,lb=TRUE,xlb=0,gamma=gamma.hat[r],lambda=lambda.hat[r-1],xi=xi.hat[r-1],y=data)
    lambda.hat[r]<-lambda.sampler(gamma.hat[r],delta.hat[r],xi.hat[r-1],data,1000)
    xi.hat[r]<-xi.sampler(gamma.hat[r],delta.hat[r],lambda.hat[r],data,1000)
  }
  deltahat<-mean(delta.hat[n.burn:n.simul])
  gammahat<-mean(gamma.hat[n.burn:n.simul])
  lambdahat<-mean(lambda.hat[n.burn:n.simul])
  xihat<-mean(xi.hat[n.burn:n.simul])
  out1<-cbind(deltahat,gammahat,lambdahat,xihat)
  out3<-stat.JSB(data,c(deltahat,gammahat,lambdahat,xihat))
  out2<-cbind(out3[1],out3[2],out3[3],out3[4])
  colnames(out1)<-c("delta","gamma","lambda","xi")
  colnames(out2)<-c("AD","CVM","KS","log.likelihood")
  list("estimate"=out1,"measures"=out2)
}
fitbayesWeibull<-function(data,n.burn=8000,n.simul=10000){
  mu.sampler<-function(a,b,x,nrep=1000){
    targ.den<-function(mu,a,b,x){(a-1)*sum(log((x-mu)))-sum(((x-mu)/b)^a)}
    Mu<-c()
    Min<-min(x)
    Mu[1]<-Min-1/length(x)
    for (i in 2:nrep){
      Mu.star<-runif(1,(Min-b/2),Min)
      alpha<-exp(targ.den(Mu.star,a,b,x)-targ.den(Mu[i-1],a,b,x))
      if (runif(1)<alpha){
        Mu[i]<-Mu.star
      }
      else{
        Mu[i]<-Mu[i-1]
      }
    }
    return(Mu[nrep])
  }
  stat.W<-function(y,hat){
    pdf0=function(par,x){shape=par[1]; scale=par[2]; dweibull(x-par[3],shape,scale)}
    cdf0=function(par,x){shape=par[1]; scale=par[2]; pweibull(x-par[3],shape,scale)}
    cramer<-u<-anderson<-ks<-log.like<-c()
    n<-length(y)
    sx<-sort(y)
    for(i in 1:n){
      u[i]<-ifelse(cdf0(hat,sx[i])==1,0.99999999,cdf0(hat,sx[i]))
      cramer[i]=(cdf0(hat,sx[i])-(2*i-1)/(2*n))^2
      anderson[i]=(2*i-1)*log(cdf0(hat,sx[i]))+(2*n+1-2*i)*log(1-cdf0(hat,sx[i]))
      log.like[i]<-log(round(pdf0(hat,sx[i]),digits=20))
      ks[i]<-max(i/n-cdf0(hat,sx[i]),cdf0(hat,sx[i])-(i-1)/n)
    }
    anderson.stat=suppressWarnings(-n-mean(anderson))
    cramer.stat=suppressWarnings(sum(cramer)+1/(12*n))
    ks.stat=suppressWarnings(max(ks))
    loglike.stat=suppressWarnings(sum(log.like))
    out<-c(anderson.stat,cramer.stat,ks.stat,loglike.stat)
    return(out)
  }
  if(n.burn>n.simul){stop ("n.burn must be less than n.simul")}
  h<-function(u){(gamma(1+2/u)-(gamma(1+1/u))^2)/(gamma(1+1/u))^2-var(data)/(mean(data))^2}
  alpha0<-uniroot(h,lower=0.02,upper=500000)$root
  beta0<-mean(data)/gamma(1+1/alpha0)
  f.log<-function(x,b,mu,y){
    n<-length(y)
    ss<-function(u){
      s<-0
      for(i in 1:n){s<-s+(y[i]-mu)^u/b^u}
      return(s)
    }
    (n-1)*log(x)+(x-1)*sum(log((y-mu)/b))-ss(x)
  }
  f.logprim<-function(x,b,mu,y){
    n<-length(y)
    s.log<-function(u){
      s<-0
      for(i in 1:n){s<-s+((y[i]-mu)/b)^u*log((y[i]-mu)/b)}
      return(s)
    }
    (n-1)/x+sum(log((y-mu)/b))-s.log(x)
  }
  alpha.hat<-beta.hat<-mu.hat<-c()
  alpha.hat[1]<-alpha0
  beta.hat[1]<-beta0
  mu.hat[1]<-min(data)-1/length(data)
  for (r in 1:n.simul){
    z<-rgamma(1,length(data),1)
    beta.hat[r+1]<-(sum((data-mu.hat[r])^alpha.hat[r])/z)^(1/alpha.hat[r])
    alpha.hat[r+1]<-ars(1,f.log,f.logprim,x=c(.05,.2,2,4,40),ns=200,m=5,emax=64,lb=TRUE,xlb=0,b=beta.hat[r+1],mu=mu.hat[r],y=data)
    mu.hat[r+1]<-mu.sampler(alpha.hat[r+1],beta.hat[r+1],data)
  }
  alphahat<-mean(alpha.hat[n.burn:n.simul])
  betahat<-mean(beta.hat[n.burn:n.simul])
  muhat<-mean(mu.hat[n.burn:n.simul])
  out1<-cbind(alphahat,betahat,muhat)
  out3<-suppressWarnings(stat.W(data,c(alphahat,betahat,muhat)))
  out2<-cbind(out3[1],out3[2],out3[3],out3[4])
  colnames(out1)<-c("alpha","beta","mu")
  colnames(out2)<-c("AD","CVM","KS","log.likelihood")
  list("estimate"=out1,"measures"=out2)
}
fitgsm<-function(data,K){
  N<-20000
  n<-length(data)
  eps<-1
  r<-2
  cri<-10e-6
  n.p<-K
  I<-seq(1,n)
  anderson<-beta.hat<-cdf0<-pdf0<-u<-von<-c()
  omega.hat<-matrix(NA,nrow=N,ncol=K)
  comp.pdf<-Z<-matrix(NA,nrow=n,ncol=K)
  beta0<-sum(data)*.35/.65
  beta.hat[1]<-rgamma(1,shape=K/max(data)*beta0,rate=beta0)
  omega.hat[1,]<-rep(1/K,K)
  while(eps>cri && (r<N))
  {
    for(k in 1:K) comp.pdf[,k]<-dgamma(data,shape=k,rate=beta.hat[r-1])
    pdf<-omega.hat[r-1,]%*%t(comp.pdf)
    for (i in 1:length(data)) Z[i,]<-omega.hat[r-1,]*comp.pdf[i,]/pdf[i]
    omega.hat[r,]<-apply(Z,2,sum)/n
    beta.hat[r]<-sum(Z%*%seq(1,K))/sum(data%*%Z)
    dif<-sum(abs(omega.hat[r,]-omega.hat[r-1,]))/K
    if (dif <= cri)
    {
    eps<-0
    }
    else
    {
    r<-r+1
    }
  }
  cdf0<-pgsm(sort(data),omega.hat[r,],beta.hat[r])
  pdf0<-dgsm(data,omega.hat[r,],beta.hat[r])
  for(i in 1:n)
  {
    u[i]<-ifelse(cdf0[i]==1,0.99999999,cdf0[i])
    von[i]<-(cdf0[i]-(2*i-1)/(2*n))^2
    anderson[i]<-suppressWarnings((2*i-1)*log(cdf0[i])+(2*n+1-2*i)*log(1-cdf0[i]))
  }
  von.stat<-suppressWarnings(sum(von)+1/(12*n))
  log.likelihood<-suppressWarnings(sum(log(pdf0)))
  ks.stat<-suppressWarnings(max(I/n-cdf0,cdf0-(I-1)/n))
  anderson.stat<-suppressWarnings(-n-mean(anderson))
  von.stat<-sum(von)+1/(12*n)
  CAIC<--2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC<--2*log.likelihood + 2*n.p
  BIC<--2*log.likelihood + n.p*log(n)
  HQIC<--2*log.likelihood + 2*log(log(n))*n.p
  out1<-beta.hat[r]
  out2<-as.vector(omega.hat[r,])
  out3<-cbind(AIC, CAIC, BIC, HQIC, anderson.stat, von.stat, ks.stat, log.likelihood)
  colnames(out3)<-c("AIC","CAIC","BIC","HQIC","AD","CVM","KS","log.likelihood")
  list("beta"=out1,"omega"=out2,"measures"=out3)
}

dgsm<-function(data,omega,beta,log = FALSE){
  K<-length(omega)
  comp.pdf<-matrix(NA,nrow=length(data),ncol=K)
  for(k in 1:K) comp.pdf[,k]<-dgamma(data, shape = k, rate = beta)
  pdf<-comp.pdf%*%omega
  if(log==TRUE){pdf<-log(pdf)}
  return(as.vector(pdf))
}

pgsm<-function(data,omega,beta,log.p = FALSE,lower.tail = TRUE){
  K<-length(omega)
  comp.cdf<-matrix(NA,nrow=length(data),ncol=K)
  for(k in 1:K) comp.cdf[,k]<-pgamma(data, shape = k, rate = beta)
  cdf<-comp.cdf%*%omega
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(as.vector(cdf))
}

rgsm<-function(n,omega,beta){
  label <- rep(NA, n)
      y <- rep(NA, n)
  label <- apply( rmultinom(n, 1, omega), 2, which.max )
  for (i in 1:n) y[i] <- rgamma(1, label[i], rate = beta)
  return(y)
}

skewtreg<-function(y, x, Fisher=FALSE){
  if( any(is.na(y)) ) warning('y contains missing values')
  if( any(is.na(x)) ) warning('x contains missing values')
  p<-dim(cbind(y,x))[2]
  OFI<-matrix(NA, p, p)
  n<-length(y)
  DST<- function(x, mu, sigma, lambda, nu=Inf){
    PDF<-c()
    for (i in 1:length(x)){
      d <- ((x[i]-mu)/sigma);A<-lambda*(x[i]-mu)/sigma
      PDF[i] <- 2*exp(lgamma((nu+1)/2)-lgamma(nu/2))/(sqrt(pi*nu)*sigma)*
        (1+d^2/nu)^(-(nu+1)/2)*pt(sqrt((1+nu)/(d^2+nu))*A,1+nu)
    }
    return(PDF)
  }
  streg<-function(Y, X){
    n<-length(Y);
    N<-4000;cri<-10e-5;j<-2;eps<-1
    s2<-s3<-M<-Mu<-K<-A<-d<-PDF<-tau<-Y1<-Y2<-rep(NA,n)
    k<-dim(cbind(Y,X))[2]
    x<-cbind(rep(1,n),X)
    out<-matrix(NA,ncol=(4+k),nrow=N)
    X1<-subset(cbind(Y,X),(Y<quantile(Y,0.8) | Y>quantile(Y,0.2)))
    coeff<-lm(X1[,1]~X1[,2:k],data=data.frame(X1))$coefficients
    Beta<-coeff[1:k]
    out[1,1:k]<-Beta
    mu<-mean(Y)
    m3 <- 1/n*sum((Y-mu)^3)
    sigma<-sqrt(var(Y))[1]
    lambda<-sign(m3/sigma^3)[1]
    nu<-2
    out[1,(k+1):(k+4)]<-c(mu,sigma,lambda,nu)
    Del<-lambda/sqrt(1+lambda^2)*sigma
    Gam<-(1-lambda^2/(1+lambda^2))*sigma^2
    while (eps>0.5 & j<N){
      Y1<-Y-x%*%Beta
      for(i in 1:n){
        d[i]<-((Y1[i]-mu)/sigma)^2
        A[i]<-lambda*(Y1[i]-mu)/sigma
        Mu[i]<-Del/(Gam+Del^2)*(Y1[i]-mu)
        M[i]<-sqrt(Gam/(Gam+Del^2))
        PDF[i]<-DST(Y1[i],mu,sigma,lambda,nu)
        K[i]<-4*nu^(nu/2)*gamma((nu+3)/2)*(d[i]+nu)^(-(nu+3)/2)*pt(sqrt((nu+3)/(d[i]+nu))*A[i],nu+3)/(PDF[i]*gamma(nu/2)*sqrt(pi)*sigma)
        tau[i]<-2*nu^(nu/2)*gamma((nu+2)/2)*(d[i]+nu+A[i]^2)^(-(nu+2)/2)/(PDF[i]*gamma(nu/2)*(pi)*sigma)
        s2[i]<-K[i]*Mu[i]+M[i]*tau[i]
        s3[i]<-K[i]*Mu[i]^2+M[i]^2+Mu[i]*M[i]*tau[i]
      }
      mu<-sum(K*Y1-Del*s2)/sum(K)
      Del<-sum((Y1-mu)*s2)/sum(s3)
      Gam<-(1/n*sum(K*(Y1-mu)^2-2*(Y1-mu)*Del*s2+Del^2*s3))
      sigma<-sqrt(Del^2+Gam)
      lambda<-Del/sqrt(Gam)
      obj<-function(w){sum(log(DST(Y1,mu,sigma,lambda,w)))}
      nu<-suppressWarnings(optimize(obj, c(1,100), tol = 0.000001, maximum = TRUE)$maximum)
      Y2<-Y*K;
      T<-matrix(0,k,k)
      for (i in 1:n){T<-T+((x[i,])%*%t(x[i,])*K[i])}
      Beta<-solve(T)%*%(t(x)%*%Y2-Del*t(x)%*%s2)
      out[j,]<-c(Beta,mu,sigma,lambda,nu)
      if (sum(abs(out[j-1,1:(k+2)]-out[j,1:(k+2)]))<cri || j>=(N-1)){
        eps<-0
      }
      else
      {j<-j+1}
    }
    return(list(Beta=out[j-1,1:k],mu=out[j-1,k+1],sigma=out[j-1,k+2],lambda=out[j-1,k+3],nu=out[j-1,k+4]))
  }
  out<-streg(y,x)
  if (Fisher==TRUE) {
    FI<-function(X,beta,sigma,lambda,nu){
      p<-dim(cbind(y,X))
      n<-p[1]
      T<-matrix(0, nrow=p[2], ncol=p[2])
      nu<-ifelse(nu<=100,nu,100)
      x<-cbind(rep(1,n),X)
      for (i in 1:n){

        f1<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;
        -(nu+1)/(nu*sigma^2*(1+d/nu))*DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f2<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;
        2*(nu+1)*d/(nu^2*sigma^2*(1+d/nu)^2)*DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }


        f3<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        lambda^2*
          (-gamma((nu+1)/2)/(nu*gamma(nu/2)*sqrt(pi*nu))*(nu+1)*(1+(A^2*(nu+1)/(d+nu))/nu)^(-nu/2-3/2)*(A*sqrt((nu+1)/(d+nu))))*
          ((nu+1)/(d+nu))/(sigma^2*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f4<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        lambda*
          (-gamma((nu+1)/2)/(nu*gamma(nu/2)*sqrt(pi*nu))*(nu+1)*(1+(A^2*(nu+1)/(d+nu))/nu)^(-nu/2-3/2)*(A*sqrt((nu+1)/(d+nu))))*
          A*sqrt(nu+1)*(d+nu)^(-3/2)*(w-(x[i,]%*%beta)[[1]])*sqrt((nu+1)/(d+nu))/(sigma^3*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f5<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        -lambda*dt(sqrt((nu+1)/(d+nu))*A,df=nu+1)*sqrt(nu+1)*(d+nu)^(-3/2)*(w-(x[i,]%*%beta)[[1]])/
          (sigma^3*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f6<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        -lambda*
          (-gamma((nu+1)/2)/(nu*gamma(nu/2)*sqrt(pi*nu))*(nu+1)*(1+(A^2*(nu+1)/(d+nu))/nu)^(-nu/2-3/2)*(A*sqrt((nu+1)/(d+nu))))*
          sqrt((nu+1)/(d+nu))*A*sqrt(nu+1)*(d+nu)^(-3/2)*(w-(x[i,]%*%beta)[[1]])/
          (sigma^3*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f7<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        (-gamma((nu+1)/2)/(nu*gamma(nu/2)*sqrt(pi*nu))*(nu+1)*(1+(A^2*(nu+1)/(d+nu))/nu)^(-nu/2-3/2)*(A*sqrt((nu+1)/(d+nu))))*
          A^2*(nu+1)*(d+nu)^(-3)*(w-(x[i,]%*%beta)[[1]])^2/(sigma^4*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f8<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        -sqrt(nu+1)*dt(sqrt((nu+1)/(d+nu))*A,df=nu+1)*lambda*(d+nu)^(-3/2)*(w-(x[i,]%*%beta)[[1]])/(sigma^3*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f9<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        sqrt(nu+1)*dt(sqrt((nu+1)/(d+nu))*A,df=nu+1)*3*A*(d+nu)^(-5/2)*(w-(x[i,]%*%beta)[[1]])^2/(sigma^4*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f10<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        -sqrt(nu+1)*dt(sqrt((nu+1)/(d+nu))*A,df=nu+1)*A*(d+nu)^(-3/2)/(sigma^2*pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f11<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        -lambda^2*((nu+1)/(d+nu))*(dt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2/
          (sigma^2*(pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2)*DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f12<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        lambda*sqrt((nu+1)/(d+nu))*(dt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2*
          A*sqrt(nu+1)*(d+nu)^(-3/2)*(w-(x[i,]%*%beta)[[1]])/(sigma^3*(pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2)*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f13<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        A*sqrt(nu+1)*(d+nu)^(-3/2)*(w-(x[i,]%*%beta)[[1]])*lambda*sqrt((nu+1)/(d+nu))*
          (dt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2/(sigma^3*(pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2)*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        f14<-function(w){d<-(w-(x[i,]%*%beta)[[1]])^2/sigma^2;                 A<-lambda*(w-(x[i,]%*%beta)[[1]])/sigma;
        -A^2*(nu+1)*(d+nu)^(-3)*(w-(x[i,]%*%beta)[[1]])^2*
          (dt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2/(sigma^4*(pt(sqrt((nu+1)/(d+nu))*A,df=nu+1))^2)*
          DST((w-(x[i,]%*%beta)[[1]]), 0, sigma, lambda, nu)
        }

        T1<- quadinf(f1, -Inf, Inf)$Q
        T2<- quadinf(f2, -Inf, Inf)$Q
        T3<- quadinf(f3, -Inf, Inf)$Q
        T4<- quadinf(f4, -Inf, Inf)$Q
        T5<- quadinf(f5, -Inf, Inf)$Q
        T6<- quadinf(f6, -Inf, Inf)$Q
        T7<- quadinf(f7, -Inf, Inf)$Q
        T8<- quadinf(f8, -Inf, Inf)$Q
        T9<- quadinf(f9, -Inf, Inf)$Q
        T10<- quadinf(f10, -Inf, Inf)$Q
        T11<- quadinf(f11, -Inf, Inf)$Q
        T12<- quadinf(f12, -Inf, Inf)$Q
        T13<- quadinf(f13, -Inf, Inf)$Q
        T14<- quadinf(f14, -Inf, Inf)$Q
        T<-T+x[i,]%*%t(x[i,])*(T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+T12+T13+T14)
      }
      return(-T)
    }
    OFI<-FI(x,out$Beta,out$sigma,out$lambda,out$nu)
    Std.Error<-sqrt(diag(solve(OFI)))*out$sigma
    t.value<-out$Beta/Std.Error
    tail<-cbind((1-pt(t.value,n-p)),pt(t.value,n-p))
    p.value<-2*apply(tail,1,min)
  }
  Error<-(y-cbind(1,x)%*%out$Beta)
  S.E<-sum(Error^2)
  S.T<-sum((y-mean(y))^2)

  out1<-cbind(out$Beta)
  colnames(out1)<-c("Estimate")
  rownames(out1)<-c("beta.0",rownames(out1[2:p,], do.NULL = FALSE, prefix = "beta."))

  if (Fisher==TRUE) {
    out1<-cbind(out$Beta,Std.Error,t.value,p.value)
    colnames(out1)<-c("Estimate", "Std. Error", "T statistic", "P value")
    rownames(out1)<-c("beta.0",rownames(out1[2:p,], do.NULL = FALSE, prefix = "beta."))
  }
  out2<-cbind(min(Error),quantile(Error,0.25)[[1]],quantile(Error,0.50)[[1]],mean(Error),quantile(Error,0.75)[[1]],max(Error))
  colnames(out2)<-c("Min", "1Q", "Median", "Mean", "3Q", "Max")
  F.value<-(S.T-S.E)/(p-1)*(n-p)*S.E

  out3<-cbind(F.value,p-1,n-p, 1-pf(F.value,p-1,n-p))
  colnames(out3)<-cbind("Value", "DF1", "DF2", "P value")
  rownames(out3)<-c("F-statistic")

  out4<-cbind(out$sigma, p-1)
  colnames(out4)<-cbind("Value", "DF")
  rownames(out4)<-c("Residual Std. Error")

  out5<-cbind(1-S.E/S.T, 1-(n-1)/(n-p)*(S.E/S.T))
  colnames(out5)<-cbind("Non-adjusted", "Adjusted")
  rownames(out5)<-c("Multiple R-Squared")

  out6<-cbind(out$mu,out$sigma,out$lambda,out$nu)
  colnames(out6)<-c("Location", "Scale", "Skewness", "DF")

  if (Fisher==TRUE) {
    colnames(OFI)<-NULL
    out7<-OFI
    colnames(out7)<-c("beta.0",colnames(out7[,2:p], do.NULL = FALSE, prefix = "beta."))
    rownames(out7)<-c("beta.0",rownames(out7[2:p,], do.NULL = FALSE, prefix = "beta."))
  }
  list("Coefficients:"=out1,
       "Residuals:"=out2,
       "F:"=out3,
       "MSE:"=out4,
       "R2:"=out5,
       "Estimated Parameters for Error Distribution:"=out6,
       "Observed Fisher Information Matrix:"=ifelse(Fisher==TRUE,out7,"Not requested"))
}



djsb<-function(data, param, log = FALSE){
  n<-length(data)
  pdf<-rep(NA,n)
  if(param[1]<0 || param[3]<0)
  {
    return(message ("Error: either delta or lambda is negative"))
  }
  if(any(data<param[4]) || any(data>param[3]+param[4]))
  {
    return (message ("Error: some element(s) of input data do not belong to support defined for the Johnson's SB distribution"))
  }
  else
  {
    for(i in 1:n)
    {
      pdf[i]<-suppressWarnings(param[1]*param[3]/(sqrt(2*pi)*(data[i]-param[4])*(param[3]+param[4]-data[i]))*
                                 exp(-1/2*(param[2]+param[1]*log((data[i]-param[4])/(param[3]+param[4]-data[i])))^2))
    }
    suppressWarnings(if(log==TRUE) pdf<-log(pdf))
    return(pdf)
  }
}

pjsb<-function(data, param, log.p = FALSE, lower.tail = TRUE){
  n<-length(data)
  cdf<-rep(NA,n)
  if(param[1]<0 || param[3]<0)
  {
    return(message ("Error: either delta or lambda is negative"))
  }
  if(any(data<param[4]) || any(data>param[3]+param[4]))
  {
    return (message ("Error: some element(s) of input data do not belong to support defined for the Johnson's SB distribution"))
  }
  else
  {
    f<-function(x) param[1]*param[3]/(sqrt(2*pi)*(x-param[4])*(param[3]+param[4]-x))*exp(-1/2*(param[2]+param[1]*log((x-param[4])/(param[3]+param[4]-x)))^2)
    for(i in 1:n)
    {
      cdf[i]<-suppressWarnings(quadinf(f, param[4], data[i])$Q)
    }
    if(log.p==TRUE  & lower.tail == FALSE) cdf<-suppressWarnings(log(1-cdf))
    if(log.p==TRUE  & lower.tail == TRUE)  cdf<-suppressWarnings(log(cdf))
    if(log.p==FALSE & lower.tail == FALSE) cdf<-suppressWarnings(1-cdf)
    return(cdf)
  }
}

rjsb<-function(n, param){
  data<-rep(NA,n)
  if(param[1]<0 || param[3]<0)
  {
    return(message ("Either delta or lambda is negative"))
  }
  else
  {
    z <- (rnorm(n)-param[2])/param[1]
    data <- (param[3]*exp(z) + param[4]*(exp(z)+1))/(exp(z)+1)
    return(data)
  }
}
fitJSB <- function(y, n.burn = 8000, n.simul = 10000)
{
  stat.JSB <- function(data, hat)
  {
    n  <- length(data)
    sx <- sort(data)
    cramer <- u <- anderson <- log.like <- ks <- c()
    pdf0 <- function(par, x)
    {
      del = par[1]; gam = par[2]; lam = par[3]; xi = par[4]
      del*lam/( sqrt(2*pi)*(x - xi)*(lam + xi - x) )*exp(-1/2*( gam + del*
                                                                  log( (x - xi)/(lam + xi - x) ) )^2)
    }
    cdf0 <- function(par, z)
    {
      del = par[1]; gam = par[2]; lam = par[3]; xi = par[4];
      f <- function(x) del*lam/( sqrt(2*pi)*(x - xi)*(lam + xi - x) )*
        exp(-1/2*(gam + del*log( (x - xi)/(lam + xi - x) ) )^2)
      up <- min(lam + xi , z)
      #quadinf(f, xi, up)$Q
      integrate(f, lower = xi, upper = up)$value
    }
    for(i in 1:n)
    {
      u[i] <- ifelse( cdf0(hat, sx[i]) >= 1, 0.9999999, cdf0(hat, sx[i]) )
      cramer[i]   <- ( u[i] - (2*i - 1)/(2*n) )^2
      anderson[i] <- (2*i - 1)*log( u[i] ) + (2*n + 1 - 2*i)*log( 1 - u[i] )
      log.like[i] <- log( round(pdf0(hat, sx[i]), digits = 20) )
      ks[i]       <- max( i/n - u[i], u[i] - (i - 1)/n )
    }
    anderson.stat <- suppressWarnings( -n - mean(anderson) )
    cramer.stat   <- suppressWarnings( sum(cramer) + 1/(12*n) )
    loglike.stat  <- suppressWarnings( sum(log.like) )
    ks.stat  <- suppressWarnings( max(ks) )
    out <- c(anderson.stat, cramer.stat, ks.stat, loglike.stat)
    return(out)
  }
  fitbayesJSB.restricted <- function(data, n.burn = n.burn, n.simul = n.simul)
  {
    xi <- min(data) - 1.34
    lambda <- max(data) - xi + 3.8
    if( n.burn>n.simul ) stop ("n.burn must be less than n.simul")
    r.y <- range(data)
    if(r.y[1] == xi) stop ("minimum of data must be greater than parameter xi")
    if(lambda < - r.y[2] + r.y[1] - 3) stop ("parameter lambda must be admissible")
    gamma.hat <- delta.hat <- c()
    n <- length(data)
    delta.hat[1] <- 1
    gamma.hat[1] <- log( 1/median((data - xi)/lambda) - 1 )
    for (r in 2:n.simul)
    {
      dd <- log( (data - xi)/(lambda + xi - data) )
      mu.gamma <- -delta.hat[r-1]*sum(dd)/n
      sigma.gamma <- sqrt(1/n)
      gamma.hat[r] <- rnorm(1, mu.gamma, sigma.gamma)
      f <- function(x, gamma, lambda, xi, y)
      {
        t1 <- gamma*sum( log( (y - xi)/(lambda + xi - y) ) )
        t2 <- sum( (log( (y - xi)/(lambda + xi - y) ) )^2)
        length(y)*log(x) - t2/2*(x + t1/t2)^2
      }
      fprima <- function(x, gamma, lambda, xi, y)
      {
        t1 <- gamma*sum( log((y - xi)/(lambda + xi - y)) )
        t2 <- sum( (log( (y - xi)/(lambda + xi - y) ) )^2 )
        length(y)/x - t2*(x + t1/t2)
      }
      delta.hat[r] <- ars(1, f, fprima, x = c(.001, 5, 300), m = 3,
                          lb = TRUE, xlb = 0, gamma = gamma.hat[r],
                          lambda = lambda, xi = xi, y = data)
    }
    deltahat  <- mean( delta.hat[n.burn:n.simul] )
    gammahat  <- mean( gamma.hat[n.burn:n.simul] )
    lambdahat <- lambda
    xihat     <- xi
    list("estimate" = c(deltahat, gammahat, lambdahat, xihat))
  }

  xi.mom <- min(y) - 1.34
  lambda.mom <- max(y) - xi.mom + 3.8
  mu <- (mean(y) - xi.mom)/lambda.mom
  s.d <- sd(y)/lambda.mom
  delta.mom <- mu*(1 - mu)/s.d + s.d/4*( 1/(mu*(1 - mu)) - 8 )
  gamma.mom <- delta.mom*log( (1 - mu)/mu ) + (0.5 - mu)/delta.mom
  estim.mom <- c(delta.mom, gamma.mom, lambda.mom, xi.mom)

  xi.cml <- min(y) - 1.34
  lambda.cml <- max(y)
  fi <- log( (y - xi.cml)/(lambda.cml + xi.cml - y) )
  delta.cml <- 1/sd(fi)
  gamma.cml <- -mean(fi)*delta.cml
  estim.cml <- c(delta.cml, gamma.cml, lambda.cml, xi.cml)

  xi.kb <- min(y) - 1.3
  lambda.kb <- max(y) - xi.kb + 3.8
  f.95 <- log(  ( quantile(y, 0.95)[[1]] - xi.kb )/
                  ( lambda.kb + xi.kb - quantile(y, 0.95)[[1]] )  )
  f.50 <- log(  ( quantile(y, 0.50)[[1]] - xi.kb )/
                  ( lambda.kb + xi.kb - quantile(y, 0.50)[[1]] )  )
  delta.kb <- qnorm(.95)/(f.95 - f.50)
  gamma.kb <- -delta.kb*f.50
  estim.kb <- c(delta.kb, gamma.kb, lambda.kb, xi.kb)

  out <- fitbayesJSB.restricted(y, n.burn = 8000, n.simul = 10000)$estimate
  estim.bayes  <- out

  stat.mom <- stat.JSB(y, estim.mom)
  stat.cml <- stat.JSB(y, estim.cml)
  stat.kb  <- stat.JSB(y, estim.kb)
  stat.bayes <- stat.JSB(y, estim.bayes )
  out.stat   <- matrix(c(stat.mom, stat.cml, stat.kb, stat.bayes),
                       nrow = 4, ncol = 4, byrow = TRUE)
  out.estim  <- matrix(c(estim.mom, estim.cml, estim.kb, estim.bayes),
                       nrow = 4, ncol = 4, byrow = TRUE)
  colnames(out.estim) <- c("delta", "gamma", "lambda", "xi")
  rownames(out.stat)  <- c("MM", "CML", "KB", "Bayesian")
  rownames(out.estim) <- c("MM", "CML", "KB", "Bayesian")
  colnames(out.stat)  <- c("AD", "CVM", "KS", "log.likelihood")
  out <- list( out.estim, out.stat )
  names(out) <- c( paste0("estimate"), paste0("statistic") )
  return( out )
}
