`el.test.bfgs` <-
function(x,mu,lam,maxit = 100,tol=1e-07){
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    mu <- as.vector(mu)
    if (length(mu) != p) 
        stop("mu must have same dimension as observation vectors.")
    if (n <= p) 
        stop("Need more observations than length(mu) in el.test().")
    z <- t(t(x) - mu)
    TINY <- sqrt(.Machine$double.xmin)
    scale <- mean(abs(z)) + TINY
    z <- z/scale
    if (!missing(lam)) {
        lam <- as.vector(lam)
        lam <- lam * scale
        if (logelr(z, lam,n) > 0) 
            lam <- rep(0, p)
    }
    if (missing(lam)) 
        lam <- rep(0, p)
    if (tol < TINY) 
        tol <- TINY
    func<-function(lam){
        -sum(logstar(1+z%*%t(t(lam)),n))
    }
    dfunc<-function(lam){ 
        p=length(lam)
        temp1<-1+z%*%lam
        temp2<-sqrt(-ddlogstar(temp1,n))
        J<-matrix(rep(temp2,p),ncol=p)*z
        y<-dlogstar(temp1,n)/temp2  
        -t(J)%*%y
    } 
    bfgsresult<-bfgsmin(lam,tol,func,dfunc,maxit)
    lamnew<-bfgsresult$x 
    w<-as.vector(1+z%*%lamnew)
    llr<-2*sum(ans<-logstar(as.vector(1+z%*%lamnew),n))
    mu0=apply(z*scale,2,function(dum) {sum(dum*(w/n))})
    list("-2LLR"=llr,Pval = 1 - pchisq(llr, df = p),lambda=as.vector(lamnew)/scale,wts=w/n,nits=bfgsresult$iter,mu0=mu0)
}

