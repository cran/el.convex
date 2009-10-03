`el.test.damped` <-
function(x,mu,lam,maxit = 200,tol=1e-07){
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
    lamold<-lam
    lold<- logelr(z,lamold,n)
    iit=0
    diff.lam=tol+1
    diff.l=tol+1        
    while ((iit<=maxit)&((diff.lam>tol)|(diff.l>tol))){
        arg <- 1 + z %*% lamold
        wts1 <- as.vector(dlogstar(arg,n))
        wts2 <- as.vector(-ddlogstar(arg,n))^0.5
        grad <- as.matrix(-z * wts1)
        grad <- as.vector(rowsum(grad, rep(1, nrow(grad))))
        hess <- z * wts2
        svdh <- svd(hess)
        if (min(svdh$d) < max(svdh$d) * tol) 
            svdh$d <- svdh$d + max(svdh$d) * tol
        nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
        nstep <- as.vector(nstep %*% matrix(wts1/wts2, n, 1))
        a<-steplength(lamold,nstep,func,dfunc,.0001,.9)
        lamnew<- lamold+a*nstep
        lnew<- logelr(z,lamnew,n)
        diff.lam<-max(abs(lamnew-lamold))
        diff.l<-abs(lnew-lold)
        lold<-lnew
        lamold<-lamnew
        iit<-iit+1
    }
    list("-2LLR"=-2*lnew,Pval = 1 - pchisq(-2 * lnew, df = p),lambda=lamnew/scale,wts=wts1/n,nits=iit,mu=apply(z*scale,2,function(x) {sum(x*(wts1/n))}))
}

