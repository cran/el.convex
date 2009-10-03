logstar<-function(z,n){ # pseudo-logarithm function, corresponding to llog in emplik
    ans<-rep(NA, length(z))
    sel<-z>(1/n)
    ans[sel]<-log(z[sel])
    ans[!sel]<-log(1/n)-1.5+2*n*z[!sel]-n^2*z[!sel]^2/2
    ans
}

dlogstar<-function(z,n){ # first derivative of pseudo-logarithm function, corresponding to llogp in emplik
    ans<-rep(NA, length(z))
    sel<-z>(1/n)
    ans[sel]<-1/(z[sel])
    ans[!sel]<-2*n-n^2*z[!sel]
    ans
}  
  
ddlogstar<-function(z,n){ # second derivative of pseudo-logarithm function, corresponding to llogpp in emplik
    ans<-rep(NA, length(z))
    sel<-z>(1/n)
    ans[sel]<- -1/((z[sel])^2)
    ans[!sel]<- -n^2
    ans
}

logelr<-function(z,lambda,n){ #log of empirical likelihood ratio
    -sum(logstar(1+z%*%t(t(lambda)),n))
}



dfpmin<-function(x, gtol, func, dfunc,ITMAX){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  n<-length(x)
  xnew<-rep(NA,n)
  gold<-rep(NA,n)
  dg<-rep(NA,n)
  hdg<-rep(NA,n)

  fx=func(x)
  g=dfunc(x)
  p=-g
  hessin=diag(n)
  sum1=sum(x^2)
  stpmax=STPMX*max(sqrt(sum1),n);
  iter=1

  while(iter<=ITMAX){
    ln<-lnsrch(x, fx, g, p, func, stpmax)
    fx<-ln$f
    xnew<-ln$x
    p<-xnew-x
    x<-xnew
    temp<-abs(x)
    temp[temp<1]<-1
    temp2<-abs(p)/temp
    test<-max(temp2)
    if(test<TOLX) break
    gold<-g
    g<-dfunc(x)
    temp2<-abs(g)*temp
    test<-max(temp2)/max(fx,1)
    if(test<gtol) break
    dg<-g-gold
    hdg<-hessin%*%dg
    fac<-sum(dg*p)
    fae<-sum(dg*hdg)
    sumdg<-sum(dg^2)
    sump<-sum(p^2)
    if(fac>sqrt(EPS*sumdg*sump)){
      hessin<-hessin+p%*%t(p)/fac-hdg%*%t(hdg)/fae
    }
    p=-hessin%*%g
    iter=iter+1
  }
  list(x=x,f=fx,iter=iter)
}


lnsrch<-function(xold, fold, g, p, func, stpmax){ # g is gradient, p is direction
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  n=length(xold)
  sum1<-sqrt(sum(p^2))
  if(sum1 > stpmax) p<-p*stpmax/sum1 # Scale if attempted step is too big.
  slope<-sum(p*g)
  if(slope>=0) stop("slope error.")
  temp<-abs(xold) #compute lambda_min
  temp[temp<1]<-1
  test<-max(abs(p)/temp)
  alamin<-TOLX/test
  alam<-1 #Always try full Newton step first.
  while(TRUE){
    x<-xold+alam*p
    f<-func(x)
    if(alam<alamin){
      x<-xold
      break
    } else {
        if(f<=fold+ALF*alam*slope){
          break
        } else {
            if(alam==1){
              tmplam<- -slope/2/(f-fold-slope)
            } else{
                rhs1<-f-fold-alam*slope
                rhs2<-f2-fold-alam2*slope
                a<-(rhs1/(alam^2)-rhs2/(alam2^2))/(alam-alam2)
                b<-(-alam2*rhs1/(alam^2)+alam*rhs2/(alam2^2))/(alam-alam2)
                if(a==0){
                  tmplam<- -slope/2/b
                } else{
                    disc<-b^2-3*a*slope
                    if(disc<0){
                      tmplam<-.5*alam
                    } else{
                        if( b<=0){
                          tmplam<-(-b+sqrt(disc))/3/a
                        } else tmplam<--slope/(b+sqrt(disc))
                      }
                  }
                if(tmplam>.5*alam)
                  tmplam<-.5*alam
              }
          }
      }
    alam2<-alam
    f2<-f
    alam<-max(tmplam,.1*alam)                   
  }
  list(x=x,f=f)
}


bfgsmin<-function(x, gtol, func, dfunc,ITMAX){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  n<-length(x)
  xnew<-rep(NA,n)
  gold<-rep(NA,n)
  dg<-rep(NA,n)
  hdg<-rep(NA,n)

  fx=func(x)
  g=dfunc(x)
  p=-g
  hessin=diag(n)
  sum1=sum(x^2)
  stpmax=STPMX*max(sqrt(sum1),n);
  iter=1

  while(iter<=ITMAX){
    ln<-lnsrch(x, fx, g, p, func, stpmax)
    fx<-ln$f
    xnew<-ln$x
    p<-xnew-x
    x<-xnew
    temp<-abs(x)
    temp[temp<1]<-1
    temp2<-abs(p)/temp
    test<-max(temp2)
    if(test<TOLX) break
    gold<-g
    g<-dfunc(x)
    temp2<-abs(g)*temp
    test<-max(temp2)/max(fx,1)
    if(test<gtol) break
    dg<-g-gold
    hdg<-hessin%*%dg
    fac<-sum(dg*p)
    fae<-sum(dg*hdg)
    sumdg<-sum(dg^2)
    sump<-sum(p^2)
    if(fac>sqrt(EPS*sumdg*sump)){
      fac=1/fac
      fad=1/fae
      dg=fac*p-fad*hdg   
      hessin<-hessin+p%*%t(p)*fac-hdg%*%t(hdg)*fad+dg%*%t(dg)*fae
    }
    p=-hessin%*%g
    iter=iter+1
  }
  list(x=x,f=fx,iter=iter)
}





mnbrak<-function(f,a,b,glimit=100){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  fa<-f(a)
  fb<-f(b)
  if(fb>fa){
    tempt<-b
    b<-a
    a<-tempt
  }    
  c<-b+gold*(b-a)
  fc<-f(c)
  if(fb<fc)  return(c(a,b,c,fa,fb,fc))
  while(fb>fc){
    r<-(b-a)*(fb-fc)
    q<-(b-c)*(fb-fa)
    if((q-r)>0){
      denom<-max(q-r,EPS)
    } else{
      denom<--max(r-q,EPS)
    }
    u<-b-((b-c)*q-(b-a)*r)/(2*denom)
    ulim<-b+glimit*(c-b)
    if(((b-u)*(u-c))>0){
      fu<-f(u)
      if(fu<fc){
        a<-b
        b<-u
        fa<-fb
        fb<-fu
        return(c(a,b,c,fa,fb,fc))
      } else if(fu>fb){
        c<-u
        fc<-fu
        return(c(a,b,c,fa,fb,fc))
      }
      u<-c+gold*(c-b)
      fu<-f(u)
    } else if((c-u)*(u-ulim)>0){
      fu<-f(u)
      if(fu<fc){
        b<-c
        c<-u
        u<-c+gold*(c-b)
        fb<-fc
        fc<-fu
        fu<-f(u)
      }
    } else if((u-ulim)*(ulim-c)>=0){
      u<-ulim
      fu<-f(u)
    } else {
      u<-c+gold*(c-b)
      fu<-f(u)
    }
    a<-b
    b<-c
    c<-u
    fa<-fb
    fb<-fc
    fc<-fu
    cat(a,b,c,fa,fb,fc,"\n")
  }
  return(c(a,b,c,fa,fb,fc))
}

dbrent<-function(ax,bx,cx, f,df, tol,ITMAX=100){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  e<-0
  a<-min(ax,cx)
  b<-max(ax,cx)
  x<-w<-v<-bx
  fx<-fw<-fv<-f(x)
  dx<-dw<-dv<-df(x)
  iter=1
  while(iter<=ITMAX){
    xm=(a+b)/2
    tol1=tol*abs(x)+EPS
    tol2=2*tol1
    if(abs(x-xm)<=(tol2-(b-a)/2)){
      return(c(x,fx))
    }
    if(abs(e)>tol1){
      d1=2*(b-a)
      d2=d1
      if(dw != dx) d1=(w-x)*dx/(dx-dw) 
      if(dv != dx) d2=(v-x)*dx/(dx-dv) 
      u1=x+d1
      u2=x+d2
      ok1 = (a-u1)*(u1-b) > 0 & dx*d1 <= 0
      ok2 = (a-u2)*(u2-b) > 0 & dx*d2 <= 0
      olde=e 
      e=d
      if (ok1 | ok2) { 
      #Take only an acceptable d, and if both are acceptable, then take the smallest one.
        if (ok1 && ok2){
          d=ifelse(abs(d1)< abs(d2) , d1 , d2)
        } else{
          if(ok1){
            d=d1
          } else{
            d=d2
          }
        }
        if(abs(d)<=abs(olde/2)){
          u=x+d
          if(u-a<tol2|b-u<tol2) d=sign(xm-x)*abs(tol1)   
        } else{
          e=ifelse(dx>=0,a-x,b-x)
          d=e/2
        } 
      } else{
        e=ifelse(dx>=0,a-x,b-x)
        d=e/2
      }
    } else{
      e=ifelse(dx>=0,a-x,b-x)
      d=e/2
    }
    if(abs(d)>=tol1){
      u=x+d
      fu=f(u)
    } else{
      u=x+sign(d)*abs(tol1)
      fu=f(u)
      if(fu>fx){
        return(c(x,fx))
      }
    }
    du=df(u)
    if(fu<=fx){
      if(u>=x){
        a=x
      } else{
        b=x
      }
      v=w
      fv=fw
      dv=dw
      w=x
      fw=fx
      dw=dx
      x=u
      fx=fu
      dx=du
    } else{
      if(u<x){
        a=u
      } else{
        b=u
      }
      if(fu<=fw|w==x){
        v=w
        fv=fw
        dv=dw
        w=u
        fw=fu
        dw=du
      } else{
        if(fu<fv|v==x|v==w){
          v=u
          fv=fu
          dv=du
        }
      }
    }
  iter=iter+1
  }
  print("no converge in dlinmin")  
  return(c(x,fx))  
}

dlinmin<-function(p, xi , func ,dfunc,TOL=2e-4){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  f1dim<-function(x){
    xt<-p+x*xi
    func(xt)
  }
  df1dim<-function(x){
    xt<-p+x*xi
    df<-dfunc(xt)
    sum(df*xi)
  }  
  ax=0.0
  xx=1.0
  bracket<-mnbrak(f1dim,ax,xx)
  ax<-bracket[1]
  xx<-bracket[2]
  bx<-bracket[3]
  dbrent.result<-dbrent(ax,xx,bx, f1dim,df1dim, TOL)
  xmin<-dbrent.result[1]
  fmin<-dbrent.result[2]
  p<-p+xmin*xi
  return(list(x=p,f=fmin))
}

frprmin<-function(x, ftol, func,dfunc,ITMAX){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  f<-func(x)
  g<-dfunc(x)
  h<- -g
  iter<-1

  while(iter<=ITMAX){
    ln<-dlinmin(x,h,func,dfunc)
    fnew<-ln$f
    xnew<-ln$x
    if(2*abs(fnew-f)<=ftol*(abs(fnew)+abs(f)+EPS)) break
    f<-fnew
    gnew<-dfunc(xnew)
    gg<-sum(g^2)
    dgg<-sum((gnew-g)*gnew)
    g<-gnew
    if(gg==0) break
    gam<-dgg/gg
    h<-g+gam*h
    x<-xnew
    iter<-iter+1
  }
  list(x=xnew,f=fnew,iter=iter)
}


steplength<-function(x,p,func,dfunc,u1,u2){
  EPS<-1e-10 # Machine precision.
  TOLX<-4*EPS # Convergence criterion on x values.
  STPMX<-100.0 # Scaled maximum step length allowed in line searches.
  ALF<-1e-4 #Ensures sufficient decrease in function value.
  gold<-1.618034
  e0=.1
  e1=.2
  t0=0
  t1=Inf
  success=FALSE
  a=1
  tau<-function(t0,t1,q0,q0d,q1){
    if(t1==Inf){
      return(2*t0)
    }else{
      te=t0-q0d*(t1-t0)^2/(q1-q0-q0d*(t1-t0))/2
      if(t0==0){
        return(max(te,e0))
      }else{
        return (min(max(te,(1-e1)*t0+e1*t1),e1*t0+(1-e1)*t1))
      }
    }
  }   
  while(success==FALSE){
    if((func(x+a*p)-func(x))>(u1*a*sum(dfunc(x)*p))){
      t1=a
      a=tau(t0,t1,func(x-t0*p),sum(dfunc(x-t0*p)*p),func(x-t1*p))
    }else{
      if(sum(dfunc(x+a*p)*p)<(u2*sum(dfunc(x)*p))){
        t0=a
        a=tau(t0,t1,func(x-t0*p),sum(dfunc(x-t0*p)*p),func(x-t1*p))  
      }else{
        success=TRUE
      }
    }
    if(a==e0|a==e1) success=TRUE
  }
  a
}



