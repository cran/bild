      subroutine integ(logL,bt2,beta2,lpsi,omega,npar,
     *x2,y2,theta2,work2,n,li,ls,epsabs,epsrel,key,limit)
      implicit double precision (a-h,o-z)    
      External f2, f2pij

      DOUBLE PRECISION logL, lpsi, lpsi1,li,ls,
     *sigma,sqr2pi
      INTEGER n,npar,y2,key,limit,neval,ier,iord,last,
     *m,mpar,y1,k,k1,k2

      DIMENSION x1(4500,10),theta1(4500),work1(4500),
     *y1(4500),lpsi1(2),beta1(10),bt1(10),
     *x2(n,npar-2),beta2(npar-2),
     *theta2(n),lpsi(2),y2(n),work2(n),bt2(npar-2)

      COMMON/param/x1,theta1,work1,
     *y1,lpsi1,beta1,bt1,m,mpar,omega1

      DATA sqr2pi/2.506628274631D0/

      do 10 k=1,(npar-2)
      bt1(k)=bt2(k)
      beta1(k)=beta2(k)
   10 continue 
      do 20 k1=1,2
      lpsi1(k1)=lpsi(k1)
   20 continue
      do 30 k2=1,n 
          do 40 k=1,(npar-2)
           x1(k2,k)=x2(k2,k)
   40 continue 
      y1(k2)=y2(k2)
      theta1(k2)=theta2(k2)
      work1(k2)=work2(k2)
   30 continue 
      m=n
      mpar=npar 
      omega1=omega
      sigma=dexp(omega/(2.0d0))
      a=li*sigma
      b=ls*sigma

      CALL dqager(f2,a,b,epsabs,epsrel,key,limit,result,abserr,
     *neval,ier,alist,blist,rlist,elist,iord,last,1)
      logL=dlog(result/(sqr2pi*sigma))

      RETURN
      END
      
