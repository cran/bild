      subroutine gint1(gbeta,glpsi1,gvar,bt2,
     *beta2,lpsi,omega,npar,x2,y2,theta2,work2,n,
     *dbt,dbt1,dder,ddb,li,ls,epsabs,epsrel,key,limit)

      implicit double precision (a-h,o-z)
      External f1b,fpsi,f1var
      DOUBLE PRECISION lpsi,lpsi1,li,ls
      INTEGER n,npar,y2,key,limit,neval,ier,iord,last,
     *m,mpar,y1,k,k1,k2,iaux,kaux

      DIMENSION x1(5000,10),theta1(5000),
     *work1(5000),y1(5000),beta1(10),
     *bt1(10),dbeta(10),dbeta1(10),der(10),db(3,10),
     * x2(n,npar-1),beta2(npar-1),
     *theta2(n),y2(n),work2(n),bt2(npar-1),
     *gbeta(npar-1),dbt(npar-1),dbt1(npar-1),
     *dder(npar-1),ddb(3,npar-1)

      COMMON/param/x1,theta1,work1,
     *y1,beta1,bt1,m,mpar,omega1,lpsi1
      COMMON/grad/ dbeta,dbeta1,der,db

      do 10 k=1,(npar-1)
      bt1(k)=bt2(k)
      beta1(k)=beta2(k)
      dbeta(k)=dbt(k)
      dbeta1(k)=dbt1(k)
      der(k)=dder(k)
      do 15 kaux=1,3
      db(kaux,k)=ddb(kaux,k)
   15 continue 
   10 continue 
      do 30 k2=1,n 
          do 40 k=1,(npar-1)
           x1(k2,k)=x2(k2,k)
   40 continue 
      y1(k2)=y2(k2)
      theta1(k2)=theta2(k2)
      work1(k2)=work2(k2)
   30 continue 
      m=n
      mpar=npar 
      omega1=omega
      lpsi1=lpsi
      a=li*dexp(omega/(2.0d0))
      b=ls*dexp(omega/(2.0d0))

      do 50 iaux=1,(npar-1)
      CALL dqager(f1b,a,b,epsabs,epsrel,key,limit,res1,abserr,
     *neval,ier,alist,blist,rlist,elist,iord,last,iaux)
      gbeta(iaux)=res1
   50  continue   

      CALL dqager(fpsi,a,b,epsabs,epsrel,key,limit,res3,abserr,
     * neval,ier,alist,blist,rlist,elist,iord,last,1)
       glpsi1=res3
      CALL dqager(f1var,a,b,epsabs,epsrel,key,limit,res5,abserr,
     * neval,ier,alist,blist,rlist,elist,iord,last,1)
      gvar=res5
      RETURN
      END
      
