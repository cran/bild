      subroutine gint(gbeta,glpsi1,glpsi2,gvar,
     *x2,theta2,work2,y2,lpsi,beta2,bt2,
     *dbt,dbt1,dbt2,dder,ddb,ddb1,ddb2,
     *omega,npar,n,li,ls,epsabs,epsrel,key,limit)

      implicit double precision (a-h,o-z)
      External f2b,fps1,fps2,f2var
      DOUBLE PRECISION lpsi,lpsi1,li,ls
      INTEGER n,npar,y2,key,limit,neval,ier,iord,last,
     *m,mpar,y1,k,k1,k2,iaux,kaux

      DIMENSION x1(4500,10),theta1(4500),work1(4500),
     *y1(4500),lpsi1(2),beta1(10),bt1(10),
     *dbeta(10),dbeta1(10),dbeta2(10),
     *der(10),db(3,10),db1(4,10),db2(5,10),
     *x2(n,npar-2),theta2(n),work2(n),
     *y2(n),lpsi(2),beta2(npar-2),bt2(npar-2),
     *dbt(npar-2),dbt1(npar-2),dbt2(npar-2),
     *dder(npar-2),ddb(3,npar-2),ddb1(4,npar-2),
     *ddb2(5,npar-2),gbeta(npar-2),
     *alist(limit),blist(limit),elist(limit),iord(limit),
     *rlist(limit)

      COMMON/param/x1,theta1,work1,
     *y1,lpsi1,beta1,bt1,m,mpar,omega1
      COMMON/grad/ dbeta,dbeta1,dbeta2,
     *der,db,db1,db2

      do 10 k=1,(npar-2)
      bt1(k)=bt2(k)
      beta1(k)=beta2(k)
      dbeta(k)=dbt(k)
      dbeta1(k)=dbt1(k)
      dbeta2(k)=dbt2(k)
      der(k)=dder(k)
      do15 kaux=1,3
      db(kaux,k)=ddb(kaux,k)
      db1(kaux,k)=ddb1(kaux,k)
      db2(kaux,k)=ddb2(kaux,k)
   15 continue
      db1(4,k)=ddb1(4,k)
      db2(4,k)=ddb2(4,k)
      db2(5,k)=ddb2(5,k)
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

      a=li*dexp(omega/(2.0d0))
      b=ls*dexp(omega/(2.0d0))

      do 50 iaux=1,(npar-2)
      CALL dqager(f2b,a,b,epsabs,epsrel,key,limit,res1,abserr,
     *neval,ier,alist,blist,rlist,elist,iord,last,iaux)
      gbeta(iaux)=res1
   50  continue   

      CALL dqager(fps1,a,b,epsabs,epsrel,key,limit,res3,abserr,
     * neval,ier,alist,blist,rlist,elist,iord,last,1)
       glpsi1=res3
      CALL dqager(fps2,a,b,epsabs,epsrel,key,limit,res4,abserr,
     * neval,ier,alist,blist,rlist,elist,iord,last,1)
      glpsi2=res4
      CALL dqager(f2var,a,b,epsabs,epsrel,key,limit,res5,abserr,
     * neval,ier,alist,blist,rlist,elist,iord,last,1)
      gvar=res5

      RETURN
      END
      
