C # Calculates the log-likelihood for the first
C# order dependence model with random effects

      subroutine mlik1i(logL,pij,npar,n)
      implicit double precision (a-h,o-z)
      DIMENSION x1(5000,10),theta1(5000),
     *work1(5000),y1(5000),beta1(10),bt1(10),
     *pij(n),tpr(2),P0(2,2), P1(2,2), P2(2,2)
      double precision logL,lpsi,pij
      integer y1,i0,i1,npar,n,n0,k,m,mpar

      COMMON/param/x1,theta1,work1,
     *y1,beta1,bt1,m,mpar,omega1,lpsi

      psi = dexp(lpsi)
      psi1 = psi
      ps1 = psi1-1
      call mati(x1,beta1,work1,5000,10,1,n,npar+1)

      do 10 i=1,n
        theta1(i) = 1/(1+dexp(-work1(i)))
   10 continue
      i0 = 1  
   20   if (y1(i0).eq.(-1)) then
      i0=i0+1
      go to 20
      end if

      n0 = n
      go to 30
   30   if (y1(n0).eq.(-1)) then
      n0=n0-1
      go to 30
      end if

      logL = 0
      p  = theta1(i0)
      pij(i0)=p
      logL = y1(i0)*dlog(p/(1-p))+dlog(1-p)
      if (i0.eq.n0) return
      i = i0+1
   40 if (i.le.n0) then
      i1 = i
   50   if (y1(i1).eq.(-1)) then
      i1=i1+1
      go to 50
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
      if (i1.eq.i) then
      th1 = theta1(i)
      th2 = theta1(i-1)
      call mcpj(th1,th2,psi1,tpr)
      else 
C    (one or more intermediate missing datum)
      call mat2 (0.0D0,1.0D0,P0)
      do 60 k=(i0+1),i1
        th1 = theta1(k)
        th2 = theta1(k-1)
        call mcpj (th1,th2,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)   
        call matp(P0,P1,P2,2,2,2)
        call matc(P2,P0,2,2)
   60 continue
      tpr(1)= P0(1,2)
      tpr(2)= P0(2,2)
      end if
      p=tpr(y1(i0)+1)
      pij(i1)=p
      logL = logL+y1(i1)*dlog(p/(1-p))+dlog(1-p)
        i0=i1
        i=i0+1
      go to 40
      end if
      return
      end
