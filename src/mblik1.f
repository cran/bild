      subroutine mblik1(logL,pij,beta,lpsi,npar,x,y,theta,work,n)
      implicit double precision (a-h,o-z)
      dimension x(n,npar-1),beta(npar-1),work(n),theta(n),pij(n),
     *tpr(2),y(n),P0(2,2), P1(2,2), P2(2,2)
      double precision logL,lpsi,pij
      integer y,i0,i1,npar,n,n0,k
      psi= dexp(lpsi)
      psi1 = psi
      ps1 = psi1-1
      call matp(x,beta,work,n,npar-1,1)
      do 10 i=1,n
        theta(i) = 1/(1+dexp(-work(i)))
   10 continue

      i0 = 1  
   20   if (y(i0).eq.(-1)) then
      i0=i0+1
      go to 20
      end if

      n0 = n
   30   if (y(n0).eq.(-1)) then
      n0=n0-1
      go to 30
      end if

        logL = 0
      p  = theta(i0)
      pij(i0)=p
      logL = y(i0)*dlog(p/(1-p))+dlog(1-p)
      if (i0.eq.n0) return
      i = i0+1
   40 if (i.le.n0) then
      i1 = i
   50   if (y(i1).eq.(-1)) then
      i1=i1+1
      go to 50
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
      if (i1.eq.i) then
      th1 = theta(i)
      th2 = theta(i-1)
      call mcpj(th1,th2,psi1,tpr)
      else 
C    (one or more intermediate missing datum)
      call mat2 (0.0D0,1.0D0,P0)
      do 60 k=(i0+1),i1
        th1 = theta(k)
        th2 = theta(k-1)
        call mcpj (th1,th2,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)   
        call matp(P0,P1,P2,2,2,2)
        call matc(P2,P0,2,2)
   60 continue
      tpr(1)= P0(1,2)
      tpr(2)= P0(2,2)
      end if
      p=tpr(y(i0)+1)
      pij(i1)=p
      logL = logL+y(i1)*dlog(p/(1-p))+dlog(1-p)
        i0=i1
        i=i0+1
      go to 40
      end if
      return
      end
