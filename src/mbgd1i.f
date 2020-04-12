C    loglik.odds.gradient(param,dataset,print.level=0)

      subroutine mbgd1i(gbeta,glpsi,npar,n)
      implicit double precision (a-h,o-z)
      DIMENSION x1(5000,10),theta1(5000),
     *work1(5000),y1(5000),beta1(10),bt1(10),
     *gbeta(10),dbeta(10),dbeta1(10),
     *tpr(2),der1(3),
     *P0(2,2), P1(2,2), P2(2,2),
     *dth(3),db(3,10),der(10)

      double precision lpsi1
      integer y1,k,k1,k2,i,i0,i1,j,npar,m,n0,n,mpar

      COMMON/param1/x1,theta1,work1,
     *y1,beta1,bt1,m,mpar,omega1,lpsi1
      COMMON/grad1/ dbeta,dbeta1,der,db

      psi = dexp(lpsi1)
      psi1 = psi
      ps1 = psi1-1
      call mati(x1,beta1,work1,5000,10,1,m,mpar+1)
      do 10 k=1,n
        theta1(k) = 1/(1+dexp(-work1(k)))
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

      p  = theta1(i0)
      do 40 k2=1,(npar-1)
        dbeta1(k2) =  p*(1-p)*x1(i0,k2)
        gbeta(k2)=(y1(i0)/(p*(1-p))-1/(1-p))*dbeta1(k2)    
   40  continue     
      glpsi = 0
      if (i0.eq.n0) return
      gpsi = 0

      i = i0+1
   50 if (i.le.n0) then
      i1=i
   60   if (y1(i1).eq.(-1)) then
      i1=i1+1
      go to 60
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
      if (i1.eq.i) then
      j  = y1(i0)
      th1 = theta1(i1)
      th2 = theta1(i1-1)
      call mcpj(th1,th2,psi1,tpr)
      p=tpr(y1(i0)+1)
      call deriv(theta1,psi,i1,j,der1)
          dpth  = der1(1)
          dpth1 = der1(2)
          dppsi = der1(3)
          dp = y1(i1)/(p*(1-p))-1/(1-p) 
          do 70 k2=1,(npar-1)
            dbeta(k2) = theta1(i1)*(1-theta1(i1))*x1(i,k2)
            dbeta1(k2) = theta1(i1-1)*(1-theta1(i1-1))*x1(i1-1,k2)
            gbeta(k2) = gbeta(k2) +
     *       dp*(dpth*dbeta(k2)+dpth1*dbeta1(k2))
   70  continue
          gpsi = gpsi+dp*dppsi
        else 
C    (exactly one intermediate missing datum between i0 and i1)
        call mcpj (theta1(i0+1),theta1(i0),psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)   
        call mcpj (theta1(i1),theta1(i0+1),psi1,tpr)
        call mat2 (tpr(1),tpr(2),P2)  
        call matp(P1,P2,P0,2,2,2)
      j  = y1(i0)
      tpr(1)= P0(1,2)
      tpr(2)= P0(2,2)
      p=tpr(y1(i0)+1)
      dp = y1(i1)/(p*(1-p))-1/(1-p) 
      do 80 k1=0,2
         do 90 k2=1,(npar-1)
           db(k1+1,k2)= theta1(k1+i0)*(1-theta1(k1+i0))*x1(k1+i0,k2)
   90  continue
   80  continue
      call deriv(theta1,psi1,i0+1,j,der1)
      dth(1)=der1(2)*(-P2(1,2)+P2(2,2))      
      dth(2)=der1(1)*(-P2(1,2)+P2(2,2))
      dpsi=der1(3)*(-P2(1,2)+P2(2,2))
      
      call deriv(theta1,psi1,i1,0,der1)
      dth(2)=dth(2)+der1(2)*(1-P1(j+1,2))
      dth(3)=der1(1)*(1-P1(j+1,2))
      dpsi=dpsi+der1(3)*(1-P1(j+1,2))

      call deriv(theta1,psi1,i1,1,der1)
      dth(2)=dth(2)+der1(2)*P1(j+1,2)
      dth(3)=dth(3)+der1(1)*P1(j+1,2)
      dpsi=dpsi+der1(3)*P1(j+1,2)

      do 100 k2=1,(npar-1)
         der(k2)=0.0D0
         do 110 k1=1,3
           der(k2)= der(k2)+dth(k1)*db(k1,k2)
  110  continue
  100  continue
      do 120 k2=1,(npar-1)
           gbeta(k2)= gbeta(k2)+dp*der(k2)
  120  continue
          gpsi = gpsi+dp*dpsi
      end if

        i0=i1
        i=i0+1
      go to 50
      end if
      glpsi=gpsi*psi
      return
      end

