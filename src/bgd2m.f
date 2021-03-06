
C    loglik.odds.gradient(param,dataset,print.level=0)

      subroutine bgd2m(gbeta,glpsi1,glpsi2,beta,lpsi,npar,x,y,
     *theta,work,dbeta,dbeta1,dbeta2,n,der,db,db1,db2)
      implicit double precision (a-h,o-z)
      dimension x(n,npar-2),beta(npar-2),work(n),theta(n),tpr(2),
     *tpr1(4),der1(5),der2(5),psi(2),gbeta(npar-2),dbeta(npar-2),
     *dbeta1(npar-2),dbeta2(npar-2),y(n),lpsi(2),
     *P0(2,2), P1(2,2), P2(2,2),P4(4,4), P5(4,4), P6(4,4),
     *Pc(4,1),Pr(2,1),Pr0(4,1),Pr1(2,1),
     *db(3,npar-2),db1(4,npar-2),db2(5,npar-2),der(npar-2),
     *P7(4,4), P8(4,4),dth(3),dth1(4),dth2(5)
      double precision lpsi
      integer y,k,k1,k2,i,i0,i1,i2,j,j1,j2,npar,n0,n

      psi(1) = dexp(lpsi(1))
      psi(2) = dexp(lpsi(2))
      psi1 = psi(1)
      psi2 = psi(2)   
      ps1 = psi1-1
      ps2 = psi2-1  
      call matp(x,beta,work,n,npar-2,1)
      do 10 k=1,n
        theta(k) = 1/(1+dexp(-work(k)))
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

      p  = theta(i0)
      do 40 k2=1,(npar-2)
        dbeta2(k2) =  p*(1-p)*x(i0,k2)
        gbeta(k2)=(y(i0)/(p*(1-p))-1/(1-p))*dbeta2(k2)    
   40  continue     
      glpsi1 = 0
      glpsi2 = 0
      if (i0.eq.n0) return
      gpsi1 = 0
      gpsi2 = 0

      i = i0+1
      i1 = i
   45   if (y(i1).eq.(-1)) then
      i1=i1+1
      go to 45
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
      if (i1.eq.i) then
      j  = y(i0)
      th1 = theta(i1)
      th2 = theta(i1-1)
      call mcpj(th1,th2,psi1,tpr)
      p=tpr(y(i0)+1)
      call deriv1(theta,psi1,psi2,n,i,j,der1)
          dpth  = der1(1)
          dpth1 = der1(2)
          dpth2 = der1(3)         
          dppsi1 = der1(4)
          dppsi2 = der1(5)
          dp = y(i1)/(p*(1-p))-1/(1-p) 
          do 50 k2=1,(npar-2)
            dbeta1(k2) = theta(i)*(1-theta(i))*x(i,k2)
            dbeta2(k2) = theta(i-1)*(1-theta(i-1))*x(i-1,k2)
            gbeta(k2) = gbeta(k2) +
     *       dp*(dpth1*dbeta1(k2)+dpth2*dbeta2(k2))
   50  continue
          gpsi1 = gpsi1+dp*dppsi1
          gpsi2 = gpsi2+dp*dppsi2
      
        else 
C    (exactly one intermediate missing datum between i0 and i1)
        th01 = theta(i0+1)
        th0  = theta(i0)
        th1 = theta(i1)
        call mcpj (th01,th0,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)   
        call mcpj (th1,th01,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P2)  
        call matp(P1,P2,P0,2,2,2)
      j  = y(i0)
      tpr(1)= P0(1,2)
      tpr(2)= P0(2,2)
      p=tpr(y(i0)+1)
      dp = y(i1)/(p*(1-p))-1/(1-p) 
      do 60 k1=0,2
         do 70 k2=1,(npar-2)
           db(k1+1,k2)= theta(k1+i0)*(1-theta(k1+i0))*x(k1+i0,k2)
   70  continue
   60  continue
      call deriv1(theta,psi1,psi2,n,i0+1,j,der1)
      dth(1)=der1(3)*(-P2(1,2)+P2(2,2))      
      dth(2)=der1(2)*(-P2(1,2)+P2(2,2))
      dpsi1=der1(4)*(-P2(1,2)+P2(2,2))
      dpsi2=der1(5)*(-P2(1,2)+P2(2,2))
      
      call deriv1(theta,psi1,psi2,n,i1,0,der1)
      dth(2)=dth(2)+der1(3)*(1-P1(j+1,2))
      dth(3)=der1(2)*(1-P1(j+1,2))
      dpsi1=dpsi1+der1(4)*(1-P1(j+1,2))
      dpsi2=dpsi2+der1(5)*(1-P1(j+1,2))

      call deriv1(theta,psi1,psi2,n,i1,1,der1)
      dth(2)=dth(2)+der1(3)*P1(j+1,2)
      dth(3)=dth(3)+der1(2)*P1(j+1,2)
      dpsi1=dpsi1+der1(4)*P1(j+1,2)
      dpsi2=dpsi2+der1(5)*P1(j+1,2)

      do 80 k2=1,(npar-2)
         der(k2)=0.0D0
         do 90 k1=1,3
           der(k2)= der(k2)+dth(k1)*db(k1,k2)
   90  continue
   80  continue
      do 95 k2=1,(npar-2)
           gbeta(k2)= gbeta(k2)+dp*der(k2)
   95  continue

          gpsi1 = gpsi1+dp*dpsi1
          gpsi2 = gpsi2+dp*dpsi2
      end if

      i=i1+1
  100 if (i.le.n0) then
      i2 = i
  110   if (y(i2).eq.(-1)) then
      i2=i2+1
      go to 110
      end if

      if (i2.eq.i.and.i1.eq.(i0+1)) then
        j  = y(i0)
        j1  = y(i1)
        th = theta(i2)
        th1 = theta(i2-1)
        th2 = theta(i2-2)
        call mcpij(th,th1,th2,psi1,psi2,tpr1)
        call deriv2(theta,psi1,psi2,n,i2,j,j1,der2)
          dpth  = der2(1)
          dpth1 = der2(2)
          dpth2 = der2(3)         
          dppsi1 = der2(4)
          dppsi2 = der2(5)
          p=tpr1(y(i0)+2*y(i1)+1)
          dp = y(i2)/(p*(1-p))-1/(1-p) 
          do 120 k2=1,(npar-2)
            dbeta (k2) = theta(i2)*(1-theta(i2))*x(i2,k2) 
            dbeta1(k2) = theta(i2-1)*(1-theta(i2-1))*x(i2-1,k2)
            dbeta2(k2) = theta(i2-2)*(1-theta(i2-2))*x(i2-2,k2)
            gbeta(k2) = gbeta(k2) +
     *      dp*(dpth*dbeta(k2)+dpth1*dbeta1(k2)+dpth2*dbeta2(k2))
  120   continue
          gpsi1 = gpsi1+dp*dppsi1
          gpsi2 = gpsi2+dp*dppsi2
        i0=i1
        i1=i2
        i=i1+1      

      else if (i1.ne.(i0+1).and.i2.eq.(i1+1)) then
C    (exactly one intermediate missing datum between i0 and i1)
        th01 = theta(i0+1)
        th0 = theta(i0)
        th = theta(i2)
        th1 = theta(i2-1)
        th2 = theta(i2-2)
        call mcpj (th01,th0,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)
        call mcpij (th,th1,th2,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P6)
        Pr(1,1)=tpr1(2*y(i1)+1)
        Pr(2,1)=tpr1(2*y(i1)+2)
        call matp(P1,Pr,Pr1,2,2,1)
        p=Pr1(y(i0)+1,1)
       dp = y(i2)/(p*(1-p))-1/(1-p) 

        j  = y(i0)
        j1  = y(i1)
      do 130 k1=0,3
         do 140 k2=1,(npar-2)
           db1(k1+1,k2)= theta(k1+i0)*(1-theta(k1+i0))*x(k1+i0,k2)
  140  continue
  130  continue

      call deriv1(theta,psi1,psi2,n,i0+1,j,der1)
      dth1(1)=der1(3)*(-P6(j1+1,2*j1+2)+P6(j1+3,2*j1+2))
      dth1(2)=der1(2)*(-P6(j1+1,2*j1+2)+P6(j1+3,2*j1+2))    
      dpsi1=der1(4)*(-P6(j1+1,2*j1+2)+P6(j1+3,2*j1+2))       

      call deriv2(theta,psi1,psi2,n,i2,0,j1,der2)
      dth1(2)=dth1(2)+der2(3)*(1-P1(j+1,2))     
      dth1(3)=der2(2)*(1-P1(j+1,2))     
      dth1(4)=der2(1)*(1-P1(j+1,2))     
      dpsi1=dpsi1+der2(4)*(1-P1(j+1,2))     
      dpsi2=der2(5)*(1-P1(j+1,2))     

      call deriv2(theta,psi1,psi2,n,i2,1,j1,der2)
      dth1(2)=dth1(2)+der2(3)*P1(j+1,2)
      dth1(3)=dth1(3)+der2(2)*P1(j+1,2)
      dth1(4)=dth1(4)+der2(1)*P1(j+1,2)
      dpsi1=dpsi1+der2(4)*P1(j+1,2)
      dpsi2=dpsi2+der2(5)*P1(j+1,2) 

      do 150 k2=1,(npar-2)
         der(k2)=0.0D0
         do 160 k1=1,4
           der(k2)= der(k2)+dth1(k1)*db1(k1,k2)
  160  continue
  150  continue

      do 170 k2=1,(npar-2)
           gbeta(k2)= gbeta(k2)+dp*der(k2)
  170  continue
          gpsi1 = gpsi1+dp*dpsi1
          gpsi2 = gpsi2+dp*dpsi2
        i0=i1
        i1=i2
        i=i1+1   

      else if (i2.ne.(i1+1).and.i2.ne.n0) then
C    (exactly one intermediate missing datum between i1 and i2)
        th02 = theta(i0+2)
	      th01 = theta(i0+1)        
	      th0 = theta(i0)
        th12 = theta(i1+2)
        th11 = theta(i1+1)
        th22 = theta(i1)
 	call mcpij (th02,th01,th0,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P5)
        call mcpij (th12,th11,th22,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P6)
        call matp(P5,P6,P4,4,4,4)

        i3=i2+1
      	th = theta(i3)
        th1 = theta(i3-1)
        th2 = theta(i3-2)
        call mcpij (th,th1,th2,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P7)
        call matp(P4,P7,P8,4,4,4)
        pstar=P8(2*y(i0)+y(i1)+1,2*y(i2)+y(i3)+1)
       dp = 1/pstar

        j  = y(i0)
        j1  = y(i1)
        j2  = y(i2)
      do 230 k1=0,4
         do 240 k2=1,(npar-2)
           db2(k1+1,k2)= theta(k1+i0)*(1-theta(k1+i0))*x(k1+i0,k2)
  240  continue
  230  continue

       call deriv2(theta,psi1,psi2,n,i0+2,j,j1,der2)
      dth2(1)=der2(3)*((P6(2*j1+1,2)-1+(1-2*P6(2*j1+1,2))*j2)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3))+
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*j2)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dth2(2)=der2(2)*((P6(2*j1+1,2)-1+(1-2*P6(2*j1+1,2))*j2)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3))+
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*j2)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dth2(3)=der2(1)*((P6(2*j1+1,2)-1+(1-2*P6(2*j1+1,2))*j2)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3))+
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*j2)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dpsi1=der2(4)*((P6(2*j1+1,2)-1+(1-2*P6(2*j1+1,2))*j2)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3))+
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*j2)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dpsi2=der2(5)*((P6(2*j1+1,2)-1+(1-2*P6(2*j1+1,2))*j2)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3))+
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*j2)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))

      call deriv2(theta,psi1,psi2,n,i2,j1,0,der2)
      dth2(2)=dth2(2)+der2(3)*((1-P5(2*j+j1+1,2*j1+2))*(2*j2-1)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3)))
      dth2(3)=dth2(3)+der2(2)*((1-P5(2*j+j1+1,2*j1+2))*(2*j2-1)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3)))
      dth2(4)=der2(1)*((1-P5(2*j+j1+1,2*j1+2))*(2*j2-1)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3)))
      dpsi1=dpsi1+der2(4)*((1-P5(2*j+j1+1,2*j1+2))*(2*j2-1)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3)))
      dpsi2=dpsi2+der2(5)*((1-P5(2*j+j1+1,2*j1+2))*(2*j2-1)*
     *(1-P7(j2+1,2*j2+2)+(2*P7(j2+1,2*j2+2)-1)*y(i3)))

      call deriv2(theta,psi1,psi2,n,i2,j1,1,der2)
      dth2(2)=dth2(2)+der2(3)*(P5(2*j+j1+1,2*j1+2)*(2*j2-1)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dth2(3)=dth2(3)+der2(2)*(P5(2*j+j1+1,2*j1+2)*(2*j2-1)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dth2(4)=dth2(4)+der2(1)*(P5(2*j+j1+1,2*j1+2)*(2*j2-1)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dpsi1=dpsi1+der2(4)*(P5(2*j+j1+1,2*j1+2)*(2*j2-1)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))
      dpsi2=dpsi2+der2(5)*(P5(2*j+j1+1,2*j1+2)*(2*j2-1)*
     *(1-P7(j2+3,2*j2+2)+(2*P7(j2+3,2*j2+2)-1)*y(i3)))

      call deriv2(theta,psi1,psi2,n,i2+1,0,1,der2)
      dth2(3)=dth2(3)+der2(3)*((1-P5(2*j+j1+1,2*j1+2))*(2*y(i3)-1)*
     *(1-P6(2*j1+1,2)+(2*P6(2*j1+1,2)-1)*y(i2)))
      dth2(4)=dth2(4)+der2(2)*((1-P5(2*j+j1+1,2*j1+2))*(2*y(i3)-1)*
     *(1-P6(2*j1+1,2)+(2*P6(2*j1+1,2)-1)*y(i2)))
      dth2(5)=der2(1)*((1-P5(2*j+j1+1,2*j1+2))*(2*y(i3)-1)*
     *(1-P6(2*j1+1,2)+(2*P6(2*j1+1,2)-1)*y(i2)))
      dpsi1=dpsi1+der2(4)*((1-P5(2*j+j1+1,2*j1+2))*(2*y(i3)-1)*
     *(1-P6(2*j1+1,2)+(2*P6(2*j1+1,2)-1)*y(i2)))
      dpsi2=dpsi2+der2(5)*((1-P5(2*j+j1+1,2*j1+2))*(2*y(i3)-1)*
     *(1-P6(2*j1+1,2)+(2*P6(2*j1+1,2)-1)*y(i2)))

      call deriv2(theta,psi1,psi2,n,i2+1,1,1,der2)
      dth2(3)=dth2(3)+der2(3)*(P5(2*j+j1+1,2*j1+2)*(2*y(i3)-1)*
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*y(i2)))
      dth2(4)=dth2(4)+der2(2)*(P5(2*j+j1+1,2*j1+2)*(2*y(i3)-1)*
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*y(i2)))
      dth2(5)=dth2(5)+der2(1)*(P5(2*j+j1+1,2*j1+2)*(2*y(i3)-1)*
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*y(i2)))
      dpsi1=dpsi1+der2(4)*(P5(2*j+j1+1,2*j1+2)*(2*y(i3)-1)*
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*y(i2)))
      dpsi2=dpsi2+der2(5)*(P5(2*j+j1+1,2*j1+2)*(2*y(i3)-1)*
     *(1-P6(2*j1+2,4)+(2*P6(2*j1+2,4)-1)*y(i2)))

      do 250 k2=1,(npar-2)
         der(k2)=0.0D0
         do 260 k1=1,5
           der(k2)= der(k2)+dth2(k1)*db2(k1,k2)
  260  continue
  250  continue

      do 270 k2=1,(npar-2)
           gbeta(k2)= gbeta(k2)+dp*der(k2)
  270  continue

          gpsi1 = gpsi1+dp*dpsi1
          gpsi2 = gpsi2+dp*dpsi2

        i0=i2
        i1=i3
        i=i1+1

      else if (i2.ne.(i1+1).and.i2.eq.n0) then
C    (exactly one intermediate missing datum between i1 and i2(last value))
        th02 = theta(i0+2)
        th01= theta(i0+1)
        th0 = theta(i0)
     	th12 = theta(i1+2)
        th11 = theta(i1+1)
        th1 = theta(i1)
	call mcpij (th02,th01,th0,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P5)
        call mcpij (th12,th11,th1,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P6)
        call matp(P5,P6,P4,4,4,4)
        Pc(1,1)=0
        Pc(2,1)=1
        Pc(3,1)=0
        Pc(4,1)=1
        call matp(P4,Pc,Pr0,4,4,1)
        p=Pr0(2*y(i0)+y(i1)+1,1)
       dp = y(i2)/(p*(1-p))-1/(1-p) 

        j  = y(i0)
        j1  = y(i1)
      do 330 k1=0,3
         do 340 k2=1,(npar-2)
           db1(k1+1,k2)= theta(k1+i0)*(1-theta(k1+i0))*x(k1+i0,k2)
  340  continue
  330  continue

      call deriv2(theta,psi1,psi2,n,i0+2,j,j1,der2)
      dth1(1)=der2(3)*(-P6(2*j1+1,2)+P6(2*j1+2,4))      
      dth1(2)=der2(2)*(-P6(2*j1+1,2)+P6(2*j1+2,4))      
      dth1(3)=der2(1)*(-P6(2*j1+1,2)+P6(2*j1+2,4))    
      dpsi1=der2(4)*(-P6(2*j1+1,2)+P6(2*j1+2,4))    
      dpsi2=der2(5)*(-P6(2*j1+1,2)+P6(2*j1+2,4))    
      
      call deriv2(theta,psi1,psi2,n,i2,j1,0,der2)
      dth1(2)=dth1(2)+der2(3)*(1-P5(2*j+j1+1,2*j1+2))      
      dth1(3)=dth1(3)+der2(2)*(1-P5(2*j+j1+1,2*j1+2))  
      dth1(4)=der2(1)*(1-P5(2*j+j1+1,2*j1+2))  
      dpsi1=dpsi1+der2(4)*(1-P5(2*j+j1+1,2*j1+2))  
      dpsi2=dpsi2+der2(5)*(1-P5(2*j+j1+1,2*j1+2))  

      call deriv2(theta,psi1,psi2,n,i2,j1,1,der2)
      dth1(2)=dth1(2)+der2(3)*P5(2*j+j1+1,2*j1+2)     
      dth1(3)=dth1(3)+der2(2)*P5(2*j+j1+1,2*j1+2)     
      dth1(4)=dth1(4)+der2(1)*P5(2*j+j1+1,2*j1+2)     
      dpsi1=dpsi1+der2(4)*P5(2*j+j1+1,2*j1+2)     
      dpsi2=dpsi2+der2(5)*P5(2*j+j1+1,2*j1+2)  

      do 350 k2=1,(npar-2)
         der(k2)=0.0D0
         do 360 k1=1,4
           der(k2)= der(k2)+dth1(k1)*db1(k1,k2)
  360  continue
  350  continue

      do 370 k2=1,(npar-2)
           gbeta(k2)= gbeta(k2)+dp*der(k2)
  370  continue

          gpsi1 = gpsi1+dp*dpsi1
          gpsi2 = gpsi2+dp*dpsi2
        i0=i1
        i1=i2
        i=i1+1   
       end if

      go to 100
      end if

      glpsi1=gpsi1*psi1
      glpsi2=gpsi2*psi2

      return
      end

