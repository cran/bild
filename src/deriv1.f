
C devivatives of p(t-1,t-2|j), i.e. prob{Y(t-1)=1|y(t-2)=j}
C wr to theta(t-1), theta(t-2), psi1
C t=3
      subroutine deriv1(theta,psi1,psi2,n,t,j,der)
      implicit double precision (a-h,o-z)
      integer t,j,j1
      dimension theta(n),der(5)
      th1 = theta(t) 
      th2 = theta(t-1) 
      ps1 = psi1-1 
      ps2 = psi2-1 
      if(dabs(ps1) .gt. 1.0e-6) then
        j1 = 2*j-1
        delta=dsqrt(1+ps1*(psi1*(th1-th2)**2-(th1+th2)**2+2*(th1+th2)))
        dth=0
        dth1=ps1*(psi1*(th1-th2)-(th1+th2)+1)/delta
        dth2=ps1*(-psi1*(th1-th2)-(th1+th2)+1)/delta
        dpsi1=((2*psi1-1)*(th1-th2)**2-(th1+th2)**2+2*(th1+th2))/
     *  (2*delta)
        dpsi2=0
        A=(2*ps1*(1-j+j1*th2))
        B=(1-delta+ps1*th2)*j1+th1*ps1
        dpth=0
        dpth1=(-j1*dth1+ps1)/A
        dpth2=(j1*(ps1-dth2)*A-2*ps1*j1*B)/A**2
        dppsi1=((j1*(-dpsi1+th2)+th1)*A-2*B*(1-j+j1*th2))/A**2
        dppsi2=0
        der(1)=dpth
        der(2)=dpth1
        der(3)=dpth2
        der(4)=dppsi1
        der(5)=dppsi2
      else
        der(1)=0
        der(2)=1
        der(3)=0
        der(4)=(th2-j)*(th1**2-th1) 
        der(5)=0
      end if
      return
      end  
