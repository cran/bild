
C devivatives of p(t,t-1|j), i.e. prob{Y(t)=1|y(t-1)=j}
C wr to theta(t), theta(t-1), psi
C t=3
      subroutine deriv(theta,psi,t,j,der)
      implicit double precision (a-h,o-z)
      integer t,j,j1
      dimension theta(t),der(3)
      th = theta(t) 
      th1 = theta(t-1) 
      ps1 = psi-1 
      if(dabs(ps1) .gt. 1.0e-6) then
        j1 = 2*j-1
        delta=dsqrt(1+ps1*(psi*(th-th1)**2-(th+th1)**2+2*(th+th1)))
        dth=ps1*(psi*(th-th1)-(th+th1)+1)/delta
        dth1=ps1*(-psi*(th-th1)-(th+th1)+1)/delta
        dpsi=((2*psi-1)*(th-th1)**2-(th+th1)**2+2*(th+th1))/
     *  (2*delta)
        A=(2*ps1*(1-j+j1*th1))
        B=(1-delta+ps1*th1)*j1+th*ps1
        dpth=(-j1*dth+ps1)/A
        dpth1=(j1*(ps1-dth1)*A-2*ps1*j1*B)/A**2
        dppsi=((j1*(-dpsi+th1)+th)*A-2*B*(1-j+j1*th1))/A**2
        der(1)=dpth
        der(2)=dpth1
        der(3)=dppsi
      else
        der(1)=1
        der(2)=0
        der(3)=(th1-j)*(th**2-th) 
      end if
      return
      end  
