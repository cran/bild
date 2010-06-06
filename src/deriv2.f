C devivatives of pij=prob{Y(t)=1|y(t-2)=i,y(t-1)=j}
C wr to theta(t),theta(t-1), theta(t-2), psi1, psi2
C j=y(i0) (j=y_(t-2)), j1=y(i1)  (j1=y_(t-1))

      subroutine deriv2(theta,psi1,psi2,n,t,j,j1,der)
      implicit double precision (a-h,o-z)
      integer t,j,j1,j2,j3,j4
      dimension theta(n),der(5)
      th = theta(t)
      th1 = theta(t-1) 
      th2 = theta(t-2)
      ps1 = psi1-1 
      ps2 = psi2-1 

  
      if (dabs(ps1).gt.1.0e-6.and.dabs(ps2).gt.1.0e-6) then

      d0=dsqrt(1+ps1*(psi1*(th-th1)**2-(th+th1)**2+2*(th+th1)))
      d0th=(ps1*(psi1*(th-th1)-(th+th1)+1))/d0
      d0th1=(ps1*(-psi1*(th-th1)-(th+th1)+1))/d0
      d0th2=0
      d0psi1=((2*psi1-1)*(th-th1)**2-(th+th1)**2+2*(th+th1))/(2*d0)
      d0psi2=0   
      d=dsqrt(1+ps1*(psi1*(th1-th2)**2-(th1+th2)**2+ 2*(th1+th2)))
      dth=0
      dth1=ps1*(psi1*(th1-th2)-(th1+th2)+1)/d
      dth2=ps1*(-psi1*(th1-th2)-(th1+th2)+1)/d
      dpsi1=((2*psi1-1)*(th1-th2)**2-(th1+th2)**2+2*(th1+th2))/(2*d)
      dpsi2=0

      d1=dsqrt(th1**2+(ps2/(4*ps1**2))*(ps2*(d0-d+ps1*(th2-th))**2-
     *4*(d-1)*(d0-1)+ 4*ps1*((d-1)*th+(d0-1)*th2)+
     * 4*ps1**2*(th1**2-th*th2)))
      d1th=(1/(2*d1))*(1/(4*ps1**2))*(2*ps2**2*(d0th-ps1)*
     *(d0-d+ps1*(th2-th))+4*ps1**2*ps2*(-th2)+4*ps1*ps2*
     *   (d-1+d0th*th2)-4*ps2*(d0th*(d-1)))
      d1th1=(th1/d1)+(1/(2*d1))*(1/(4*ps1**2))*(2*ps2**2*(d0th1-dth1)*
     *  (d0-d+ps1*(th2-th))+ 4*ps1**2*ps2*(2*th1)+4*ps1*ps2*
     *  (dth1*th+d0th1*th2)-4*ps2*(dth1*(d0-1)+d0th1*(d-1)))
      d1th2=(1/(2*d1))*(1/(4*ps1**2))*(2*ps2**2*(-dth2+ps1)*
     *   (d0-d+ps1*(th2-th)) +4*ps1**2*ps2*(-th)+4*ps1*ps2*
     * (dth2*th+d0-1)-4*ps2*(dth2*(d0-1)))
      d1psi1=(1/(2*d1))*(1/(4*ps1**2))*(2*ps2**2*(d0psi1-dpsi1+th2-th)*
     *  (d0-d+ps1*(th2-th))+8*ps1*ps2*(th1**2-th*th2)+4*ps2*
     * ((d-1)*th+(d0-1)*th2)+4*ps1*ps2*(dpsi1*th+d0psi1*th2)-4*ps2*
     * (dpsi1*(d0-1)+d0psi1*(d-1)))-(1/(2*d1))*(d1**2-th1**2)*(2/ps1)
      d1psi2=(1/(2*d1))*(1/(4*ps1**2))*(2*ps2*(d0-d+ps1*(th2-th))**2+
     * 4*ps1**2*(th1**2-th*th2)+4*ps1*((d-1)*th+(d0-1)*th2)-
     * 4*(d-1)*(d0-1))
      
      d2=dsqrt((1-th1)**2+(ps2/(4*ps1**2))*(ps2*(d-d0+ps1*(th2-th))**2- 
     *4*(1-d)*(1-d0)+ 4*ps1*((d-1)*(1-th)+(d0-1)*(1-th2))+
     *  4*ps1**2*((1-th1)**2+(1-th)*(th2-1)) ) )
      d2th=(1/(2*d2))*(1/(4*ps1**2))*(2*ps2**2*(-d0th-ps1)*
     *(d-d0+ps1*(th2-th))+ 4*ps1**2*ps2*(1-th2)+
     *4*ps1*ps2*(1-d+d0th*(1-th2))- 4*ps2*(d0th*(d-1)) )      
      d2th1=((th1-1)/d2)+(1/(2*d2))*(1/(4*ps1**2))*(2*ps2**2*
     *  (dth1-d0th1)*(d-d0+ps1*(th2-th))+4*ps1**2*ps2*(2*(th1-1))+ 
     *  4*ps1*ps2*(dth1*(1-th)+d0th1*(1-th2))-
     * 4*ps2*(dth1*(d0-1)+d0th1*(d-1)))
      d2th2=(1/(2*d2))*(1/(4*ps1**2))*(2*ps2**2*(dth2+ps1)*
     *(d-d0+ps1*(th2-th))+ 4*ps1**2*ps2*(1-th)+4*ps1*ps2*
     *(dth2*(1-th)+1-d0)- 4*ps2*(dth2*(d0-1)) )
      d2psi1=(1/(2*d2))*(1/(4*ps1**2))*(2*ps2**2*(dpsi1-d0psi1+th2-th)*
     *   (d-d0+ps1*(th2-th))+8*ps1*ps2*((th1-1)**2+ (1-th)*(th2-1))+
     *   4*ps2*((d-1)*(1-th)+(d0-1)*(1-th2))+4*ps1*ps2*
     *(dpsi1*(1-th)+d0psi1*(1-th2))-   4*ps2*(dpsi1*(d0-1)+
     * d0psi1*(d-1)) )-(1/(2*d2))*(d2**2-(1-th1)**2)*(2/ps1)
      d2psi2=(1/(2*d2))*(1/(4*ps1**2))*(2*ps2*(d-d0+ps1*(th2-th))**2+
     * 4*ps1**2*((th1-1)**2+(1-th)*(th2-1))+ 4*ps1*((d-1)*(1-th)+
     * (d0-1)*(1-th2))-4*(1-d)*(1-d0))
C80      
      j2 = 2*j1-1
      j3 = 2*j-1
      j4 = 2*(j+j1-2*j*j1)-1
      A=2*ps2*(j4*(d-1)+ps1*(2-2*(j+j1-j*j1)+j2*th1+j3*th2))
      B=ps2*(2*j*j2-j2*d0+j4*d+ps1*(th+2*j*j2*th1+j3*th2))-
     *    2*j1*j3*ps1*(d1-th1)+2*(j1-1)*j3*ps1*(th1-1+d2)
       dpth=(1/A)*(ps2*(-j2*d0th+ps1)-2*j1*j3*ps1*d1th+
     * 2*(j1-1)*j3*ps1*d2th)
       dpth1=(1/A)*(ps2*(-j2*d0th1+j4*dth1+2*j*j2*ps1)-
     * 2*j1*j3*ps1*(d1th1-1)+2*(j1-1)*j3*ps1*(1+d2th1))- 
     *(B/A**2)*(2*ps2*(j4*dth1+j2*ps1))
       dpth2=(1/A)*(ps2*(j4*dth2+j3*ps1)-2*j1*j3*ps1*d1th2+
     * 2*(j1-1)*j3*ps1*d2th2)-(B/A**2)*(2*ps2*(j4*dth2+j3*ps1))      
      dppsi1=(1/A)*(ps2*(-j2*d0psi1+j4*dpsi1+(th+2*j*j2*th1+j3*th2))-
     *  2*j1*j3*(d1-th1)-2*j1*j3*ps1*d1psi1+2*(j1-1)*j3*(th1-1+d2)+
     *  2*(j1-1)*j3*ps1*d2psi1)-
     *  (B/A**2)*2*ps2*(j4*dpsi1+(2-2*(j+j1-j*j1)+j2*th1+j3*th2))
      dppsi2=(1/A)*(2*j*j2-j2*d0+j4*d+ps1*(th+2*j*j2*th1+j3*th2)-
     *  2*j1*j3*ps1*d1psi2+2*(j1-1)*j3*ps1*d2psi2)-
     *  (B/A**2)*2*(j4*(d-1)+ps1*(2-2*(j+j1-j*j1)+j2*th1+j3*th2))


      else if (dabs(ps1).gt.1.0e-6.and.dabs(ps2).lt.1.0e-6) then
     
      d0=dsqrt(1+ps1*(psi1*(th-th1)**2-(th+th1)**2+2*(th+th1)))
      d0th=(ps1*(psi1*(th-th1)-(th+th1)+1))/d0
      d0th1=(ps1*(-psi1*(th-th1)-(th+th1)+1))/d0
      d0th2=0
      d0psi1=((2*psi1-1)*(th-th1)**2-(th+th1)**2+2*(th+th1))/(2*d0)
      d0psi2=0   
      d=dsqrt(1+ps1*(psi1*(th1-th2)**2-(th1+th2)**2+ 2*(th1+th2)))
      dth=0
      dth1=ps1*(psi1*(th1-th2)-(th1+th2)+1)/d
      dth2=ps1*(-psi1*(th1-th2)-(th1+th2)+1)/d
      dpsi1=((2*psi1-1)*(th1-th2)**2-(th1+th2)**2+2*(th1+th2))/(2*d)
      dpsi2=0

      j2 = 2*j1-1
      j3 = 2*j-1
      j4 = 2*(j+j1-2*j*j1)-1
      A= 2*ps1*(1-j1+j2*th1)
      B=j2*(1-d0+ps1*th1)+ps1*th
C  d2d1=second derivative of delta1 in respect to psi2
      d2d1=(1/(512*ps1**6*th1**3))*(128*ps1**4*th1**2*
     *(d0-d+ps1*(th2-th))**2 -8*ps1**2*(4*ps1**2*(th1**2-
     *th*th2)+4*ps1*((d-1)*th+(d0-1)*th2)-
     *4*(1-d)*(1-d0))**2)
C  d2d2=second derivative of delta2 in respect to psi2
      d2d2=(1/(512*ps1**6*(1-th1)**3))*(128*ps1**4*(1-th1)**2*
     *(d-d0+ps1*(th2-th))**2 -8*ps1**2*(4*ps1**2*((th1-1)**2+
     *(1-th)*(th2-1))+4*ps1*((d-1)*(1-th)+(d0-1)*(1-th2))-
     *4*(1-d)*(1-d0))**2)
      dpth=(1/A)*(-j2*d0th+ps1)
      dpth1=(1/A)*(j2)*(-d0th1+ps1)-(B/A**2)*(j2*2*ps1)
      dpth2=0
      dppsi1=(1/A)*(j2*(-d0psi1+th1)+th)-
     *(B/A**2)*(2*(1-j1+j2*th1))
      dppsi2=(-j1*j3*2*ps1*d2d1+(j1-1)*j3*2*ps1*d2d2)/
     * (4*(j4*(d-1)+ps1*(2-2*(j+j1-j*j1)+j2*th1+j3*th2)))
C140

      else if (dabs(ps1).lt.1.0e-6.and.dabs(ps2).gt.1.0e-6) then
   
      d3=dsqrt(1+ps2**2*(th2-th)**2+ps2*(2*th+2*th2-4*th*th2))
      d3th=(-ps2**2*(th2-th)+ps2*(1-2*th2))/d3
      d3th1=0
      d3th2=(ps2**2*(th2-th)+ps2*(1-2*th))/d3
      d3psi1=0
      d3psi2=(ps2*(th2-th)**2+(th+th2-2*th*th2))/d3 

     
C  d1d=first derivative of delta in respect to psi1;  d1d0=first derivative of delta0 in respect to psi1
C  d2d=second derivative of delta in respect to psi1; d2d0=second derivative of delta0 in respect to psi1
C  d1d1=first derivative of delta1 in respect to psi1;  d1d2=first derivative of delta2 in respect to psi1

      d1d= -2*th1*th2+th1+th2
      d1d0=-2*th*th1+th+th1
      d2d=4*th1*th2*(1-th2)*(th1-1)
      d2d0=4*th*th1*(1-th1)*(th-1)

      d1=dsqrt(th1**2+(2*ps2**2*(d1d0-d1d+th2-th)**2+
     *    8*ps2*(th1**2-th*th2)+8*ps2*(d1d*th+d1d0*th2)-
     *    8*ps2*d1d*d1d0)/8)

      d2=dsqrt((1-th1)**2+(2*ps2**2*(d1d-d1d0+th2-th)**2+
     *    8*ps2*((1-th1)**2+(1-th)*(th2-1))+
     *    8*ps2*(d1d*(1-th)+d1d0*(1-th2))-8*ps2*d1d*d1d0)/8) 

      d1d1= (6*ps2**2*(d2d0-d2d)*(d1d0-d1d+th2-th)+
     *    12*ps2*(d2d*th+d2d0*th2)-
     *    12*ps2*(d2d*d1d0+d2d0*d1d) )/(48*d1)

      d1d2=(6*ps2**2*(d2d-d2d0)*(d1d-d1d0+th2-th)+
     *    12*ps2*(d2d*(1-th)+d2d0*(1-th2))-
     *    12*ps2*(d2d*d1d0+d2d0*d1d))/(48*d2)

      j2 = 2*j1-1
      j3 = 2*j-1
      j4 = 2*(j+j1-2*j*j1)-1
      A=2*ps2*(1-j+j3*th2)
      B=j3*(1-d3+ps2*th2)+ps2*th

      dpth=(1/A)*(-j3*d3th+ps2)
      dpth1=0
      dpth2=(1/A)*(j3)*(-d3th2+ps2)-(B/A**2)*(j3*2*ps2)
      dppsi1=0 
 
      dppsi1=(ps2*(-j2*d2d0+j4*d2d)-4*j1*j3*d1d1+
     * 4*(j1-1)*j3*d1d2)/
     * (4*ps2*(j4*d1d+2-2*(j+j1-j*j1)+j2*th1+j3*th2))-
     * ((ps2*(-j2*d1d0+j4*d1d+th+2*j*j2*th1+j3*th2)-
     * 2*j1*j3*(d1-th1)+2*(j1-1)*j3*(d2+th1-1))*2*j4*ps2*d2d)/
     * (2*(2*ps2*(j4*d1d+2-2*(j+j1-j*j1)+j2*th1+j3*th2))**2)
	
      dppsi2=(1/A)*(j3*(-d3psi2+th2)+th)-
     *(B/A**2)*(2*(1-j+j3*th2))

      else
C      if (ps1.eq.0.0.and.ps2.eq.0.0) then
C  # when psi(1)=1 and psi(2)=1      

      dpth=1
      dpth1=0
      dpth2=0
      dppsi1=(th1-j1)*(th**2-th)
      dppsi2=(th2-j)*(th**2-th)

      end if

       der(1)=dpth
       der(2)=dpth1
       der(3)=dpth2
       der(4)=dppsi1
       der(5)=dppsi2
      return
      end  
