
C#######calculate the probabilities of pij depending on psi1 and psi2    

      subroutine mcpij(th,th1,th2,psi1,psi2,p)
      implicit double precision(a-h,o-z)
      dimension p(4)

      ps1=psi1-1
      ps2=psi2-1

C----- p=(p00,p10,p01,p11)

      if (dabs(ps1).gt.1.0e-10.and.dabs(ps2).gt.1.0e-10) then

      d0=dsqrt(1+ps1*(psi1*(th-th1)**2-(th+th1)**2+2*(th+th1)))

      d=dsqrt(1+ps1*(psi1*(th1-th2)**2-(th1+th2)**2+2*(th1+th2)))

      d1=dsqrt(th1**2+(ps2/(4*ps1**2))*(ps2*(d0-d+ps1*(th2-th))**2-
     *  4*(d-1)*(d0-1)+4*ps1*((d-1)*th+(d0-1)*th2)+
     *  4*ps1**2*(th1**2-th*th2)))

      d2=dsqrt((1-th1)**2+(ps2/(4*ps1**2))*(ps2*(d-d0+ps1*(th2-th))**2-
     *   4*(1-d)*(1-d0)+4*ps1*((d-1)*(1-th)+(d0-1)*(1-th2))+
     *   4*ps1**2*((1-th1)**2+(1-th)*(th2-1)) ) )

	p(1)=(ps2*((d0-d)+ps1*(th-th2))+2*ps1*(th1-1+d2))/
     *		(2*ps2*(1-d+ps1*(2-th1-th2)))

	p(2)=(ps2*((d+d0-2)+ps1*(th-2*th1+th2))+2*ps1*(1-th1-d2))/ 
     *		(2*ps2*(d-1+ps1*(th2-th1)))

	p(3)=(ps2*(d-d0+ps1*(th-th2))+2*ps1*(d1-th1))/
     *		(2*ps2*(d-1+ps1*(th1-th2)))
	
	p(4)=(ps2*(2-d0-d+ps1*(th+2*th1+th2))-2*ps1*(d1-th1))/
     *	     	(2*ps2*(1-d+ps1*(th1+th2)))

      else if (dabs(ps1).gt.1.0e-10.and.dabs(ps2).lt.1.0e-10) then
C      else if (dabs(ps1).gt.1.0e-10.and.ps2.eq.0.0) then

      d0=dsqrt(1+ps1*(psi1*(th-th1)**2-(th+th1)**2+2*(th+th1)))

	p(1)=((d0-1)+ps1*(th-th1))/(2*ps1*(1-th1))

	p(2)=((d0-1)+ps1*(th-th1))/(2*ps1*(1-th1))
	
	p(3)=((1-d0)+ps1*(th+th1))/(2*ps1*th1)

	p(4)=((1-d0)+ps1*(th+th1))/(2*ps1*th1)


      else if (dabs(ps1).lt.1.0e-10.and.dabs(ps2).gt.1.0e-10) then
C      else if (ps1.eq.0.0.and.dabs(ps2).gt.1.0e-10) then

      	d3=dsqrt(1+ps2**2*(th2-th)**2+ps2*(2*th+2*th2-4*th*th2))

	p(1)=((d3-1)+ps2*(th-th2))/(2*ps2*(1-th2))

	p(2)=((1-d3)+ps2*(th+th2))/(2*ps2*th2)

	p(3)=((d3-1)+ps2*(th-th2))/(2*ps2*(1-th2))

	p(4)=((1-d3)+ps2*(th+th2))/(2*ps2*th2)
       else

      p(1)=th
      p(2)=th
      p(3)=th
      p(4)=th

      end if
      return
      end
