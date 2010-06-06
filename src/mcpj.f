
      subroutine mcpj(th1,th2,psi1,p)
      implicit double precision(a-h,o-z)
      dimension p(2)

      ps1=psi1-1

C# calculate the probabilities of pj depending on psi1 pj=(p0,p1)

      if (dabs(ps1).gt.1.e-10) then

      d=dsqrt(1+ps1*(psi1*(th1-th2)**2-(th1+th2)**2+2*(th1+th2)))

        p(1)=(d-1+ps1*(th1-th2))/(2*ps1*(1-th2))

        p(2)=p(1)+(1-d+ps1*(th1+th2-2*th1*th2))/(2*ps1*th2*(1-th2))

      else
      p(1)=th1
      p(2)=th1

      end if
      return
      end
