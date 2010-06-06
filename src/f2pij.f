C CCCC
C  function to calculate the integral of the log-likelihood
C  for the first order dependence model with random effects
CCCCC

      DOUBLE PRECISION FUNCTION f2pij(v,i)
      DOUBLE PRECISION logL,v,prob,
     *bt1,beta1, lpsi1,x1,theta1,work1,omega1
      INTEGER m,mpar,y1,i


      DIMENSION x1(4500,10),theta1(4500),
     *work1(4500),y1(4500),prob(4500),lpsi1(2),
     *beta1(10),bt1(10)

      COMMON/param/x1,theta1,work1,
     *y1,lpsi1,beta1,bt1,m,mpar,omega1,logL
      beta1(1) = v +bt1(1)

      CALL mlik2i(logL,prob,mpar,m)
      f2pij= prob(i)*dexp(-v**2/(2*dexp(omega1)))
      return
      end
