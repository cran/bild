C CCCC
C  function to calculate the integral of the log-likelihood
C  for the first order dependence model with random effects
CCCCC

      DOUBLE PRECISION FUNCTION fpij(v,i)
      DOUBLE PRECISION logL,v,prob,
     *bt1,beta1, lpsi1,x1,theta1,work1,omega1
      INTEGER m,mpar,y1,i,k

      DIMENSION x1(5000,10),theta1(5000),
     *work1(5000),y1(5000),prob(5000),beta1(10),bt1(10)

      COMMON/param/x1,theta1,work1,
     *y1,beta1,bt1,m,mpar,omega1,lpsi1
      beta1(1) = v +bt1(1)

      CALL mlik1i(logL,prob,mpar,m)
      fpij= prob(i)*dexp(-v**2/(2.0d0*dexp(omega1)))
      return
      end
