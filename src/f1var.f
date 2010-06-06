C CCCC
C  function to calculate the integral for omega for the first order 
C dependence model with random effects
C dimension of x and vectors assume the maximum values
CCCCC

      DOUBLE PRECISION FUNCTION f1var(v,i)
      DOUBLE PRECISION logL,prob,beta1,bt1,lpsi1,x1,
     *theta1,work1,omega1,v
      INTEGER m, mpar, y1,i
      DIMENSION x1(5000,10),theta1(5000),
     *work1(5000),y1(5000),prob(5000),beta1(10),bt1(10)

      COMMON/param/x1,theta1,work1,
     *y1,beta1,bt1,m,mpar,omega1,lpsi1

      beta1(1) = v +bt1(1)

      CALL mlik1i(logL,prob,mpar,m)

      f1var = dexp(logL-(v**2)/(2*dexp(omega1)))*
     *((v**2-dexp(omega1))/(2*dexp(2*omega1)))
      return
      end
