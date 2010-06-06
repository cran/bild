C CCCC
C  function to calculate the integral for psi for the first order 
C dependence model with random effects
C dimension of x and vectors assume the maximum values
CCCCC

      DOUBLE PRECISION FUNCTION fpsi(v,i)
      DOUBLE PRECISION logL,prob,beta1,bt1,lpsi1,x1,
     *theta1,work1,omega1,v,d2,gbt,glps1,
     *dbeta,dbeta1,der,db
      INTEGER m, mpar, y1,i
      DIMENSION x1(5000,10),theta1(5000),
     *work1(5000),y1(5000),prob(5000),beta1(10),bt1(10),
     *gbt(10),dbeta(10),dbeta1(10),der(10),db(3,10)

      COMMON/param/x1,theta1,work1,
     *y1,beta1,bt1,m,mpar,omega1,lpsi1
      COMMON/grad/ dbeta,dbeta1,der,db

      beta1(1) = v +bt1(1)

      CALL mlik1i(logL,prob,mpar,m)
      
      CALL mbgd1i(gbt,glps1,mpar,m)
      d2=glps1
      fpsi = dexp(logL-(v**2)/(2*dexp(omega1)))*d2
      return
      end
