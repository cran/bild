C CCCC
C  function to calculate the integral for beta for the first order 
C dependence model with random effects
CCCCC

      DOUBLE PRECISION FUNCTION f1b(v,i)
      DOUBLE PRECISION logL,beta1,bt1,lpsi1,x1,
     *theta1,work1,omega1,v,prob,d0,gbt,glps1,
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
      d0=gbt(i)
      f1b= dexp(logL-(v**2)/(2*dexp(omega1)))*d0
      return
      end
