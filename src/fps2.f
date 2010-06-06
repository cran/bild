C CCCC
C  function to calculate the integral for psi2 for the second order 
C dependence model with random effects
C dimension of x and vectors assume the maximum values
CCCCC

      DOUBLE PRECISION FUNCTION fps2(v,i)
      DOUBLE PRECISION logL,prob,beta1,bt1,lpsi1,x1,
     *theta1,work1,omega1,v,d3,gbt,glps1,glps2,
     *dbeta,dbeta1,dbeta2,der,db,db1,db2
      INTEGER m, mpar, y1,i
      DIMENSION x1(4500,10),theta1(4500),
     *work1(4500),y1(4500),prob(4500),lpsi1(2),beta1(10),
     *bt1(10),dbeta(10),dbeta1(10),dbeta2(10),
     *der(10),db(3,10),db1(4,10),db2(5,10),gbt(10)

      COMMON/param/x1,theta1,work1,
     *y1,lpsi1,beta1,bt1,m,mpar,omega1
      COMMON/grad/ dbeta,dbeta1,dbeta2,
     *der,db,db1,db2

      beta1(1) = v +bt1(1)

      CALL mlik2i(logL,prob,mpar,m)
  
      CALL mbgd2i(gbt,glps1,glps2,mpar,m)
      d3=glps2
      fps2 = dexp(logL-(v**2)/(2*dexp(omega1)))*d3
      return
      end
