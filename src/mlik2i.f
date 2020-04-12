C # Calculates the log-likelihood for the second
C# order dependence model

      subroutine mlik2i(logL,pij,npar,n)
      implicit double precision (a-h,o-z)
      DIMENSION x1(4500,10),theta1(4500),
     *work1(4500),y1(4500),lpsi1(2),
     *beta1(10),bt1(10),pij(n),psi(2),tpr(2),tpr1(4),
     *P0(2,2), P1(2,2), P2(2,2),P3(4,4), P4(4,4), P5(4,4), 
     *P6(4,4),Pc(4,1),Pr(2,1),Pr0(4,1),Pr1(2,1)

      double precision logL,lpsi1,pij
      integer y1,i0,i1,i2,npar,n,n0,k,m,mpar

      COMMON/param/x1,theta1,work1,
     *y1,lpsi1,beta1,bt1,m,mpar,omega1

      psi(1) = dexp(lpsi1(1))
      psi(2) = dexp(lpsi1(2))
      psi1 = psi(1)
      psi2 = psi(2)
      ps1 = psi1-1
      ps2 = psi2-1

      call mati(x1,beta1,work1,4500,10,1,n,npar)

      do 20 i=1,n
        theta1(i) = 1/(1+dexp(-work1(i)))
   20 continue

      i0 = 1  
   25   if (y1(i0).eq.(-1)) then
      i0=i0+1
      go to 25
      end if

      n0 = n
   30   if (y1(n0).eq.(-1)) then
      n0=n0-1
      go to 30
      end if

      logL = 0
      p  = theta1(i0)
      pij(i0)=p
      logL = y1(i0)*dlog(p/(1-p))+dlog(1-p)
      if (i0.eq.n0) return
      i = i0+1
      i1 = i
   40   if (y1(i1).eq.(-1)) then
      i1=i1+1
      go to 40
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
      if (i1.eq.i) then
      th1 = theta1(i)
      th2 = theta1(i-1)
      call mcpj(th1,th2,psi1,tpr)
      else 
C    (one or more intermediate missing datum)
      call mat2 (0.0D0,1.0D0,P0)
      do 50 k=(i0+1),i1
        th1 = theta1(k)
        th2 = theta1(k-1)
        call mcpj (th1,th2,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)   
        call matp(P0,P1,P2,2,2,2)
        call matc(P2,P0,2,2)
   50 continue
      tpr(1)= P0(1,2)
      tpr(2)= P0(2,2)
      end if
      p=tpr(y1(i0)+1)
      pij(i1)=p
      logL = logL+y1(i1)*dlog(p/(1-p))+dlog(1-p)

      if (i1.eq.n0) return
      i=i1+1
   60 if (i.le.n0) then
      i2 = i
   70   if (y1(i2).eq.(-1)) then
      i2=i2+1
      go to 70
      end if

      if (i2.eq.i.and.i1.eq.(i0+1)) then
        th = theta1(i2)
        th1 = theta1(i1)
        th2 = theta1(i0)
      call mcpij(th,th1,th2,psi1,psi2,tpr1)
        p=tpr1(y1(i0)+2*y1(i1)+1)
        pij(i2)=p
        logL = logL+y1(i2)*dlog(p/(1-p))+dlog(1-p)
        i0=i1
        i1=i2
        i=i1+1

      else if (i1.ne.(i0+1).and.i2.eq.(i1+1)) then
C    (one intermediate missing datum between i0 and i1)
      call mat2 (0.0D0,1.0D0,P0)
      do 75 k=(i0+1),(i1-1)
        th1 = theta1(k)
        th2 = theta1(k-1)
        call mcpj (th1,th2,psi1,tpr)
        call mat2 (tpr(1),tpr(2),P1)   
        call matp(P0,P1,P2,2,2,2)
        call matc(P2,P0,2,2)
   75 continue
        th = theta1(i2)
        th1 = theta1(i1)
        th2 = theta1(i1-1)
        call mcpij (th,th1,th2,psi1,psi2,tpr1)
        Pr(1,1)=tpr1(2*y1(i1)+1)
        Pr(2,1)=tpr1(2*y1(i1)+2)
        call matp(P0,Pr,Pr1,2,2,1)
        p=Pr1(y1(i0)+1,1)
        pij(i2)=p
        logL = logL+y1(i2)*dlog(p/(1-p))+dlog(1-p)
        i0=i1
        i1=i2
        i=i1+1

      else if (i2.ne.(i1+1).and.i2.ne.n0) then
C    (one or more intermediate missing datum between i1 and i2)
      call mat4 (0.0D0,0.0D0,1.0D0,1.0D0,P3)
      call matp (P3,P3,P4,4,4,4)
      do 80 k=(i1+1),i2
        th = theta1(k)
        th1 = theta1(k-1)
        th2 = theta1(k-2)
        call mcpij (th,th1,th2,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P5)
        call matp(P4,P5,P6,4,4,4)
        call matc(P6,P4,4,4)
   80 continue

        i3=i2+1
        th = theta1(i3)
        th1 = theta1(i3-1)
        th2 = theta1(i3-2)
        call mcpij (th,th1,th2,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P5)
        call matp(P4,P5,P6,4,4,4)
        pstar=P6(2*y1(i0)+y1(i1)+1,2*y1(i2)+y1(i3)+1)
        pij(i2)=pstar
        logL = logL+dlog(pstar)
        i0=i2
        i1=i3
        i=i1+1

      else if (i2.ne.(i1+1).and.i2.eq.n0) then
C    (one intermediate missing datum between i1 and i2(last value))
      call mat4 (0.0D0,0.0D0,1.0D0,1.0D0,P3)
      call matp (P3,P3,P4,4,4,4)
      do 90 k=(i1+1),i2
        th = theta1(k)
        th1 = theta1(k-1)
        th2 = theta1(k-2)
        call mcpij (th,th1,th2,psi1,psi2,tpr1)
        call mat4 (tpr1(1),tpr1(3),tpr1(2),tpr1(4),P5)
        call matp(P4,P5,P6,4,4,4)
        call matc(P6,P4,4,4)
   90 continue
        Pc(1,1)=0
        Pc(2,1)=1
        Pc(3,1)=0
        Pc(4,1)=1
        call matp(P4,Pc,Pr0,4,4,1)
        p=Pr0(2*y1(i0)+y1(i1)+1,1)
        pij(i2)=p
        logL = logL+y1(i2)*dlog(p/(1-p))+dlog(1-p)
        i0=i1
        i1=i2
        i=i1+1
       end if
      go to 60
      end if
      return
      end
