
      subroutine mati(a,b,c,n1,n2,n3,n,npar)
C      ! computes c=a*b
      implicit double precision (a-h,o-z)
      dimension a(n1,n2),b(n2,n3),c(n1,n3)    
      do 10 i=1,n
        do 20 j=1,n3
          s=0
            do 30 k=1,(npar-2)
              s=s+a(i,k)*b(k,j)
   30        continue
          c(i,j)=s
   20   continue
   10 continue 
      return
      end
