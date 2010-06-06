      subroutine matc(a,b,n1,n2)
C   b=a
      implicit double precision (a-h,o-z)
      dimension a(n1,n2),b(n1,n2)
      do 10 i=1,n1
        do 20 j=1,n2
          b(i,j)=a(i,j)
   20   continue
   10 continue 
      return
      end
