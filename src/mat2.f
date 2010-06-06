      subroutine mat2(p1,p2,P)
C       creates 2x2 MC transition matrix, given probs of 1's
      implicit double precision (a-h,o-z)
      dimension P(2,2)
      P(1,1)=1-p1
      P(1,2)=p1
      P(2,1)=1-p2
      P(2,2)=p2
      return
      end

 
