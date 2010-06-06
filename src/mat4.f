      subroutine mat4(p1,p2,p3,p4,P)
C   creates 4x4 MC transition matrix, given probs of 1's
      implicit double precision (a-h,o-z)
      dimension P(4,4)
      P(1,1)=1-p1
      P(1,2)=p1
      P(1,3)=0
      P(1,4)=0
      P(2,1)=0
      P(2,2)=0
      P(2,3)=1-p2
      P(2,4)=p2
      P(3,1)=1-p3
      P(3,2)=p3
      P(3,3)=0
      P(3,4)=0
      P(4,1)=0
      P(4,2)=0
      P(4,3)=1-p4
      P(4,4)=p4
      return
      end
