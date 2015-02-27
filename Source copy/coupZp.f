
c input: vector and axial Zp couplings to up and down quarks
c output: left and right chiral couplings to up and down quarks 

      subroutine coupZp(model)
      implicit real*8 (a-h,o-z)
      common/coupZp/gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
      common/ZpLRcoup/gZpd(2,5),gZpu(2,5)

c local variables
      real*8 fac(5)   

      do i=1,5
        fac(i) = gp(i)/2.d0  
        gZpd(1,i) = fac(i)*(gV_d(i)+gA_d(i))
        gZpd(2,i) = fac(i)*(gV_d(i)-gA_d(i))
        gZpu(1,i) = fac(i)*(gV_u(i)+gA_u(i))
        gZpu(2,i) = fac(i)*(gV_u(i)-gA_u(i))
      enddo

!       print*, gp
!       print*, gV_d ,gA_d
!       print*, gV_u ,gA_u
!       print*, gZpd(1) ,gZpd(2)
!       print*, gZpu(1) ,gZpu(2)
 
      return
      end