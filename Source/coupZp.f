
! input: vector and axial Zp couplings to up and down quarks
! output: left and right chiral couplings to up and down quarks

    subroutine coupZpx!(o_NWA,model_name)
    implicit real (a-h,o-z)
    common/coupZpVA/gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
    common/coupZp/gZpd(2,5),gZpu(2,5)

    do i=1,5
        gZpd(1,i) = gp(i)*(gV_d(i)+gA_d(i))/2.d0
        gZpd(2,i) = gp(i)*(gV_d(i)-gA_d(i))/2.d0
        gZpu(1,i) = gp(i)*(gV_u(i)+gA_u(i))/2.d0
        gZpu(2,i) = gp(i)*(gV_u(i)-gA_u(i))/2.d0
    enddo
     
    return
    end subroutine coupZpx
