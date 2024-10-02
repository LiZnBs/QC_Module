program decay
    use parameters
    use functions
    implicit none
    call init
    do i=1,nk
        do j=1, nz
            quark_condensate=quark_condensate+trp(kp(i,1), p1, zp(j,1),&
            dse_a(i,nz+1-j),dse_a(i,j),dse_b(i,nz+1-j),dse_b(i,j),&
            f(i,j,1,1), f(i,j,2,1), f(i,j,3,1), f(i,j,4,1))*&
            (2*pi)**(-3)*kp(i,2)*kp(i,1)*sqrt(1-zp(j,1)**2)*zp(j,2)
        end do
    end do
    quark_condensate=quark_condensate*z4*3
    print*, 'The rp is ', quark_condensate
end program decay