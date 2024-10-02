program decay
    use parameters
    use functions
    implicit none
    call init
    do i=1,nk
        do j=1, nz
            pai1=pai1+trpi(kp(i,1),p1, p3, zp(j,1),&
            dse_a1(i,nz+1-j), dse_a1(i,j), dse_b1(i,nz+1-j), dse_b1(i,j),&
            f(i,j,1,1), f(i,j,2,1), f(i,j,3,1), f(i,j,4,1))&
            *(2*pi)**(-3)*kp(i,2)*kp(i,1)*sqrt(1-zp(j,1)**2)*zp(j,2)
            pai2=pai2+trpi(kp(i,1),p2, p3, zp(j,1),&
            dse_a2(i,nz+1-j), dse_a2(i,j), dse_b2(i,nz+1-j), dse_b2(i,j),&
            f(i,j,1,1), f(i,j,2,1), f(i,j,3,1), f(i,j,4,1))&
            *(2*pi)**(-3)*kp(i,2)*kp(i,1)*sqrt(1-zp(j,1)**2)*zp(j,2)
        end do
    end do
    c=((pai1-pai2)/(p1-p2))**(-0.5)
    f=c*f
    do i=1,nk
        do j=1, nz
            decay_constant=decay_constant+trf(kp(i,1), p3, zp(j,1),&
            dse_a1(i,nz+1-j),dse_a1(i,j),dse_b1(i,nz+1-j),dse_b1(i,j),&
            f(i,j,1,1), f(i,j,2,1), f(i,j,3,1), f(i,j,4,1))*&
            (2*pi)**(-3)*kp(i,2)*kp(i,1)*sqrt(1-zp(j,1)**2)*zp(j,2)
        end do
    end do
    decay_constant=decay_constant*z2*3/p3
    print*, 'The decay constant of the pion is ', decay_constant
    call export
end program decay