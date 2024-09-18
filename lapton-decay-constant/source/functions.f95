module functions
use parameters
    implicit none
    private
    public :: gauleg
    public :: init
    public :: trpi
    public :: trf
    real(8), parameter :: EPS = 3.d-16
contains
! variable k, p are all momentum square
! P={0,0,0,i*sqrt(-p)}
! k={0,0,sqrt(1-kz^2),kz}*sqrt(k)
! p<0
function trf(k, p, kz, ams, apu, bms, bpu, f1, f2, f3, f4) result(res)
    implicit none
    real(8), intent(in) :: k, p, kz
    complex(8), intent(in) :: ams, apu, bms, bpu, f1, f2, f3, f4
    real(8) :: pairpp, pairkk
    complex(8) :: pairkp
    complex(8) :: res
    ! initial variables
    pairpp=p
    pairkp=(0,1)*kz*sqrt(k)*sqrt(-p)
    pairkk=k
    ! calculate the result
    res=(16*(4*(apu*bms - ams*bpu)*f1*pairkp + pairpp*(2*apu*bms*f1&
    + 2*ams*bpu*f1 - 4*bms*bpu*f2 + 4*(ams*apu*f2 + apu*bms*f4 + &
    ams*bpu*f4)*pairkk + ams*apu*f2*pairpp) - pairkp**2*(8*ams*apu*f2&
    + 4*bms*bpu*f3 + 4*apu*bms*f4 + 4*ams*bpu*f4 + 4*ams*apu*f3*&
    pairkk - ams*apu*f3*pairpp)))/((4*bms**2 + 4*ams**2*pairkk -&
    4*ams**2*pairkp + ams**2*pairpp)*(4*bpu**2 + 4*apu**2*pairkk +&
    4*apu**2*pairkp + apu**2*pairpp))
end function trf

! variable k, p, pz are all momentum square
! P={0,0,0,i*sqrt(-p)}
! k={0,0,sqrt(1-kz^2),kz}*sqrt(k)
! Z={0,0,0,i*sqrt(-pz)}
! p,pz<0
function trpi(k, p, pz, kz, ams, apu, bms, bpu, f1, f2, f3, f4) result(res)
    implicit none
    real(8), intent(in) :: k, p, pz, kz
    complex(8), intent(in) :: ams, apu, bms, bpu, f1, f2, f3, f4
    real(8) :: pairpp, pairpz, pairzz, pairkk
    complex(8) :: pairkz, pairkp
    complex(8) :: res
    ! initial variables
    pairpp=p
    pairkz=(0d0,1d0)*kz*sqrt(k)*sqrt(-pz)
    pairpz=-sqrt(-p)*sqrt(-pz)
    pairzz=pz
    pairkp=(0d0,1d0)*kz*sqrt(k)*sqrt(-p)
    pairkk=k
    ! calculate the result
    res=(96*(-4*bms*bpu*f1**2 + ams*apu*f1**2*pairpp + pairkz**2*&
    (8*ams*apu*f2**2 + 8*bms*bpu*f2*f3 + 8*apu*bms*f2*f4 + 8*ams*&
    bpu*f2*f4 + 4*bms*bpu*f4**2 + 4*(-(apu*bms) + ams*bpu)*f3*f4*&
    pairkp - 2*ams*apu*f3**2*pairkp**2 + ams*apu*(2*f2*f3 - f4**2)&
    *pairpp) - 4*apu*bms*f1*f2*pairpz - 4*ams*bpu*f1*f2*pairpz - 2&
    *ams*apu*f2**2*pairpz**2 - 4*pairkz*((apu*bms - ams*bpu)*f2*(2&
    *f1 - f4*pairpz) + pairkp*(f1*(apu*bms*f3 + ams*bpu*f3 + 2*ams*&
    apu*f4) + ams*apu*(f2*f3 - f4**2)*pairpz)) + 4*bms*bpu*f2**2*&
    pairzz - 4*apu*bms*f2*f4*pairkp*pairzz + 4*ams*bpu*f2*f4*pairkp*&
    pairzz - 2*ams*apu*f4**2*pairkp**2*pairzz + ams*apu*f2**2*pairpp*&
    pairzz + 4*ams*apu*pairkk**2*(f3**2*pairkz**2 + f4**2*pairzz) + &
    pairkk*(-4*ams*apu*f1**2 + pairkz**2*(8*ams*apu*f2*f3 + 4*bms*bpu*&
    f3**2 - 4*ams*apu*f4**2 + ams*apu*f3**2*pairpp) + 8*ams*apu*f1*&
    f4*pairpz - 2*ams*apu*f4**2*pairpz**2 - 4*(apu*bms - ams*bpu)*&
    f3*pairkz*(2*f1 - f4*pairpz) - 4*ams*apu*f2**2*pairzz - 8*apu*&
    bms*f2*f4*pairzz - 8*ams*bpu*f2*f4*pairzz - 4*bms*bpu*f4**2*&
    pairzz + ams*apu*f4**2*pairpp*pairzz)))/((4*bms**2 + 4*ams**2*&
    pairkk - 4*ams**2*pairkp + ams**2*pairpp)*(4*bpu**2 + 4*apu**2*&
    pairkk + 4*apu**2*pairkp + apu**2*pairpp))
end function trpi
subroutine init
        print*, 'Start calculating the decay constant of the pion'
        ! get points
        call gauleg(o1, o2, kp(:,1), kp(:,2), nk)
        call gauleg(o1, o2, zp(:,1), zp(:,2), nz)
        kp(:,1) = lambda**kp(:,1)
        kp(:,2) = kp(:,2)*kp(:,1)*log(lambda)
        !> Get dse_a1,dse_b1 from file
    open(dse_a1unit, file='./DSE-BSE_results/Complex-dse_A1.txt', status='old', action='read')
    open(dse_b1unit, file='./DSE-BSE_results/Complex-dse_B1.txt', status='old', action='read')
    do j=1, nz
        do i=1, nk
            read(dse_a1unit, *) dse_a1(i,j)
            read(dse_b1unit, *) dse_b1(i,j)
        end do
    end do
    close(dse_a1unit)
    close(dse_b1unit)
    ! get dse_a2,dse_b2 from file
    open(dse_a2unit, file='./DSE-BSE_results/Complex-dse_A2.txt', status='old', action='read')
    open(dse_b2unit, file='./DSE-BSE_results/Complex-dse_B2.txt', status='old', action='read')
    do j=1, nz
        do i=1, nk
            read(dse_a2unit, *) dse_a2(i,j)
            read(dse_b2unit, *) dse_b2(i,j)
        end do
    end do
    close(dse_a2unit)
    close(dse_b2unit)
    !> Get f from file
    open(funit, file='./DSE-BSE_results/Complex-f.txt', status='old', action='read')
    read (funit, *) f
    close(funit)
    ! get z2 from file
    open(z2unit, file='./DSE-BSE_results/Z2.txt', status='old', action='read')
    read(z2unit, *) z2
    close(z2unit)
end subroutine init

    ! gauss get points
SUBROUTINE gauleg(x1,x2,x,w,n)
    INTEGER n
    REAL(8) x1,x2,x(n),w(n)
    DOUBLE PRECISION EPS
    PARAMETER (EPS=3.d-16)
    INTEGER i,j,m
    DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do 12 i=1,m
        z=cos(3.14159265358979d0*(i-.25d0)/(n+.5d0))
        1       continue
        p1=1.d0
        p2=0.d0
        do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        11        continue
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
    12    continue
    return
END  

end module functions