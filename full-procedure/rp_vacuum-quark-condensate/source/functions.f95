module functions
use parameters
    implicit none
    private
    public :: gauleg
    public :: init
    public :: trp
    real(8), parameter :: EPS = 3.d-16
contains

! variable k, p, pz are all momentum square
! P={0,0,0,i*sqrt(-p)}
! k={0,0,sqrt(1-kz^2),kz}*sqrt(k)
! Z={0,0,0,i*sqrt(-pz)}
! p,pz<0
function trp(k, p, kz, ams, apu, bms, bpu, f1, f2, f3, f4) result(res)
    implicit none
    real(8), intent(in) :: k, p, kz
    complex(8), intent(in) :: ams, apu, bms, bpu, f1, f2, f3, f4
    real(8) :: pairpp, pairkk
    complex(8) :: pairkp
    complex(8) :: res
    ! initial variables
    pairpp=p
    pairkp=(0d0,1d0)*kz*sqrt(k)*sqrt(-p)
    pairkk=k
    ! calculate the result
    res=(16*(4*bms*bpu*f1 + 4*(apu*bms - ams*bpu)*f2*pairkp + &
    2*(apu*bms*f3 + ams*bpu*f3 + 2*ams*apu*f4)*pairkp**2 - &
    ams*apu*f1*pairpp + 2*apu*bms*f2*pairpp + 2*ams*bpu*f2*pairpp + &
    4*pairkk*((apu*bms - ams*bpu)*f3*pairkp + ams*apu*(f1 - f4*pairpp))))&
    /((4*bms**2 + 4*ams**2*pairkk - 4*ams**2*pairkp + ams**2*pairpp)*&
    (4*bpu**2 + 4*apu**2*pairkk + 4*apu**2*pairkp + apu**2*pairpp))
end function trp

subroutine init
        print*, 'Start calculating the rp (vacuum-quark-condensate) !'
        ! get parameters from file
        open(input_file, file='../parameters/rp-parameters.txt', status='old', action='read')
        read(input_file, *) m, lambda, nk, nz
        close(input_file)
        p1=-m**2
        ! Allocate memory 
        allocate(f(nk,nz,4,1))
        allocate(kp(nk,2), zp(nz,2))
        allocate(dse_a(nk,nz), dse_b(nk,nz))
        f=(1,1) 
        ! get points
        call gauleg(o1, o2, kp(:,1), kp(:,2), nk)
        call gauleg(o1, o2, zp(:,1), zp(:,2), nz)
        kp(:,1) = lambda**kp(:,1)
        kp(:,2) = kp(:,2)*kp(:,1)*log(lambda)
        !> Get dse_a1,dse_b1 from file
    open(dse_aunit, file='./DSE-BSE_results/Complex-dse_A.txt', status='old', action='read')
    open(dse_bunit, file='./DSE-BSE_results/Complex-dse_B.txt', status='old', action='read')
    do j=1, nz
        do i=1, nk
            read(dse_aunit, *) dse_a(i,j)
            read(dse_bunit, *) dse_b(i,j)
        end do
    end do
    close(dse_aunit)
    close(dse_bunit)
    !> Get f from file
    open(funit, file='./DSE-BSE_results/normalized-f.txt', status='old', action='read')
    read (funit, *) f
    close(funit)
    ! get z4 from file
    open(z4unit, file='./DSE-BSE_results/Z4.txt', status='old', action='read')
    read(z4unit, *) z4
    close(z4unit)
end subroutine init

    ! gauss get points
SUBROUTINE gauleg(x1,x2,x,w,n)
    INTEGER n
    REAL(8) x1,x2,x(n),w(n)
    DOUBLE PRECISION EPS
    PARAMETER (EPS=3.d-16)
    INTEGER i,j,mm
    DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
    mm=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do 12 i=1,mm
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