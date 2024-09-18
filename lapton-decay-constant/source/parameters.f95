module parameters
    implicit none
    ! integrate limit
    real(8), parameter :: o1=-1, o2=1, lambda=100**2
    ! necessary constant
    real*8, parameter :: pi=3.14159265358979
    ! bse mass
    real(8), parameter :: m=0.134
    ! differential need
    ! real(8), parameter :: p1=-m**2, p2=-(m*1.0001)**2, p3=-m**2
    real(8), parameter :: p1=-m**2, p2=-0.1341**2, p3=-m**2
    ! number of gauss points
    integer, parameter :: nk=40, nz=24
    ! DSE and BSE results
    complex(8) ::dse_a1(nk,nz), dse_b1(nk,nz), dse_a2(nk,nz), dse_b2(nk,nz)
    complex(8) :: f(nk,nz,4,1)
    real(8) :: z2
    ! pai
    complex(8) :: pai1=0, pai2=0
    ! renomalization constant
    complex(8) :: c=0
    ! decay constant
    complex(8) :: decay_constant=(0,0)
    ! gauss get point
    real(8) :: kp(nk,2), zp(nz,2)
    ! sum variables
    integer :: i, j, k, l
    !> For imput file， output file and private variables
    integer, parameter :: dse_a1unit=10, dse_b1unit=11, funit=12, z2unit=20, dse_a2unit=13, dse_b2unit=14
end module parameters