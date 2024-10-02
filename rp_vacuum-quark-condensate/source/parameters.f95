module parameters
    implicit none
    ! integrate limit
    real(8), parameter :: o1=-1, o2=1, lambda=100**2
    ! necessary constant
    real*8, parameter :: pi=3.14159265358979
    ! bse mass
    real(8), parameter :: m=0.134
    ! differential need
    real(8), parameter :: p1=-m**2
    ! number of gauss points
    integer, parameter :: nk=40, nz=24
    ! DSE and BSE results
    complex(8) ::dse_a(nk,nz), dse_b(nk,nz)
    complex(8) :: f(nk,nz,4,1)
    real(8) :: z4
    ! decay constant
    complex(8) :: quark_condensate=(0,0)
    ! gauss get point
    real(8) :: kp(nk,2), zp(nz,2)
    ! sum variables
    integer :: i, j, k, l
    !> For imput file， output file and private variables
    integer, parameter :: dse_aunit=10, dse_bunit=11, funit=12, z4unit=20
end module parameters