module parameters
    implicit none
    ! integrate limit
    real(8), parameter :: o1=-1, o2=1
    ! necessary constant
    real*8, parameter :: pi=3.14159265358979
    ! bse mass and differential need
    real(8) :: m, lambda, p3
    ! number of gauss points
    integer :: nk, nz
    ! DSE and BSE results
    complex(8),  dimension(:,:), allocatable :: dse_a, dse_b
    complex(8), dimension(:,:,:,:), allocatable :: f
    real(8) :: z2
    ! decay constant
    complex(8) :: decay_constant=(0,0)
    ! gauss get point
    real(8), dimension(:,:), allocatable ::  kp, zp
    ! sum variables
    integer :: i, j, k, l
    !> For imput file， output file and private variables
    integer, parameter :: dse_aunit=10, dse_bunit=11, funit=12, &
    z2unit=20, input_file=30
end module parameters