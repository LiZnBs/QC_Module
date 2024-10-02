module parameters
    implicit none
    ! integrate limit
    real(8), parameter :: o1=-1, o2=1
    ! necessary constant
    real*8, parameter :: pi=3.14159265358979
    ! bse mass and differential need
    real(8) :: m, lambda, p1, p2, p3
    ! number of gauss points
    integer :: nk, nz
    ! DSE and BSE results
    complex(8),  dimension(:,:), allocatable :: dse_a1, dse_b1, dse_a2, dse_b2
    complex(8), dimension(:,:,:,:), allocatable :: f
    real(8) :: z2
    ! pai
    complex(8) :: pai1=0, pai2=0
    ! renomalization constant
    complex(8) :: c=0
    ! decay constant
    complex(8) :: decay_constant=(0,0)
    ! gauss get point
    real(8), dimension(:,:), allocatable ::  kp, zp
    ! sum variables
    integer :: i, j, k, l
    !> For imput file， output file and private variables
    integer, parameter :: dse_a1unit=10, dse_b1unit=11, funit=12, &
    z2unit=20, dse_a2unit=13, dse_b2unit=14, input_file=30, export_file=40
end module parameters