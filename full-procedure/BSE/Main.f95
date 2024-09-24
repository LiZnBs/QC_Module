module parameters
    implicit none
    !> Calculate paremeters
    integer :: nk, nz, ny
    real*8 :: pi=3.14159265358979, erf
    real*4 :: lambda
    real*8 :: o1=-1.0 ,o2=1.0
    !> Paramters in QCM
    real(4) :: d, w, mt, tau, lambda_qcd, gamma_m
    !> Data from input
    real(4) :: m
    !> Data from DSE
    real(8) :: z2
    real(8), dimension(:,:), allocatable ::  kp, yp, zp
    complex(8),  dimension(:,:), allocatable ::dse_a, dse_b
    !> Basic calculation parameters in BSE
    complex(8), dimension(:,:,:,:), allocatable :: f, f0
    complex(8) :: eigen, eigen0
    !> For imput file， output file and private variables
    integer :: i, j, k
    integer, parameter :: z2_unit=10, dse_aunit=20, dse_bunit=21, input_file=50
    integer, parameter :: output_file=40
    private :: i, j, k, z2_unit, dse_aunit, dse_bunit, input_file
contains
    !> Init the parameters
    subroutine init
    print*, 'Start calculating BSE'
        !> Get Z2 from file
    open(z2_unit, file='./DSE_Results/Z2.txt', status='old', action='read')
    read(z2_unit, *) z2
    close(z2_unit)
        !> Get d,w,mt,tau,Nf,lambda_qcd,gamma_m from file
    open(input_file, file='../parameters/BSE-parameters.txt', status='old', action='read')
    read(input_file, *) d, w, mt, tau, lambda_qcd, gamma_m,&
    m, lambda, erf, nk, nz, ny
    close(input_file)
    ! Allocate memory 
    allocate(f(nk,nz,4,1), f0(nk,nz,4,1))
    allocate(kp(nk,2), yp(ny,2), zp(nz,2))
    allocate(dse_a(nk,nz), dse_b(nk,nz))
    f=(1,1) 
    f0=(1,1)
    eigen=(1,1)
    eigen0=(0,0)
        !> Get kp,yp,zp 
    call gauleg(o1, o2, kp(:,1), kp(:,2), nk)
    call gauleg(o1, o2, zp(:,1), zp(:,2), nz)
    call gauleg(o1, o2, yp(:,1), yp(:,2), ny)
    kp(:,1) = lambda**kp(:,1)
        !> Get dse_a,dse_b from file
    open(dse_aunit, file='./DSE_Results/Complex-dse_A.txt', status='old', action='read')
    open(dse_bunit, file='./DSE_Results/Complex-dse_B.txt', status='old', action='read')
    do j=1, nz
        do i=1, nk
            read(dse_aunit, *) dse_a(i,j)
            read(dse_bunit, *) dse_b(i,j)
        end do
    end do
    close(dse_aunit)
    close(dse_bunit)
    end subroutine init
end module parameters

program bse
        !> Use module parameters
    use parameters
    implicit none
        !> Interace for functions
    interface
        function k_matrix(k, q, yq, zk, zq, p, ams, apu, bms, bpu) result(kernel)
        implicit none
        !> p={0,0,0,p=IM};
        !> k=k{0,0,Sqrt[1 - zk^2], zk};
        !> q=q{0, Sqrt[1 - zq^2] Sqrt[1 - yq^2], yq Sqrt[1 - zq^2], zq}
        ! k,q are the square of the momentum
        real(8), intent(in) :: k, q, yq, zk, zq
        real(4), intent(in) :: p
        complex(8), intent(in) :: ams, apu, bms, bpu
        complex(8), dimension(4,4) :: kernel
        end function k_matrix       
    end interface

    integer :: qn, qzn, qyn, pn, pzn
    call init
    do while (abs(eigen-eigen0)>erf)
        eigen0=eigen
        f0=f
    !$omp parallel do
    do pn=1, nk
        do pzn=1, nz
        f(pn,pzn,:,1)=(0,0)
            do qn = 1, nk
                do qzn = 1, nz
                    do qyn = 1, ny
                        f(pn,pzn,:,1)=f(pn,pzn,:,1)+(2*pi)**(-3)*kp(qn,2)*kp(qn,1)*log(lambda)&
                        *yp(qyn,2)*zp(qzn,2)/2*kp(qn,1)*sqrt(1-zp(qzn,1)**2)&
                        *matmul(k_matrix(kp(pn,1),kp(qn,1),&
                        yp(qyn,1),&
                        zp(pzn,1),zp(qzn,1),&
                        m,&
                        dse_a(qn,nz+1-qzn),dse_a(qn,qzn),&
                        dse_b(qn,nz+1-qzn),dse_b(qn,qzn)),f0(qn,qzn,:,1))
                    end do
                end do
            end do
        end do
    end do
    eigen=f(1,1,1,1)/f0(1,1,1,1)
    print*, eigen 
    end do
    
    open(output_file, file='./Complex-f.txt', status='replace', action='write')
    write(output_file,*) f
    close(output_file)
end program bse

!> Function to create matrix K
function k_matrix(k, q, yq, zk, zq, p, ams, apu, bms, bpu) result(kernel)
    use parameters
    implicit none
    real(8), intent(in) :: k, q, yq, zk, zq
    real(4), intent(in) :: p
    complex(8), intent(in) :: ams, apu, bms, bpu
    real(8) :: pairkk, pairkq, pairqq, pairpp, k_q
    complex(8) :: pairkp, pairpq
    complex(8) :: glpp
    complex(8), dimension(4,4) :: kernel
    !> Initial variables
    pairkk = k
    pairkq = sqrt(k*q)*(zk*zq + yq*Sqrt(1 - zk**2)*Sqrt(1 - zq**2))
    pairqq = q
    pairpp = -p**2
    pairkp = (0,1)*sqrt(k)*p*zk
    pairpq = (0,1)*p*sqrt(q)*zq
    k_q = pairkk + pairqq - 2*pairkq
    ! Z2^4
    glpp=z2**2*8*pi**2*(d/w**4*exp(-k_q/w**2)+gamma_m/log(tau+(k_q/lambda_qcd**2+1)**2)*(1-exp(-k_q/(4*mt**2)))/k_q)
    !> Initialize kernel
    kernel(1,1)=(16*glpp*(4*bms*bpu - ams*apu*pairpp + 4*ams*apu*pairqq))/&
        ((ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(1,2)= (32*glpp*((apu*bms + ams*bpu)*pairpp + 2*(apu*bms - ams*bpu)*pairpq))/&
        ((ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(1,3)= (32*glpp*pairpq*((apu*bms + ams*bpu)*pairpq + 2*(apu*bms - ams*bpu)*pairqq))/&
        ((ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(1,4)=(64*ams*apu*glpp*(pairpq**2 - pairpp*pairqq))/&
        ((ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(2,1)=(32*glpp*(-(pairkk**2*((apu*bms + ams*bpu)*pairpp + 2*(apu*bms - ams*bpu)*pairpq)) + &
        pairkk*((apu*bms + ams*bpu)*pairkP**2 - 2*apu*bms*pairpq**2 - 2*ams*bpu*pairpq**2 + &
        2*pairkq*((apu*bms + ams*bpu)*pairpp + 4*(apu*bms - ams*bpu)*pairpq) + &
        2*pairkP*((apu*bms - ams*bpu)*pairkq + (apu*bms + ams*bpu)*pairpq) - &
        apu*bms*pairpp*pairqq - ams*bpu*pairpp*pairqq - 6*apu*bms*pairpq*pairqq + &
        6*ams*bpu*pairpq*pairqq) + &
        pairkP*(-((apu*bms + ams*bpu)*pairkP*(4*pairkq - pairqq)) + &
        2*pairkq*((-4*apu*bms + 4*ams*bpu)*pairkq + (apu*bms + ams*bpu)*pairpq + &
        3*(apu*bms - ams*bpu)*pairqq))))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(2,2)=(16*glpp*(pairkk**2*(-(ams*apu*pairpp**2) + 8*ams*apu*pairpq**2 + &
        4*pairpp*(bms*bpu - ams*apu*pairqq)) + &
        pairkP*(2*pairkq*pairpq*(-4*bms*bpu + 16*ams*apu*pairkq + ams*apu*pairpp - &
        8*ams*apu*pairqq) - pairkP*(4*pairkq - pairqq)*&
        (-4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq)) + &
        pairkk*(8*bms*bpu*pairpq**2 - 2*ams*apu*pairpp*pairpq**2 + 4*bms*bpu*pairpp*pairqq - &
        ams*apu*pairpp**2*pairqq + 16*ams*apu*pairpq**2*pairqq - &
        4*ams*apu*pairpp*pairqq**2 + &
        pairkP**2*(-4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq) + &
        2*pairkP*pairpq*(-4*bms*bpu - 4*ams*apu*pairkq + ams*apu*pairpp + &
        4*ams*apu*pairqq) + pairkq*&
        (2*ams*apu*pairpp**2 - 32*ams*apu*pairpq**2 + &
        pairpp*(-8*bms*bpu + 8*ams*apu*pairqq)))))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(2,3)=(16*glpp*pairpq*(pairkk**2*pairpq*(4*bms*bpu - ams*apu*pairpp + 4*ams*apu*pairqq) + &
        pairkP*(2*ams*apu*pairkP*pairpq*pairqq + &
        4*pairkq**2*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq) - &
        pairkq*(8*ams*apu*pairkP*pairpq - 4*ams*apu*pairpq**2 + &
        3*pairqq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq))) + &
        pairkk*(2*ams*apu*pairkP**2*pairpq - &
        pairkP*(-4*ams*apu*pairpq**2 + &
        pairkq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq)) + &
        pairpq*(-4*ams*apu*pairpq**2 - 16*pairkq*(bms*bpu + ams*apu*pairqq) + &
        pairqq*(12*bms*bpu + ams*apu*pairpp + 12*ams*apu*pairqq)))))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(2,4)=(32*glpp*(2*(apu*bms + ams*bpu)*pairkk**2*(pairpq**2 - pairpp*pairqq) + &
        pairkP*(4*pairkq**2*((apu*bms - ams*bpu)*pairpp + 2*(apu*bms + ams*bpu)*pairpq) + &
        pairkP*pairqq*((apu*bms - ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq) - &
        pairkq*((-2*apu*bms + 2*ams*bpu)*pairpq**2 + 3*(apu*bms - ams*bpu)*pairpp*pairqq + &
        2*(apu*bms + ams*bpu)*pairpq*pairqq + &
        4*pairkP*((apu*bms - ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq))) + &
        pairkk*(pairkP**2*((apu*bms - ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq) - &
        2*((apu*bms - ams*bpu)*pairpq - (apu*bms + ams*bpu)*pairqq)*&
        (pairpq**2 - pairpp*pairqq) + &
        pairkP*(-(pairkq*((apu*bms - ams*bpu)*pairpp + 2*(apu*bms + ams*bpu)*pairpq)) + &
        2*pairpq*((apu*bms - ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq)) + &
        2*pairkq*(-4*(apu*bms + ams*bpu)*pairpq**2 + &
        pairpp*((-(apu*bms) + ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq)))))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(3,1)=(64*glpp*((apu*bms + ams*bpu)*pairkP**3 + &
        2*pairkP**2*((apu*bms - ams*bpu)*pairkq - (apu*bms + ams*bpu)*pairpq + &
        (-(apu*bms) + ams*bpu)*pairqq) + &
        pairkP*(-(pairkk*((apu*bms + ams*bpu)*pairpp + (-(apu*bms) + ams*bpu)*pairpq)) + &
        pairkq*((apu*bms + ams*bpu)*pairpp + 4*(-(apu*bms) + ams*bpu)*pairpq) + &
        pairpq*((apu*bms + ams*bpu)*pairpq + 3*(apu*bms - ams*bpu)*pairqq)) + &
        pairpp*(pairkk*((-3*apu*bms + 3*ams*bpu)*pairkq + (apu*bms + ams*bpu)*pairpq + &
        2*(apu*bms - ams*bpu)*pairqq) + &
        pairkq*(4*(apu*bms - ams*bpu)*pairkq - (apu*bms + ams*bpu)*pairpq + &
        3*(-(apu*bms) + ams*bpu)*pairqq))))/&
        (3.*pairkP*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(3,2)=(32*glpp*(2*pairkP**2*(4*bms*bpu - 4*ams*apu*pairkq - ams*apu*pairpp)*pairpq + &
        pairkP**3*(-4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq) + &
        pairpp*pairpq*(pairkk*(-4*bms*bpu + 12*ams*apu*pairkq + ams*apu*pairpp - &
        4*ams*apu*pairqq) + pairkq*&
        (4*bms*bpu - 16*ams*apu*pairkq - ams*apu*pairpp + 8*ams*apu*pairqq)) + &
        pairkP*(pairpq**2*(-4*bms*bpu + ams*apu*pairpp - 8*ams*apu*pairqq) + &
        pairkk*(-(ams*apu*pairpp**2) - 4*ams*apu*pairpq**2 + &
        4*pairpp*(bms*bpu - ams*apu*pairqq)) + &
        pairkq*(ams*apu*pairpp**2 + 16*ams*apu*pairpq**2 + &
        pairpp*(-4*bms*bpu + 4*ams*apu*pairqq)))))/&
        (3.*pairkP*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(3,3)=(16*glpp*pairpq*(4*ams*apu*pairkP**3*pairpq - &
        2*pairkP**2*(4*ams*apu*pairpq**2 + &
        pairkq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq) - &
        pairqq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq)) - &
        pairkP*pairpq*(-4*ams*apu*pairpq**2 + 12*bms*bpu*pairqq + 3*ams*apu*pairpp*pairqq + &
        12*ams*apu*pairqq**2 - 8*pairkq*(2*bms*bpu + ams*apu*pairpp + 2*ams*apu*pairqq) + &
        pairkk*(4*bms*bpu + 5*ams*apu*pairpp + 4*ams*apu*pairqq)) + &
        pairpp*(pairkk*(4*ams*apu*pairpq**2 + &
        3*pairkq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq) - &
        2*pairqq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq)) + &
        pairkq*(-4*ams*apu*pairpq**2 - &
        4*pairkq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq) + &
        3*pairqq*(4*bms*bpu + ams*apu*pairpp + 4*ams*apu*pairqq)))))/&
        (3.*pairkP*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(3,4)=(32*glpp*(2*pairkP**3*((apu*bms - ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq) - &
        2*pairkP**2*(2*(apu*bms - ams*bpu)*pairpq**2 + &
        pairkq*((apu*bms - ams*bpu)*pairpp + 2*(apu*bms + ams*bpu)*pairpq) + &
        (-(apu*bms) + ams*bpu)*pairpp*pairqq + 2*(apu*bms + ams*bpu)*pairpq*pairqq) - &
        pairkP*(pairpq*((-2*apu*bms + 2*ams*bpu)*pairpq**2 + &
        3*(apu*bms - ams*bpu)*pairpp*pairqq + 2*(apu*bms + ams*bpu)*pairpq*pairqq) - &
        2*pairkq*(4*(apu*bms + ams*bpu)*pairpq**2 + &
        pairpp*(3*(apu*bms - ams*bpu)*pairpq + 2*(apu*bms + ams*bpu)*pairqq)) + &
        pairkk*(2*(apu*bms + ams*bpu)*pairpq**2 + &
        pairpp*(3*(apu*bms - ams*bpu)*pairpq + 4*(apu*bms + ams*bpu)*pairqq))) + &
        pairpp*(pairkq*((-2*apu*bms + 2*ams*bpu)*pairpq**2 - &
        4*pairkq*((apu*bms - ams*bpu)*pairpp + 2*(apu*bms + ams*bpu)*pairpq) + &
        3*(apu*bms - ams*bpu)*pairpp*pairqq + 2*(apu*bms + ams*bpu)*pairpq*pairqq) + &
        pairkk*(3*pairkq*((apu*bms - ams*bpu)*pairpp + 2*(apu*bms + ams*bpu)*pairpq) + &
        2*(apu*bms - ams*bpu)*(pairpq**2 - pairpp*pairqq)))))/&
        (3.*pairkP*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(4,1)=(64*ams*apu*glpp*(2*pairkP**2*pairqq + pairkq*pairpp*pairqq - &
        pairkP*pairpq*(2*pairkq + pairqq) + &
        pairkk*(pairkq*pairpp - pairkP*pairpq + 2*pairpq**2 - 2*pairpp*pairqq)))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(4,2)=(-64*(apu*bms + ams*bpu)*glpp*(2*pairkP**2*pairqq + pairkq*pairpp*pairqq - &
        pairkP*pairpq*(2*pairkq + pairqq) + &
        pairkk*(pairkq*pairpp - pairkP*pairpq + 2*pairpq**2 - 2*pairpp*pairqq)))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(4,3)=(32*(apu*bms - ams*bpu)*glpp*pairpq*&
        (2*pairkP**2*pairqq + pairkq*pairpp*pairqq - pairkP*pairpq*(2*pairkq + pairqq) + &
            pairkk*(pairkq*pairpp - pairkP*pairpq + 2*pairpq**2 - 2*pairpp*pairqq)))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
    kernel(4,4)=(16*glpp*(4*bms*bpu + ams*apu*pairpp - 4*ams*apu*pairqq)*&
        (-2*pairkP**2*pairqq - pairkq*pairpp*pairqq + pairkP*pairpq*(2*pairkq + pairqq) + &
        pairkk*(-(pairkq*pairpp) + pairkP*pairpq - 2*pairpq**2 + 2*pairpp*pairqq)))/&
        (3.*(pairkP**2 - pairkk*pairpp)*(pairkk - 2*pairkq + pairqq)*&
        (ams**2*pairpp + 4*(bms**2 - ams**2*pairpq + ams**2*pairqq))*&
        (apu**2*pairpp + 4*(bpu**2 + apu**2*pairpq + apu**2*pairqq)))
end function k_matrix

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
