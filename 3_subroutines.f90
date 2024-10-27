!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++       S U B R O U T I N E S        ++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++    Calculo de los Coeficientes de Fourier    +++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Parte que implementa el método numérico.
! mp : m_max + 1
!
subroutine get_coef(n,mp,pa,h,p_th,coef)
    use constants
    implicit none
    integer(dp) :: i,j
    integer(dp),intent(in) :: n,mp,pa
    real(dp),intent(inout) :: h, p_th
    real(dp) :: radio,error,w
    complex(dp),dimension(mp,n),intent(inout) :: coef
    real(dp),dimension(:,:),allocatable :: A
    real(dp),dimension(:),allocatable :: sol_aux1,B
    

    radio = 1.0_dp
    error = (10.0_dp**(-6.0_dp))
    w = 1.5_dp
    p_th = 2.0_dp*pi/real(pa-1, dp)
    h = radio/real(n)

    allocate(A(n*2,n*2))
    allocate(B(n*2))
    allocate(sol_aux1(n*2))

    
    call m_zeros(A,n*2,n*2)
    call v_zeros(B,n*2)
    call m_zeros_complex(coef,mp,n)
    
    do i=0,mp-1,1
        call matrix_problem(A,B,n,i,h)
        call sor(A,B,n*2,error,w,.true.,sol_aux1)
        do j=1,n,1
            coef(i+1,j)=complex(sol_aux1(j),sol_aux1(j+n))
        end do
    end do

    deallocate(A)
    deallocate(B)
    deallocate(sol_aux1)

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++    Vector de puntos igualmente espaciados    +++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Se incluyen los extremos 'from' y 'to'.
!
subroutine linspace(from,to,vec,n)
    use constants
    implicit none
    integer(dp) :: i
    integer(dp), intent(in) :: n
    real(dp), intent(in) :: from, to
    real(dp), dimension(n), intent(inout) :: vec
    real(dp) :: step

    select case (n)
        case (0)
            return
        case (1)
            vec(1) = from
            return
    end select

    step = (to - from) / real(n-1, dp)

    do i=1,n,1
        vec(i) = from + (step*real(i-1, dp))
    end do

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++    Completa el Sistema a Resolver    +++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Depende del modo (m) que corresponde resolver
!
subroutine matrix_problem(A,B,n,m,h)
    use constants
    implicit none
    integer(dp) :: i
    integer(dp),intent(in) :: n,m
    real(dp),intent(in) :: h
    real(dp),dimension(n*2,n*2),intent(inout) :: A
    real(dp),dimension(n*2),intent(inout) :: B
    real(dp) :: alpha, aux1, aux2, aux3
    complex(dp) :: aux4, contorno_simpsom
    
    alpha=0.0_dp
    if (mod(m,4_dp)==0_dp) then
        alpha=4.0_dp
    end if
    
    A(1,1) = 4.0_dp
    A(1+n,1+n) = 4.0_dp
    A(1,2) = (-1.0_dp)*alpha
    A(1+n,2+n) = (-1.0_dp)*alpha
    do i=2,n,1
        aux1 = 1.0_dp/(h**2.0_dp)
        aux2 = 1.0_dp/(2.0_dp*(real(i-1, dp)*h)*h)
        aux3 = (real(m, dp)**2.0_dp)/((real(i-1, dp)*h)**2.0_dp)
        A(i,i-1) = (aux1-aux2)
        A(i+n,i+n-1) = (aux1-aux2)
        A(i,i) = ((-2.0_dp*aux1)-aux3)
        A(i+n,i+n) = ((-2.0_dp*aux1)-aux3)
        if (i/=n) then
            A(i,i+1) = (aux1+aux2)
            A(i+n,i+n+1) = (aux1+aux2)
        end if
        if (i==n) then
            aux4 = contorno_simpsom(m)
            B(n) = (-1.0_dp)*(aux1+aux2)*real(aux4, dp)
            B(2*n) = (-1.0_dp)*(aux1+aux2)*aimag(aux4)
        end if
    end do

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++        M e t o d o    S O R        ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
subroutine sor(A,B,n,error,w,logic,sol)
    use constants
    implicit none
    integer(dp),intent(in) :: n
    integer(dp) :: i,j,ite
    real(dp),intent(in) :: error,  w
    real(dp) :: delta,cuad,suma1,suma2,suma3
    real(dp),dimension(n,n),intent(in) :: A
    real(dp),dimension(n),intent(in) :: B
    real(dp),dimension(n),intent(inout) :: sol
    real(dp),dimension(n) :: aux,aux1,dif2
    logical,intent(in) :: logic
    ite=0 ; delta=100.0_dp

    call v_zeros(sol,n)
    call v_zeros(aux,n)
    call v_zeros(aux1,n)
    call v_zeros(dif2,n)

    do while (delta>=error)
        aux1=sol
        do i=1,n,1
            suma1 = 0.0_dp
            do j=1,i-1,1
                suma1 = suma1 + (A(i,j)*sol(j))
            end do
            suma2 = 0.0_dp
            do j=i+1,n,1
                suma2 = suma2 + (A(i,j)*sol(j))
            end do
            suma3 = 0.0_dp
            do j=1,i-1,1
                suma3 = suma3 + (A(i,j)*aux(j))
            end do
            aux(i)=(B(i)-((1.0_dp-w)*suma1)-suma2-(w*suma3))/A(i,i)
        end do
        sol=aux
        cuad = 0.0_dp
        dif2=(sol-aux1)
        do i=1,n,1
            cuad = cuad + (dif2(i)**2.0_dp)
        end do
        delta=(cuad**(0.5_dp))
    end do

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++    Completa Matriz Real con Ceros    +++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine m_zeros(A,n,m)
    use constants
    implicit none
    integer(dp),intent(in) :: n,m
    real(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j

    do i=1,n,1
        do j=1,m
            A(i,j) = 0.0_dp
        end do
    end do
    
    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++    Completa Vector Real con Ceros    +++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine v_zeros(A,n)
    use constants
    implicit none
    integer(dp),intent(in) :: n
    real(dp),dimension(n),intent(inout) :: A
    integer :: i

    do i=1,n,1
        A(i) = 0.0_dp
    end do
    
    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++    Completa Matriz Compleja con Ceros    +++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine m_zeros_complex(A,n,m)
    use constants
    implicit none
    integer(dp),intent(in) :: n,m
    complex(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j

    do i=1,n,1
        do j=1,m
            A(i,j) = (0.0_dp,0.0_dp)
        end do
    end do
    
    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++    Completa Vector Complejo con Ceros    +++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine v_zeros_complex(A,n)
    use constants
    implicit none
    integer(dp),intent(in) :: n
    complex(dp),dimension(n),intent(inout) :: A
    integer :: i

    do i=1,n,1
        A(i) = (0.0_dp,0.0_dp)
    end do
    
    return
end subroutine