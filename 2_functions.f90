!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++    R E A L   F U N C T I O N S     ++++++++++++++++
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
!++++++++++++++   Condicion de Contorno de Dirichlet   ++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Necesario para obtener una solucion.
!
function funcion_contorno(x) result(aux)
    use constants
    implicit none
    real(dp),intent(in) :: x
    real(dp) :: aux

    aux = sin(x)
    !aux = (1.0_dp+(x**2_dp))*sin(x)
    !aux = 5.0_dp*sin(x) + 10
    !aux = 20.0_dp + (0.0_dp*x)
    !aux = (2.0_dp/3.0_dp)*( ((5.0_dp*sin(x)+10.0_dp)**(1.5_dp))/(1.0_dp**0.5_dp) - 1.0_dp )

    return
end function
!
!
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++                                          +++++++++++++
!+++++++++++++    C O M P L E X   F U N C T I O N S     +++++++++++++
!+++++++++++++                                          +++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++   Coeficiente (C_m)  de Fourier    ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Aplicando metodo de integración de Simpson (1/3) compuesto.
!
function contorno_simpsom(m) result(coef)
    use constants
    implicit none
    integer(dp),intent(in) :: m
    integer(dp) :: n, i 
    real(dp) :: fun, funcion_contorno, I_re, I_im
    real(dp) :: aux, aux1_y1, aux2_y1, aux1_y2, aux2_y2
    real(dp), dimension(:), allocatable :: x_vec, y1_vec, y2_vec
    complex(dp) :: coef

    n=201
    allocate(x_vec(n))
    allocate(y1_vec(n))
    allocate(y2_vec(n))

    call linspace(0.0_dp, 2.0_dp*pi, x_vec, n)
    do i=1,n,1
        fun = funcion_contorno(x_vec(i))
        y1_vec(i) = cos(real(m, dp)*x_vec(i)) * fun 
        y2_vec(i) = sin(real(m, dp)*x_vec(i)) * fun
    end do
    
    aux1_y1 = 0.0_dp
    aux2_y1 = 0.0_dp
    aux1_y2 = 0.0_dp
    aux2_y2 = 0.0_dp

    do i = 2,n-1,1
        if (mod(i,2)==0) then
            aux1_y1 = aux1_y1 + y1_vec(i)
            aux1_y2 = aux1_y2 + y2_vec(i)
        else
            aux2_y1 = aux2_y1 + y1_vec(i)
            aux2_y2 = aux2_y2 + y2_vec(i)
        end if
    end do
    
    aux = (x_vec(n)-x_vec(1)) / (3.0_dp*real(n-1, dp))
    I_re = aux * ( y1_vec(1) &
        & + (4.0_dp*aux1_y1) &
        & + (2.0_dp*aux2_y1) &
        & + y1_vec(n) )
    I_im = aux * ( y2_vec(1) &
        & + (4.0_dp*aux1_y2) &
        & + (2.0_dp*aux2_y2) &
        & + y2_vec(n) )

    coef = complex(I_re/(2.0_dp*pi), -I_im/(2.0_dp*pi))
    !write(*,*) coef

    return
end function
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++    Función temperatura en (r_p,th)     ++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Se necesita la matriz con coeficientes (coef) y el valor de m_max.
! mp : m_max + 1
!
function temperatura(coef,mp,n,p,th) result(t)
    use constants
    implicit none
    integer(dp) :: i
    integer(dp),intent(in) :: mp,n,p
    real(dp),intent(in) :: th
    real(dp) :: t,cz,sz
    complex(dp),dimension(mp,n),intent(in) :: coef
    complex(dp) :: aux,z,C
    
    aux=coef(1,p)
    do i=2,mp,1
        C = coef(i,p)
        cz = cos(real(i-1, dp)*th)
        sz = sin(real(i-1, dp)*th)
        z = complex(cz,sz)
        aux=aux+(C*z)+(conjg(C)*conjg(z))
    end do
    t= real(aux)

    return
end function