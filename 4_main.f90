!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++      M A I N   P R O G R A M       ++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Solucion Numerica para 'n_r' y 'm_max'
!
!
program main
    use constants
    implicit none
    integer(dp) :: i,j,l
    integer(dp) :: r,n,m,mp,pa
    real(dp) :: h,rd,th,p_th
    real(dp) :: temperatura, funcion_contorno
    complex(dp),dimension(:,:),allocatable :: coef
    real(dp),dimension(:),allocatable :: sol_aux
    

    write(*,*) '---------------------------------------------'
    write(*,*) '--------------- PROGRAMA LSPC ---------------'
    write(*,*) '---------------------------------------------'
    write(*,*) ' '
    write(*,*) 'Ingrese el valor de n_r:'
    read(*,*) r
    write(*,*) 'Ingrese el valor de m_max:'
    read(*,*) m
    
    n = r+1
    mp = m+1
    pa=101 

    allocate(coef(mp,n))
    allocate(sol_aux(pa))

    call get_coef(n,mp,pa,h,p_th,coef)
    
    open(20,file='file1.dat',status='replace')
    open(30,file='file2.dat',status='replace')
    open(40,file='file3.dat',status='replace')

    write(40,*) '#  Matriz de temperaturas T[i,j] ---> T(r_i,th_j)'
    do i=1,n,1
        rd= h*real(i-1, dp)
        write(20,*) rd
        do j=0,pa-1,1
            th=p_th*real(j, dp)
            if (i==1) then
                write(30,*) th
            end if
            sol_aux(j+1) = temperatura(coef,mp,n,i,th)
        end do
        !write(40,*) sol_aux
        write(40,'(999(F11.6,X))') (sol_aux(l), l=1,pa)
    end do

    write(20,*) h*real(n, dp)
    do j=0,pa-1,1
        th=p_th*real(j, dp)
        sol_aux(j+1) = funcion_contorno(th)
    end do
        !write(40,*) sol_aux
        write(40,'(999(F11.6,X))') (sol_aux(l), l=1,pa)
    close(20)
    close(30)
    close(40)
    
end program main
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++