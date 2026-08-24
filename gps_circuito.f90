program lectura_coordenadas
  real(8), allocatable :: latitud(:),longitud(:)
  real(8), allocatable :: dist(:),dir_final(:),curvat(:),df(:)
  integer, allocatable :: ipas(:)
  
  real(8) lat1, lon1, lat2, lon2, lat3, lon3,s_total1,s_total2
  real(8) R_curv, dir_deg(3)
  

  
  open(1,file='pista_gps.dat', status='unknown')
  read(1,*)ndatos
  allocate(latitud(ndatos),longitud(ndatos))
  do i=1,ndatos
     read(1,*)latitud(i),longitud(i)
  end do
  close(1)

  allocate(dist(ndatos), dir_final(ndatos),ipas(ndatos),curvat(ndatos),df(ndatos))
  open(1,file='pista_gps.out', status='unknown')
  dist(:)=0.d0
  dir_final(:)=0.d0
  ipas(:)=0
  
  do i=1,ndatos-2,2
     lat1=latitud(i); lon1=longitud(i)
     lat2=latitud(i+1); lon2=longitud(i+1)
     lat3=latitud(i+2); lon3=longitud(i+2)

     call analizar_tramo_curva(lat1, lon1, lat2, lon2, lat3, lon3, dir_deg, &
          R_curv, s_total1,s_total2)
     
     dist=dist+s_total
     if (R_curv.gt.1000) R_curv=1000

     do j=1,3
        ipas(i+j-1)=ipas(i+j-1)+1
        if (ipas(i+j-1).gt.1) then
           if(dabs(dir_final(i+j-1)-dir_deg(j)).gt.150) then
              dir_final(i+j-1)=dir_deg(j)*2
           else
              dir_final(i+j-1)=dir_final(i+j-1)+dir_deg(j)
           end if
        else
           dir_final(i+j-1)=dir_final(i+j-1)+dir_deg(j)
        end if
     end do
     
     dist(i+1)=dist(i)+s_total1
     dist(i+2)=dist(i+1)+s_total2
     curvat(i)=R_curv; curvat(i+1)=R_curv; curvat(i+2)=R_curv
  end do

  do i=1,ndatos
     dir_final(i)=dir_final(i)/ipas(i)
  end do
  
  rc=0.d0
  df(1)=dir_final(1)
  do i=2,ndatos
     if (dir_final(i)-dir_final(i-1).lt.-100.) then
        rc=-360.d0
     end if
     if (dir_final(i)-dir_final(i-1).gt.100.) then
        rc=360.d0
     end if     
     df(i)=dir_final(i)-rc
  end do
  
  open(1,file='pista_gps.out', status='unknown')
  write(1,*) ndatos
  do i=1,ndatos
     write(1,*) dist(i),df(i) ,0.0d0, curvat(i)
  end do
  close(1)
end program lectura_coordenadas


subroutine analizar_tramo_curva(lat1, lon1, lat2, lon2, lat3, lon3, &
                               dir_deg, R_curv, s_total1, s_total2)
    implicit none
    ! Entradas (Lat/Lon de Google Maps)
    real(8), intent(in)  :: lat1, lon1, lat2, lon2, lat3, lon3
    
    ! Salidas
    real(8), intent(out) :: dir_deg(3)     ! Dirección tangente en P1, P2 y P3 (0° a 360°)
    real(8), intent(out) :: R_curv         ! Radio de curvatura (m)
    real(8), intent(out) :: s_total1,s_total2   ! Longitud total del arco P1 -> P3 (m)
    real(8) s_total
    real(8) A1,A2,A3

    
    character(len=10)  sentido ! 'Izquierda', 'Derecha', 'Recta'

    real(8), parameter :: DEG2RAD = 3.14159265358979323846_8 / 180.0_8
    real(8), parameter :: R_EARTH = 6371000.0_8

    real(8) :: lat_m, x(3), y(3)
    real(8) :: a, b, c, D, Dc
    real(8) :: xc, yc, alpha, dir_recta, theta12, theta23
    integer :: i

    ! 1. Proyección plana local con origen en P1
    lat_m = 0.5_8 * (lat1 + lat3) * DEG2RAD
    
    x(1) = 0.0_8
    y(1) = 0.0_8
    
    x(2) = R_EARTH * (lon2 - lon1) * DEG2RAD * cos(lat_m)
    y(2) = R_EARTH * (lat2 - lat1) * DEG2RAD
    
    x(3) = R_EARTH * (lon3 - lon1) * DEG2RAD * cos(lat_m)
    y(3) = R_EARTH * (lat3 - lat1) * DEG2RAD

    ! 2. Longitudes de las cuerdas (rectas)
    a = sqrt((x(2)-x(1))**2 + (y(2)-y(1))**2)  ! P1 -> P2
    b = sqrt((x(3)-x(2))**2 + (y(3)-y(2))**2)  ! P2 -> P3
    c = sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2)  ! P1 -> P3

    ! Producto cruz 2D (Doble del área orientada)
    D = (x(2) - x(1))*(y(3) - y(1)) - (y(2) - y(1))*(x(3) - x(1))
    write(*,*)  D
    ! 3. Evaluación Recta vs Curva
    if (abs(D) < 1.0e+1_8) then
        ! Tramo en Recta
        R_curv = 1.0e6_8
        sentido = 'Recta'
        s_total = a + b  ! Suma directa de segmentos
        s_total1=a; s_total2=b
        
        ! Dirección uniforme desde P1 hacia P3
        dir_recta = atan2(x(3) - x(1), y(3) - y(1)) / DEG2RAD
        dir_recta = mod(dir_recta + 360.0_8, 360.0_8)
        dir_deg(:) = dir_recta
    else
        ! Tramo en Curva
        R_curv = (a * b * c) / (2.0_8 * abs(D))
        
        if (D > 0.0_8) then
            sentido = 'Izquierda'
        else
            sentido = 'Derecha'
        endif

        ! 4. Longitud Total sobre el Arco
        theta12 = 2.0_8 * asin(min(1.0_8, a / (2.0_8 * R_curv)))
        theta23 = 2.0_8 * asin(min(1.0_8, b / (2.0_8 * R_curv)))
        s_total = R_curv * (theta12 + theta23)
        s_total1=R_curv * (theta12); s_total2=R_curv * (theta23)
        
        ! 5. Centro de Curvatura y Direcciones Tangentes
        Dc = 2.0_8 * D
        xc = ((x(1)**2 + y(1)**2)*(y(2) - y(3)) + &
              (x(2)**2 + y(2)**2)*(y(3) - y(1)) + &
              (x(3)**2 + y(3)**2)*(y(1) - y(2))) / Dc
              
        yc = ((x(1)**2 + y(1)**2)*(x(3) - x(2)) + &
              (x(2)**2 + y(2)**2)*(x(1) - x(3)) + &
              (x(3)**2 + y(3)**2)*(x(2) - x(1))) / Dc

        do i = 1, 3
            alpha = atan2(x(i) - xc, y(i) - yc) / DEG2RAD
            if (D > 0.0_8) then
                dir_deg(i) = alpha - 90.0_8  ! Izquierda
            else
                dir_deg(i) = alpha + 90.0_8  ! Derecha
            endif
            dir_deg(i) = mod(dir_deg(i) + 360.0_8, 360.0_8)
        end do
     endif
     A1=90-dir_deg(1); if (A1.lt.0.d0) a1=a1+360
     A2=90-dir_deg(2); if (A2.lt.0.d0) a2=a2+360
     A3=90-dir_deg(3); if (A3.lt.0.d0) a3=a3+360
     dir_deg(1)=A1; dir_deg(2)=A2; dir_deg(3)=A3
     
     write(*,*) sentido,a1,a2,a3
    
end subroutine analizar_tramo_curva

