module auto
  real(8) masa
  real(8) Cd, A_ref
  real(8) R_rueda,Ri
  real(8) E_disp,P_max
  integer, parameter :: n=200, m=1000

  real(8) :: kl1=2.5d0, kl2=0.002D0, kl3=8.d0    !cálculo de eficiencia del conjunto
!  real(8) :: kl1=0.d0, kl2=0.0d0, kl3=0.d0    !cálculo de eficiencia del conjunto
  real(8) :: eta_transm=0.97d0
  real(8) :: relacion=7.d0, radio=0.25d0       !cálculo de eficiencia del conjunto
  
  real(8) C_alfa,Pres                                  !calculo de C_alfa
  real(8) :: P_ref=400000d0, coef_a=3.2d0, coef_b=0.4d0 !calculo de C_alfa
  !4.5 bares de referencia son ~65PSI
  real(8) Crr                                                       !calculo de Crr
  real(8) :: Crr_ref=0.0072d0, coef_n=0.4d0, Fact_diametro=0.901d0 !calculo de Crr
  
end module auto

module viento
  real(8) V_w,phi_w,rho
end module viento

module circuito
  integer ndatos
  real(8),allocatable:: SD(:),psi(:),theta(:),RD(:)
  real(8) Lcirc,T_carrera
end module circuito


Program derivadas
  use circuito
  use auto
  
  real(8) P_bar,s(m+1),ds
  complex(8) P(n),PP(n),Time,eps,P_media,V(m+1)
  real(8) g(n),gf(n),gg(n),gq(n)
  real(8) mod_gg,gfxgg,time_old
  real(8) interpolacion,altura,distx,disty,direc
  real(8) max_up,sumg
  
  real(8) alfa, alfa_min, beta_ls
  real(8) time_actual, time_trial
  real(8) P_media_actual, P_media_trial
  complex(8) P_trial(n), V_trial(m+1)
  logical aceptado
  integer ils, max_ls
  
  call lectura()

  C_alfa=coef_a*(masa*9.81d0)/4.d0*(P_ref/Pres)**(-coef_b)
  Crr=Crr_ref*(P_ref/Pres)**coef_n*Fact_diametro
    
  ds=Lcirc/m
  do j=1,m
     S(j)=(j-1)*ds
  end do  
  S(m+1)=Lcirc
  
  fact=1
  do i=1,n
     P(i)=E_disp/T_carrera*fact
     if (real(P(i)).gt.P_max) P(i)=P_max
  end do
  
  alfa     = 100000.d0
  alfa_min = 1.d-6
  beta_ls  = 0.5d0
  max_ls   = 20
  call Funcion(Time,P_media,P_bar,P,V)
  time_actual     = real(Time)
  P_media_actual  = real(P_media)
  time_old        = time_actual
  
  do kk=1,3000
     do i=1,n
        PP(:)=P(:)
        eps=(0,1d-10)
        PP(i)=P(i)+eps
        call Funcion(Time,P_media,P_bar,PP,V)
        
        gf(i)=aimag(Time)/1d-10
        g(i)=-(Real(P_media)-P_bar)
        gg(i)=-aimag(P_media)/1d-10
     end do
     
     mod_gg=0.d0; gfxgg=0.d0 
     do i=1,n
        mod_gg=mod_gg+gg(i)**2
        gfxgg=gfxgg+gf(i)*gg(i)
     end do
     
     do i=1,n
        gq(i)=gf(i)-gfxgg/mod_gg*gg(i)
     end do
     
     !line search
     !-----------
     aceptado=.false.
     alfa=max(alfa,alfa_min)
     
     do ils=1,max_ls
        p_trial(:)=p(:)-alfa*gq(:)
        
        do i=1,n
           
           if (real(P_trial(i)).gt.P_max) then
              P_trial(i)=cmplx(P_max,0.d0,kind=8)
           else if (real(P_trial(i)).lt.0.d0) then
              P_trial(i)=cmplx(0.d0,0.d0,kind=8)
           end if
           
        end do
             
        do j=1,3
           do i=1,n
              PP(:)=P_trial(:)
              eps=(0,1d-10)
              PP(i)=P_trial(i)+eps
              call Funcion(Time,P_media,P_bar,PP,V)
              
              g(i)=-(Real(P_media)-P_bar)
              gg(i)=-aimag(P_media)/1d-10
           end do
           
           mod_gg=0.d0
           do i=1,n
              mod_gg=mod_gg+gg(i)**2
           end do
           
           do i=1,n
              P_trial(i)=P_trial(i)-g(i)/mod_gg*gg(i)
           end do
           
           do i=1,n
              if (real(P_trial(i)).gt.P_max) then
                 P_trial(i)=cmplx(P_max,0.d0,kind=8)      
              else if (real(P_trial(i)).lt.0.d0) then
                 P_trial(i)=cmplx(0.d0,0.d0,kind=8)           
              end if
           end do
           
        end do
        call Funcion(Time,P_media,P_bar,P_trial,V)
        time_trial=real(Time)
        P_media_trial=real(P_media)
        
        if (time_trial.lt.time_actual) then
           write(*,*) 'fin del cálculo'
           aceptado=.true.
           exit
        end if
        alfa=beta_ls*alfa
        if(alfa.lt.alfa_min) exit
        
     end do
     
     if(aceptado) then
        p(:)=P_trial(:)
        
        time_old=time_actual
        time_actual=time_trial
        
        P_media_actual=p_media_trial
        alfa=min(1.2d0*alfa,100000.d0)
     else
        write(*,*) 'Line search no encontro mejora'
        write(*,*) 'Iteracion:',kk,' Alfa:',alfa
        exit
     end if

     max_up=0.d0
     do i=1,n
        if (real(P(i)).gt.max_up) max_up=real(P(i))
     end do
     
     !  write(*,*) 'Tiempo:',real(Time),'  Cambio de P()',max_up
     write(*,'(A,I5,A,F14.8,A,F12.4,A,F14.6)')               &
          'Iter:',kk,                                          &
          '  Time=',time_actual,                               &
          '  Alfa=',alfa,                                      &
          '  Pmedia-Pbar=',P_media_actual-P_bar 
     
     if (dabs(time_actual-time_old).lt.1.d-4) then
        write(*,*) 'fin del calculo'
        exit        
     end if
     
     !   write(*,*) 'Pasos:',kk,'   Tiempo:',real(Time),'  Alfa=',alfa
  end do


  open(10,file='datos.json',status='unknown')
  write(10,'(a)') '['
  do j=1,m-1
     write(10,'(A,f14.2,A,f14.2,A,f14.2,A)') &
          '     {"metros":',s(j),', "porcentaje":',real(P(1 + (j-1)*n/m)),', "v_ideal:',real(V(j))*3.6d0,'},'
  end do
  j=m
  write(10,'(A,f14.2,A,f14.2,A,f14.2,A)') &
       '     {"metros":',s(j),', "porcentaje":',real(P(1 + (j-1)*n/m)),', "v_ideal:',real(V(j))*3.6d0,'}'
  write(10,'(a)') ']'
  close(10)
  
  open(8,file='Circuito.dat',status='unknown')
  open(7,file='Res_optimizado.dat',status='unknown')  
  altura=0.D0
  distx=0.d0; disty=0.d0
  direc=interpolacion(SD,psi,ndatos,0.d0)
  write(8,*) distx,disty,direc

  do j=1,m
     altura=altura+ds*dsin(interpolacion(SD,theta,ndatos,S(j)))
     direc=interpolacion(SD,psi,ndatos,S(j))
     distx=distx+ds*dcos(direc); disty=disty+ds*dsin(direc)
 
     write(8,*) distx,disty,direc,real(P(1 + (j-1)*n/m)),real(V(j))
     write(7,*) s(j),altura,real(P(1 + (j-1)*n/m)),real(V(j))
  end do
  close(8)
  close(7)
  
end Program derivadas


Subroutine Funcion(Time,P_media,P_bar,P,V)
  use circuito
  use auto
  implicit none
  
  integer i,idx,j,k,kk
  real(8) S(m+1)
  complex(8) V(m+1)
  real(8) F_rod,F_pend,F_aero_p
  complex(8) F_aero,F_tot,F_trac,F_curvas,f_j
  real(8) ds,P_bar
  complex(8) P_media,E,P_el
  real(8) interpolacion
  complex(8) P(n),Time


! posiciones en el circuito  
  ds=Lcirc/m

  do j=1,m
     S(j)=(j-1)*ds
  end do  
  S(m+1)=Lcirc

  P_bar=E_disp/T_carrera
  
  V(:)=(10.d0,0.d0)

  do  !calculo en la pista!    
     V(1)=V(m+1)
     do j=1,m
        idx = 1 + (j-1)*n/m
        call fuerzas_aero(F_aero,V(j),S(j))
        call fuerzas_rodadura(F_rod,S(j))
        call fuerzas_pendiente(F_pend,S(j))
        call fuerzas_traccion(F_trac,V(j),j,P)
        call fuerzas_curvas(F_curvas,V(j),S(j))
        
        F_tot=F_trac-(F_aero+F_rod+F_pend+F_curvas)
        f_j=1.d0/(masa*V(j))*F_tot
        
        V(j+1)=V(j)+ds*f_j
        
     end do
     
     if (dabs(Real(V(m+1))-Real(V(1))).lt.1d-5) exit
  end do
  
  Time=(0.d0,0.d0)
  do j=1,m
     Time=Time+ds/V(j)
  end do

  E=(0.d0,0.d0)
  do j=1,m
     idx = 1 + (j-1)*n/m
     P_el = P(idx)
     E=E+P_el/V(j)*ds
  end do
  P_media=E/Time
    
end Subroutine Funcion


subroutine fuerzas_aero(F_aero,V,S)
  use auto !Cd,A_ref
  use viento !V_w,phi_w, rho
  use circuito !psi(i),ndatos

  real(8) S
  complex(8) F_aero,V,V_real
  real(8) interpolacion
  real(8) phi_s
  
  phi_s=phi_w-interpolacion(SD,psi,ndatos,S)
  V_real=V-V_w*dcos(phi_s)
  F_aero=0.5d0*rho*Cd*A_ref*V_real**2
  
end subroutine fuerzas_aero

subroutine fuerzas_rodadura(F_rod,S)
  use circuito !psi(i),theta(i),ndatos
  use auto !Crr, masa
  
  real(8) F_rod,S
  real(8) :: g=9.81D0
  real(8) interpolacion
  
  F_rod=Crr*Masa*g*dcos(interpolacion(SD,theta,ndatos,S))
  
end subroutine fuerzas_rodadura

subroutine fuerzas_pendiente(F_pend,S)
  use circuito
  use auto
  
  real(8) F_pend,S
  real(8) :: g=9.81D0
  real(8) interpolacion
  
  F_pend=masa*g*dsin(interpolacion(SD,theta,ndatos,S))

end subroutine fuerzas_pendiente

subroutine fuerzas_curvas(F_curvas,V,S)
  use circuito !psi(i),theta(i),ndatos
  use auto !Crr, masa
  
  complex(8) F_curvas,V
  real(8) interpolacion2,R_eff,S
  R_eff=interpolacion2(SD,RD,ndatos,S)
  F_curvas=masa**2*V**4/(C_alfa*R_eff**2)
  
end subroutine fuerzas_curvas

subroutine fuerzas_traccion(F_trac,V,j,P)
  use auto
  use circuito

  integer j
  complex(8) F_trac,V,Pot,P(n)
  complex(8) Torque,omega,eta

  Pot=P(1+INT((j-1)*n/m))
  omega=V/radio*relacion
  Torque=Pot/omega
  eta=Torque*omega/(Torque*omega+kl1*Torque**2+kl2*omega+kl3)
  
  F_trac=eta*eta_transm*Pot/V
  
end subroutine fuerzas_traccion
        

subroutine lectura()
  use auto
  use viento
  use circuito

  open(1,file='datos.dat',status='old')
  read(1,*)
  read(1,*) masa
  read(1,*) Cd, A_ref
  read(1,*) Pres
  Pres=Pres*6894.75729 !pasa de PSI a N/m^2
!  read(1,*) E_disp,P_max
  read(1,*) P_max
  Write(*,'(A$)') "Energia disponible [Joules]:"
  read(*,*) E_disp
  
  read(1,*)
  read(1,*) rho
  read(1,*) V_w,phi_w
  phi_w=phi_w*3.14159265d0/180.d0

  read(1,*)
!  read(1,*) Lcirc,T_carrera
  read(1,*) Lcirc
  Write(*,'(A$)') "Tiempo restante [min]:"
  read(*,*) T_carrera
  T_carrera=T_carrera*60
    
  read(1,*) ndatos
  allocate(SD(ndatos),psi(ndatos),theta(ndatos),RD(ndatos))

  do i=1,ndatos
     read(1,*) SD(i),psi(i),theta(i),RD(i)
     psi(i)= psi(i)*3.14159265d0/180.d0
     theta(i)=theta(i)*3.14159265d0/180.d0
  end do
  close(1)

end subroutine lectura


function interpolacion(X,Y,ndatos,X0)
  integer ndatos
  real(8) interpolacion
  real(8) X(ndatos),Y(ndatos)
  real(8) X0

  if(X0.le.X(1)) then
     interpolacion=Y(1)
     return
  else if(X0.ge.X(ndatos)) then
     interpolacion=Y(ndatos)
     return
  end if
  do i=1,ndatos-1
     if(X0.ge.X(i).and.X0.le.x(i+1)) then
        interpolacion=Y(i)+(Y(i+1)-Y(i))/(X(i+1)-X(i))*(X0-X(i))
        return
     end if
  end do
 
end function interpolacion
  
function interpolacion2(X,Y,ndatos,X0)
  integer ndatos
  real(8) interpolacion2
  real(8) X(ndatos),Y(ndatos)
  real(8) X0

  if(X0.le.X(1)) then
     interpolacion2=Y(1)
     return
  else if(X0.ge.X(ndatos)) then
     interpolacion2=Y(ndatos)
     return
  end if
  do i=1,ndatos-1
     if(X0.ge.X(i).and.X0.le.x(i+1)) then
        interpolacion2=Y(i)
        return
     end if
  end do
 
end function interpolacion2
  
