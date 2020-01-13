!---------------------------------------------------------------------------------------------------!
!A simple model to simulate wind-driven ocean currents (Stommel,1948)
!Model structure: a main program "Main" and  1 subroutine: SWE2D
!Every change in Initial conditions or simulation setup must be decided in the Main program. 
!---------------------------------------------------------------------------------------------------!
PROGRAM Main

IMPLICIT NONE 
INTEGER, PARAMETER:: L=10E5, dx=20000, dy=20000
REAL, PARAMETER:: f0=10E-5, beta=10E-12, g=10, gm=10E-7, ro=1000, H=1000, tau_zero=0.2, pi= 3.1415927

!Setting simulation setup: 
INTEGER nx,ny,nt,dt,Sim_time,Steady_nt,i,j,t
REAL cg, C
LOGICAL forward_backward
PARAMETER (Sim_time= 60*60*24*25 , dt=100)! total length of simulation (seconds) and time step (seconds) 
PARAMETER (nx=L/dx, ny=L/dy, nt= Sim_time/dt)  !number of iteration required in space and time
REAL, DIMENSION(nx,ny)::  Eta_in, u_in, v_in
REAL,DIMENSION(nx,ny):: u_st, v_st, Eta_st
REAL,DIMENSION(nt,nx,ny)::  Eta, u, v 

write(*,*)'Running Fortran Code'
write(*,*)'... '


 cg=(g*H)**0.5 !gravity waves phase speed 
 C=cg*dt/dx !Courant number

 if(C.ge.(1/2**0.5))then !check if the CFL condition is respected 
   write(*,*) '<------ CFL condition is not respected! Simulation aborted ------>' 
   stop
 else 
  write(*,*) '<------ CFL condition is respected ------>' 
  write(*,*)  'Courant_number =', C
 endif

!Setting initial conditions (ICs)

do i=1,nx ! building a grid nx*ny at the initial time t=1 where all the values (eta, u and v) are zero 
do j=1,ny 
Eta_in(i,j)= 0
u_in(i,j)= 0
v_in(i,j)= 0
enddo
enddo

forward_backward=.true. !change it in .flase. to run the model with the semi-lagrangian scheme 

if(forward_backward)then 

   write(*,*) 'Time scheme used: forward-backward (Matsuno Sensei)'
  
   !solve the 2D Shallow Water Equations 
   call SWE2D(u_in,v_in,Eta_in,u,v,Eta,nx,ny,nt,dx,dy,dt, L,f0, beta, g, gm, ro, H, tau_zero, pi) 

   !calculate the Energy associated with the system
   call Energy(u,v,Eta,u_st, v_st, Eta_st,nx,ny,nt,Steady_nt,dx,dy,dt,g,ro,H,1) 
   
   !compute the analytical solution (to be compared with the numerical solution)
   call Analytical_solution(u_st,v_st,Eta_st,Eta,L,nt,Steady_nt,nx,ny,dx,dy, tau_zero, pi, H, gm, ro, f0, beta, g) 
 
   !calculate the Energy differences between model and analytical solution 
   call Energy(u,v,Eta,u_st,v_st,Eta_st,nx,ny,nt,Steady_nt,dx,dy,dt,g,ro,H,2) 
  

else 
  write(*,*) 'Time scheme used: Semi-lagrangian'

  !call Semilagrangian()

endif 

Stop

END PROGRAM Main 

!------------------------------------------------------------------------------------------------------------------!
!SUBROUTINE SWE2D: in this subroutine the 2D Shallow Water Equations are solved on the Arakawa-C grid, using a- 
!                  Foreward-Backward (Matsuno, 1966) time discretization scheme.
!INPUT values from Main program:physical parameters, Initial Conditions, size of the domain, grid spacing and time step
!OUTPUT variables: surface elevation Eta, u and v-component fields ( Eta(nt,nx,ny), u(nt,nx,ny), v(nt,nx,ny) ). 
!------------------------------------------------------------------------------------------------------------------!
SUBROUTINE SWE2D(u_in,v_in,Eta_in,u,v,Eta,nx,ny,nt,dx,dy,dt, L,f0, beta, g, gm, ro, H, tau_zero, pi)

IMPLICIT NONE 
INTEGER i,j,t, nx, ny, nt, dt, dx,dy, L,  ondeday_nt
REAL SumE, f0, beta, g, gm, ro, H, tau_zero, pi
REAL, DIMENSION(ny)::  tau_x, tau_y, f
REAL, DIMENSION(nx,ny)::  Eta_in, u_in, v_in, Eta_mid,  u_mid, v_mid, fv_mid, fu_mid, fv
REAL, DIMENSION(nt,nx,ny)::  Eta, u, v 

!N.B. The Foreward-Backward time scheme used requires that the equations are evaluated twice during each time step:
!at first, intermediate-values will be calculated; later, these values will be used to compute definitive-values 

DO t=1,nt  !start the TIME discretization 

!--- calculate INTERMEDIATE VALUES for the three prognostic variables (compute u before v) ---! 

!start the SPACE discretization using the staggered C-grid --> we are working on three different grids at the same time

!we are now on the Eta-grid (Eta is located at the center of each computational cell) 
  do i=1,nx-1
   do j=1,ny-1
     Eta_mid(i,j)= Eta_in(i,j)-H*dt*( ( (u_in(i+1,j)-u_in(i,j)) /dx ) + ( (v_in(i,j+1)-v_in(i,j)) /dy ) ) 
   enddo
  enddo

! we are now on the u-grid (u is located at the west boundary of each computational cell) 
  do i=2,nx
   do j=1,ny-1
    tau_x(j)= tau_zero*(-cos((pi*j*dy)/L))
    f(j)=(f0+beta*j*dy)
    fv_mid(i,j)= f(j)*((v_in(i,j)+v_in(i,j+1)+v_in(i-1,j)+v_in(i-1,j+1))/4.)
    u_mid(i,j)= u_in(i,j)+ fv_mid(i,j)*dt-g*dt*( (Eta_mid(i,j)-Eta_mid(i-1,j))/dx )-gm*dt*u_in(i,j)+ (tau_x(j)/(ro*H) )*dt
  enddo
 enddo 


! we are now on the v-grid (v is located at the south boundary of each computational cell) 
  do i=1,nx-1
   do j=2,ny
    f(j)=(f0+beta*j*dy)
    tau_y(j)= tau_zero*0.
    fu_mid(i,j)= f(j)*((u_mid(i,j)+u_mid(i+1,j)+u_mid(i,j-1)+u_mid(i+1,j-1))/4.)
    v_mid(i,j)= v_in(i,j)- fu_mid(i,j)*dt-g*dt*( (Eta_mid(i,j)-Eta_mid(i,j-1))/dy )-gm*dt*v_in(i,j)+ (tau_y(j)/(ro*H) )*dt
   enddo
  enddo


!-- now calculate DEFINITIVE VALUES for Eta u and v (compute v before u) ---!

! Eta-grid 
  do i=1,nx-1
   do j=1,ny-1
    Eta(t,i,j)= Eta_mid(i,j)-H*dt*( ((u_mid(i+1,j)-u_mid(i,j)) /dx ) + ( (v_mid(i,j+1)-v_mid(i,j)) /dy ) ) 
    Eta_in(i,j)=Eta(t,i,j)
   enddo
 enddo 

! v-grid
  do i=1,nx-1
   do j=2,ny
    if((j.eq.1).or.(j.eq.ny))then !setting no nornmal flow boundary conditions
    v(t,i,j)=0
    else
    v(t,i,j)= v_mid(i,j)- fu_mid(i,j)*dt-g*dt*( (Eta(t,i,j)-Eta(t,i,j-1))/dy )-gm*dt*v_mid(i,j)+ (tau_y(j)/(ro*H) )*dt
    endif
    v_in(i,j)=v(t,i,j)
  enddo
 enddo

! u-grid 
  do i=2,nx
   do j=1,ny-1
     if((i.eq.1).or.(i.eq.nx))then !setting no nornmal flow boundary conditions
     u(t,i,j)=0
     else
     f(j)=(f0+beta*j*dy)
     fv(i,j)= f(j)*( (v(t,i,j)+v(t,i,j+1)+v(t,i-1,j)+v(t,i-1,j+1))/4.)
     u(t,i,j)= u_mid(i,j)+ fv(i,j)*dt-g*dt*( (Eta(t,i,j)-Eta(t,i-1,j))/dx )-gm*dt*u_mid(i,j)+ (tau_x(j)/(ro*H) )*dt
     endif
     u_in(i,j)=u(t,i,j)
   enddo
 enddo


ENDDO


!Write variables in files in order to plot them using python scripts. 
!We are excluding the boundaries where the no nornmal flow boundary conditions have been applied. 

 open(unit=20, file="u_versus_x_1day") 
 open(unit=30, file="v_versus_y_1day")
 open(unit=40, file="Eta_elevation_1day") 
 open(unit=50, file="Eta_cross_1day") 

  ondeday_nt= 60*60*24/dt !we want to plot the Eta, u and v fields after one day of simulation, 
                          !thus we need to know which time step one day corresponds to 

  do i=2,nx-2  
    do j=2,ny-2
      write(40,*) i,j, Eta(ondeday_nt,i,j)
      write(50,*) (i-1)*dx,Eta(ondeday_nt,i,ny/2)
    enddo
 enddo


 do i=2,nx-1
    write(20,*)  (i-1)*dx, u( ondeday_nt,i,1)
 enddo

 do j=2,ny-1
  write(30,*)  (j-1)*dy, v( ondeday_nt,1,j)
 enddo


 close(20)
 close(30)
 close(40)
 close(50)


END SUBROUTINE SWE2D

!---------------------------------------------------------------------------------------------------------------------------------------------!
