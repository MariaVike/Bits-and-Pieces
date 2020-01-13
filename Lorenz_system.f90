!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This program solves the 3D Lorenz system using the 4th-order Runge-Kutta time integration scheme. The user can decide whether 
!running a single experiment with fixed IC, or running an ensemble of size X whose ICs are generated using the random number
!generator function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM LorenzSystem
IMPLICIT NONE

INTEGER k, nt, Sim_time
INTEGER, PARAMETER:: nk=1!10
REAL, DIMENSION(nk):: x, y, z
REAL dt, t 


INTERFACE
REAL FUNCTION Random(x,y)
 real,intent(in) :: x,y
END FUNCTION Random
END INTERFACE

write(*,*)'Running Fortran Code'
write(*,*)'... '


!Setting simulation setup: 

Sim_time=100   ! total length of simulation (seconds)
dt=0.01                 ! time step (seconds) 
nt= Sim_time/dt        ! number of iterations


! Setting Initial Conditions and calling Runge-Kutta 4th order scheme:
  
  do k=1,nk    !set (above) nk=1 to run only the given initial condition, 
               !set nk=forecast_number to run the ensemble 
   
   if(k.eq.1)then
   x(k)= 1.
   y(k)= 1.
   z(k)= 1.
   t   = 0
   else 
   x(k)=Random(-5. , 5.0)   !use the random function to generate random values 
   y(k)=Random(-10. , 10.0) !of x,y and z 
   y(k)=Random(-10. , 10.0) 
                         
   endif 

   call RK4th(x, y, z, t, dt, nt, k)
   
  enddo 

END PROGRAM LorenzSystem

!---------------------------------------------------------!

 REAL FUNCTION f1(x, y, t)
      IMPLICIT NONE
      real, intent(in) :: x,y,t
      real, parameter::  alfa_1 = 10,  alfa_2 = 20
      real R 
     
      call Random_number(R)
      f1 = alfa_1*(y-x) +R
      write(*,*) R
 END FUNCTION f1

REAL FUNCTION f2(x, y, z, t)
      IMPLICIT NONE
      real, intent(in) :: x,y,z,t
      real, parameter::  ghamma_1=28, ghamma_2=46.92
      
      f2 = ghamma_1*x-y-x*z

 END FUNCTION f2

REAL FUNCTION f3(x, y, z, t)
      IMPLICIT NONE
      real, intent(in) :: x,y,z,t
      real, parameter::  beta_1=2.66,  beta_2=4
      
      f3 = x*y-beta_1*z

 END FUNCTION f3

!------------------------------------------------

SUBROUTINE RK4th(Q1_in, Q2_in, Q3_in, t, dt, nt, k) 

IMPLICIT NONE

INTEGER  nt, i, k
REAL, DIMENSION(0:nt)::  Q1_out, Q2_out, Q3_out
REAL dt, t, Q1_in, Q2_in, Q3_in, k1, k2, k3,k4, j1, j2, j3, j4, m1, m2, m3, m4, Rnd_num
integer, dimension(1):: seed=1
     


     INTERFACE
	Real Function f1(Q1_in, Q2_in, t)
	 real,intent(in) :: Q1_in, Q2_in, t
	End Function f1
    
	Real Function f2(Q1_in, Q2_in, Q3_in, t)
	 real,intent(in) ::  Q1_in, Q2_in, Q3_in, t
	End Function f2
 
        Real Function f3(Q1_in, Q2_in, Q3_in, t)
	 real,intent(in) ::  Q1_in, Q2_in, Q3_in, t
	End Function f3

        Real Function Random(x,y)
	 real,intent(in) ::  x,y
	End Function Random
     END INTERFACE


         open(unit=20+k)

  call random_seed(PUT=seed)
  
      do i=0,nt-1


           k1 = f1(Q1_in,Q2_in, t)   
           j1 = f2(Q1_in,Q2_in,Q3_in, t)
           m1 = f3(Q1_in,Q2_in,Q3_in, t)          
       
           k2 = f1(Q1_in+0.5*k1*dt, Q2_in+0.5*j1*dt, t+ 0.5*dt) 
           j2 = f2(Q1_in+0.5*k1*dt, Q2_in+0.5*j1*dt, Q3_in+0.5*m1*dt, t+ 0.5*dt) 
           m2 = f3(Q1_in+0.5*k1*dt, Q2_in+0.5*j1*dt, Q3_in+0.5*m1*dt, t+ 0.5*dt)  
         
           k3 = f1(Q1_in+0.5*k2*dt, Q2_in+0.5*j2*dt, t+ 0.5*dt)   
           j3 = f2(Q1_in+0.5*k2*dt, Q2_in+0.5*j2*dt, Q3_in+0.5*m2*dt, t+ 0.5*dt)
           m3 = f3(Q1_in+0.5*k2*dt, Q2_in+0.5*j2*dt, Q3_in+0.5*m2*dt, t+ 0.5*dt)
          
           k4 = f1(Q1_in+k3*dt, Q2_in+j3*dt, t+ dt)
           j4 = f2(Q1_in+k3*dt, Q2_in+j3*dt, Q3_in+m3*dt, t+ dt)
           m4 = f3(Q1_in+k3*dt, Q2_in+j3*dt, Q3_in+m3*dt, t+ dt)
         
           Q1_out(i) = Q1_in+dt*(k1+2*k2+2*k3+k4)/6.
           Q2_out(i) = Q2_in+dt*(j1+2*j2+2*j3+j4)/6.
           Q3_out(i) = Q3_in+dt*(m1+2*m2+2*m3+m4)/6.
            
            write(20+k,*) t, Q1_out(i), Q2_out(i), Q3_out(i)
 
            Q1_in=Q1_out(i)
            Q2_in=Q2_out(i)
            Q3_in=Q3_out(i)
            t=t+dt
        

     enddo

        close(unit=20+k)        

END SUBROUTINE RK4th
