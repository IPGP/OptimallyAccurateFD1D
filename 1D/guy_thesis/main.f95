program main
!
!             1d FWI code using the following Forward modelling operators
!                 - optimally accurate operators O(2,2)
!                 - optimally accurate operators O(2,4)
!                 - conventional operators O(2,2)
!                 - conventional operators O(2,4)
!                 - staggered grid first order derivatives O(2,2)
!                           - Guy Mabyalaht                    
!                              12-March-2016
!
!---------------------------------------------------------------------------
!  	                  MODELLING SETTINGS
!---------------------------------------------------------------------------
  use mod_FW
  use mod_data
  use mod_pml
  !
  implicit none
  !
  !
  ! variables to store the choice of operators
  logical :: CONV2, OPT2, CONV4, OPT4
  integer :: choice

  ! variables used in the FWI-----------------------------------------------------
  integer, parameter :: dp=kind(1.0d0)
  real(dp), dimension(:,:), allocatable :: u_field, dobs
  real(dp), dimension(:), allocatable :: time_vec, x_vec, sig, cpe  
  real(dp), dimension(:), allocatable :: rho, mu, u, u1, cp, cs
  integer ::  i, NSTEP, ISOURCE, npml, k=1, NX, NX1, IT_DISPLAY, mg
  real(dp):: time=0.08d0, DELTAT, XMAX=100.0d0, XMIN=0.d0, DELTAX, &
  & xsource, mindist, courant_number1, vm, a=2.5, b=0.5, vmin, &
  & courant_number2, dx2, dt2, dist
  real(dp), dimension(:,:,:), allocatable :: A_conv, A_opt, K_opt, &
        K_conv, dA, dK,UI, A_conv4, A_opt4, K_conv4, K_opt4, dA4, dK4, UI4
  !--------------------------------------------------------------------------------

  ! Displays the Modelling operators availabe in the FWI code
  write(*,*) 'Welcome to the 1D FWI code ' 
  write(*,*) 'Here are the options for the modelling operators to be used'
  write(*,*) 'Type 1 for the Opt(2,2)'
  write(*,*) 'Type 2 for the Opt(2,2)'
  write(*,*) 'Type 3 for the Conv(2,2)'
  write(*,*) 'Type 4 for the Conv(2,4)'
  write(*,*) 'Type 5 for the staggered grid'
  write (*,*) 'the choice must be defined in the input_data.don file'
  !
   
  ! reading in the choice of the user
  read(*,*) choice
  !
  ! reading in the model parameters
  read(*,*) NX     ! model size
  read(*,*) NX1    ! model size
  read(*,*) vmin   ! vmin
  read(*,*) DELTAT ! time step
  read(*,*) mg
  !
  !
  allocate(rho(NX+1), mu(NX+1), u(NX+1), u1(NX+1), cp(NX+1), cs(NX+1))
  !

  !
  
  ! selects the operator to be used
  select case (choice)
  ! Opt(2,2)
  case (1)
    write(*,*) 'opt(2,2) has been selected'
    CONV2 = .false.
    CONV4 = .false.
    OPT2  = .true.
    OPT4  = .false.

  ! Opt(2,4)
  case(2) 
   write(*,*) 'opt(2,4) has been selected'
    CONV2 = .false.
    CONV4 = .false.
    OPT2  = .false.
    OPT4  = .true.
  
  ! Conv(2,2)
  case(3) 
   write(*,*) 'conv(2,2) has been selected'
    CONV2 = .true.
    CONV4 = .false.
    OPT2  = .false.
    OPT4  = .false.

  ! Conv(2,4)
  case(4)  
   write(*,*) 'conv(2,4) has been selected'
    CONV2 = .false.
    CONV4 = .true.
    OPT2  = .false.
    OPT4  = .false.
  
  ! No choice has been done
  case default
    write(*,*) 'Choice out of range, the code will stop'
    stop
  end select
 !
 
 allocate (A_conv(1:3,1:3,NX-1), A_opt(1:3, 1:3, NX-1), K_opt(1:3, 1:3, NX-1), &
             & K_conv(1:3, 1:3, NX-1), dA(1:3, 1:3, NX-1), dK(1:3, 1:3, NX-1), UI(1:3, 1:3, NX-1))
 

 allocate (A_conv4(1:3,1:5,NX-1), A_opt4(1:3, 1:5, NX-1), K_opt4(1:3, 1:5, NX-1), &
             & K_conv4(1:3, 1:5, NX-1), dA4(1:3, 1:5, NX-1), dK4(1:3, 1:5, NX-1), UI4(1:3, 1:5, NX-1))
 !----------------------------------------------------------------------------
 
 ! time step in seconds
 DELTAT=DELTAT/mg

 ! total number of time steps
   NSTEP= nint(time/DELTAT)

 ! time axis
   allocate(time_vec(NSTEP))
   do i=1,NSTEP
   time_vec(i)=i*DELTAT
   end do
 !
 allocate(dobs(NX+1, NSTEP))  

 ! size of a grid cell
   DELTAX=(XMAX-XMIN)/NX
   allocate(x_vec(NX+1))
   do i=1,NX+1
   x_vec(i)=(i-1)*DELTAX
   end do

 ! define i for the source
   xsource = 0.1_dp*XMAX    ! arbitrary source position in [m]
   mindist = XMAX
   do i = 2, NX
   dist = xsource - DELTAX*i
      if (abs(dist) < mindist) then
         mindist = dist
         ISOURCE = i
      end if
   end do
 
 
  !Number of iterations between each display
 IT_DISPLAY = 50

 ! velocity models (cp, cs)
   cp=(/((vmin+k*a), k=0,NX)/)
   cs=cp/1.732_dp
  
 ! rho (density)
   rho=(/((2400_dp+k*b), k=0,NX)/)
 

 ! mu (stiffness/rigidity matrix) 
   do i=1,NX+1
       mu(i)=rho(i)*cs(i)*cs(i)
   end do

 
 ! stability 
    courant_number1=2500_dp*DELTAT*DELTAX
   write(*,*)'courant number1 = ', courant_number1
    courant_number2 = 3000_dp*DELTAT*DELTAX
   write(*,*)'courant number2 = ', courant_number2
   if((courant_number1>1_dp) .or.( courant_number2>1_dp)) then
      write(*,*) ' time step is too large, simulation will be unstable'
   end if
 
 !
   dx2=DELTAX*DELTAX
   dt2=DELTAT*DELTAT


 ! predictor corrector data vectors
  allocate(u_field(1:NX+1, NSTEP)) ! save the field
  !
  allocate(sig(NX+1+2*npml))
  !
  npml=10
  !
  allocate(cpe(NX+1+2*npml))  
  !
  vm=sum(cp)/size(cp,1)
  !
  call pml (DELTAX, vm, npml, sig, NX) 
  !
  call extend (cp, npml, cpe) 
 
 write(*,*) ISOURCE, DELTAT, DELTAX
!stop
  
  ! Forward modelling O(2,2) and Conv(2,2)
 if(((OPT2).eqv.(.true.)).or.((CONV2).eqv.(.true.))) then

  call optimally_accurate_22 (UI,u,u1,u_field,rho,mu,NSTEP,NX,dx2,dt2,DELTAT,&
        ISOURCE,A_conv,A_opt,K_conv, K_opt, dA, dK, dobs, CONV2, OPT2)
    
       if ((CONV2).eqv.(.true.)) then
        write(*,*) 'Conv2 will be computed'
       end if

       if ((OPT2).eqv.(.true.))  then
        write(*,*) 'OPT2 will be computed'
       end if
   
   
 end if
 !------------------------------------------------------------------------------

 
 ! Forward modelling O(2,4) and Conv(2,4)
  if(((OPT4).eqv.(.true.)).or.((CONV4).eqv.(.true.))) then
  
  call optimally_accurate_24 (UI4,u,u1,u_field,rho,mu,NSTEP,NX,dx2,dt2,DELTAT,&
        ISOURCE,A_conv4,A_opt4,K_conv4, K_opt4, dA4, dK4, dobs, CONV4, OPT4)
     
      if ((CONV4).eqv.(.true.))   then
        write(*,*) 'Conv4 will be computed'
       end if

      if ((OPT4).eqv.(.true.))   then
        write(*,*) 'OPT4 will be computed'
      end if

  end if
  !------------------------------------------------------------------------------
  

  ! Write data 
  call write_data(u_field, dobs, NSTEP) ! ASCII 
  !call write_data_bin(u_field, dobs, NSTEP) ! BINARY
  
 !write(*,*) mu
  
end program main


