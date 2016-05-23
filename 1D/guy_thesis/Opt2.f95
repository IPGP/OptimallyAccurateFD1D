
module mod_FW
!
contains
!------------------------------------------------------------------------
!
subroutine optimally_accurate_22 (UI,u,u1,u_field,rho,mu,NSTEP,NX,&
 dx2,dt2,DELTAT,ISOURCE,A_conv,A_opt,K_conv, K_opt, dA, dK, dobs, &
 CONV2, OPT2)
  !
  implicit none
  !
  ! external variables
  integer, parameter :: dp=kind(1.0d0)
  real(dp) :: DELTAT
  integer :: NSTEP, ISOURCE, NX
  real(dp), dimension(:) :: u, u1, rho, mu
  real(dp), dimension(:,:) :: u_field, dobs
  real(dp), dimension(:,:,:) :: A_conv, A_opt, K_conv, K_opt, &
  & dA, dK
  logical :: CONV2, OPT2
  !  
  
  !internal variables
  integer :: i, it  
  real(dp), parameter ::pi=3.141592653589793238462643d0
  real(dp) ::umaxv, to, fo, a, ampl, dx2, dt2, t, source_term, &
  & value_du_dxx
  real(dp), dimension(1:NX+1) :: unm1,unm2, u1nm1, u1nm2, utotal, g
  real(dp), dimension(:,:,:) ::  UI
  !


  ! operators as well as their corresponding residuals
  do i=2, NX-1
     if ((CONV2 .eqv. (.true.)).or.((OPT2).eqv.(.true.))) then

     A_conv(:,:,i)=(reshape((/real(dp)::0,0,0 ,rho(i), -2_dp*rho(i), rho(i), 0,0,0/), (/3,3/)))*(1_dp/dt2)

     k_conv(:,:,i)=(reshape((/real(dp)::0,(mu(i-1)+mu(i)),0,0,-(mu(i-1)+2_dp*mu(i)+mu(i+1)), 0, 0,(mu(i)+mu(i+1)),0/), &
                   (/3,3/)))*1_dp/2_dp/dx2
     end if
     
   
  if ((OPT2).eqv.(.true.)) then

     A_opt(:,:,i)=(reshape((/real(dp)::rho(i), -2_dp*rho(i), rho(i),10_dp*rho(i),-20_dp*rho(i),10_dp*rho(i),rho(i),&
                  -2_dp*rho(i), rho(i)/), (/3,3/)))*1_dp/12_dp/dt2

     K_opt(:,:,i)=(reshape((/real(dp)::(mu(i-1)+mu(i)),10_dp*(mu(i-1)+mu(i)), (mu(i-1)+mu(i)),-(mu(i-1)+2_dp*mu(i)+mu(i+1)),&
                  -10_dp*(mu(i-1)+2_dp*mu(i)+mu(i+1)),-(mu(i-1)+2_dp*mu(i)+mu(i+1)), (mu(i)+mu(i+1)), 10_dp*(mu(i)+mu(i+1)),&
                   (mu(i)+mu(i+1))/), (/3,3/)))*1_dp/24_dp/dx2

     dA(:,:,i)=A_opt(:,:,i)-A_conv(:,:,i)

     dK(:,:,i)=K_opt(:,:,i)-K_conv(:,:,i)
   end if

  end do 
  !

  !parameters for the source
  fo=100._dp    ! frequency
  to=1.20_dp/fo ! time of source activation
  ampl=1.0d9
  a=pi*pi*fo*fo ! amplitude of the signal
  umaxv=(ampl*exp(-a*to**2))/2._dp*DELTAT

!--------------------------------------------------------------------------------------
!                           MODELLING CODE
!--------------------------------------------------------------------------------------

  ! predictor corrector scheme
  do it = 1, NSTEP
     u(:) = 0._dp
     
     
   if ((CONV2 .eqv. (.true.)).or.((OPT2).eqv.(.true.))) then
       ! ricker source time function (second derivative of a Gaussian)
       t = (it-1)*DELTAT
       source_term = ampl*exp(-a*(t-to)**2)/rho(ISOURCE)

    ! predictor step
      !g(:)=0._dp
      !g(ISOURCE)=(source_term*DELTAT*DELTAT)/rho(ISOURCE)

   ! wavefield extrapolation with conventional operators
      do i=2,NX-1
       value_du_dxx = unm1(i-1)*K_conv(2,1,i)+unm1(i)*K_conv(2,2,i)+unm1(i+1)*K_conv(2,3,i)
       u(i) = (value_du_dxx-unm1(i)*A_conv(2,2,i)-unm2(i)*A_conv(3,2,i))/A_conv(1,2,i)
      end do
      
   if (CONV2 .eqv. (.true.))  then
     u(ISOURCE) = source_term*DELTAT
   end if
   end if
   !


    if ((OPT2).eqv.(.true.)) then
   ! corrector step
      u1(:) = 0._dp
     UI(:,:,:) = 0._dp
     g(:)=0._dp

     do i=2,NX-1
        UI(1:3,1:3,i-1)= reshape((/real(kind=dp):: u(i-1),unm1(i-1),unm2(i-1),u(i),unm1(i),&
            unm2(i),u(i+1),unm1(i+1),unm2(i+1)/),(/3,3/))
       
        g(i)=-sum(sum((dA(:,:,i-1)-dK(:,:,i-1))*UI(:,:,i-1),dim=1),dim=1)/A_conv(1,2,i)
    end do

    !
    do i=2,NX-1
      !value_du_dxx= (u1nm1(i-1)-2*u1nm1(i)+u1nm1(i+1))/dx2
      !u1(i)= 2._dp*u1nm1(i)-u1nm2(i)+mu(i)*value_du_dxx*(DELTAT**2._dp)/rho(i)+g(i)
      u(i)= u(i) + (g(i))/A_conv(1,2,i)
    end do
   u(ISOURCE) = source_term*DELTAT
    end if
   !

    if ((OPT2).eqv.(.true.)) then
   !corrected displacement field
   !utotal = u+u1
   
   !update of u_conv
   unm2 = unm1
   unm1 = u

   !update of du
   !u1nm2 = u1nm1
   !u1nm1 = u1
   dobs(:, it) = u 
   end if 

   if ((CONV2).eqv.(.true.)) then
   unm2 = unm1 ! update of time t-1 to t
   unm1 = u  ! update of time t 
   u_field(:,it)=u
   end if

   ! Dirichlet BC
   !u(1)=0
   !u(NX+1)=0
 
   end do

  end subroutine Optimally_accurate_22
 !
 
 
subroutine optimally_accurate_24 (UI4,u,u1,u_field,rho,mu,NSTEP,NX,&
 dx2,dt2,DELTAT,ISOURCE,A_conv4,A_opt4,K_conv4, K_opt4, dA4, dK4, dobs, CONV4, OPT4)
  !
  implicit none
  !
  ! external variables
  integer, parameter :: dp=kind(1.0d0)
  real(dp) :: DELTAT
  integer :: NSTEP, ISOURCE, NX
  real(dp), dimension(:) :: u, u1, rho, mu
  real(dp), dimension(:,:) :: u_field, dobs
  real(dp), dimension(:,:,:) :: A_CONV4, A_OPT4, K_CONV4, K_OPT4, dA4, dK4
  logical :: CONV4, OPT4
  !  
  
  !internal variables
  integer :: i, it  
  real(dp), parameter ::pi=3.141592653589793238462643d0
  real(dp) ::umaxv, to, fo, a, ampl, dx2, dt2, t, source_term, value_du_dxx
  real(dp), dimension(1:NX+1) :: unm1,unm2, u1nm1, u1nm2, utotal, g
  real(dp), dimension(:,:,:) ::  UI4
  !


  ! operators as well as their corresponding residuals
 
  
  do i = 3, NX-1

  if ((CONV4 .eqv. (.true.)).or.((OPT4).eqv.(.true.))) then
  
 A_conv4(:,:,i)=(reshape((/real(dp):: 0,0,0 ,0,0,0, rho(i), -2._dp*rho(i), rho(i), 0,0,0, 0,0,0/), (/3,5/)))*(1._dp/dt2)

     k_conv4(:,:,i)=(reshape((/real(dp):: 0, -(mu(i-2)+mu(i)), 0, 0, 16._dp*(mu(i)+mu(i-1)), 0, 0, &
                                   -16*(mu(i-1)+2*mu(i)+mu(i+1)) +(mu(i-2)+2*mu(i)+mu(i+2)),&
                     0, 0, 16*(mu(i)+mu(i+1)), 0, 0,-(mu(i)+mu(i+2)),0/), (/3,5/)))*(1._dp/(24._dp*dx2))
  
  end if


   if ((OPT4).eqv.(.true.)) then

     A_opt4(:,:,i)= (reshape((/real(dp):: -rho(i), 2._dp*rho(i), -rho(i), 4._dp*rho(i), -8._dp*rho(i), 4._dp*rho(i), 84._dp*rho(i),&
                  -168._dp*rho(i), 84._dp*rho(i), 4._dp*rho(i), -8._dp*rho(i), 4._dp*rho(i), -rho(i), 2._dp*rho(i), -rho(i)/), &
                   (/3,5/)))/(90._dp*dt2)

     K_opt4(:,:,i)=(reshape((/real(dp):: -(mu(i-2)+mu(i)), -10._dp*(mu(i-2)+mu(i)), -(mu(i-2)+mu(i)), 16._dp*(mu(i)+mu(i-1)), &
              16._dp*10._dp*(mu(i)+mu(i-1)), 16._dp*(mu(i)+mu(i-1)), -16._dp*(mu(i-1)+2*mu(i)+mu(i+1))+(mu(i-2)+2*mu(i)+mu(i+2)), &
         -16*10*(mu(i-1)+2*mu(i)+mu(i+1))+(mu(i-2)+2*mu(i)+mu(i+2)), -16*(mu(i-1)+2*mu(i)+mu(i+1))+10*(mu(i-2)+2*mu(i)+mu(i+2)), &
         16*(mu(i)+mu(i+1)), 16*10*(mu(i)+mu(i+1)), 16*(mu(i)+mu(i+1)), -(mu(i)+mu(i+2)), -10*(mu(i)+mu(i+2)), -(mu(i)+mu(i+2))/),&
              (/3,5/)))*(1/(288*dx2))

     dA4(:,:,i)=A_opt4(:,:,i)-A_conv4(:,:,i)

     dK4(:,:,i)=K_opt4(:,:,i)-K_conv4(:,:,i)
   end if

  end do 
  !

  !parameters for the source
  fo=100._dp    ! frequency
  to=1.20_dp/fo ! time of source activation
  ampl=1.0d9
  a=pi*pi*fo*fo ! amplitude of the signal
  umaxv=(ampl*exp(-a*to**2))/2._dp*DELTAT

!--------------------------------------------------------------------------------------
!                           MODELLING CODE
!--------------------------------------------------------------------------------------

  ! predictor corrector scheme
  do it = 1, NSTEP
     u(:) = 0._dp
   
     
  if ((CONV4 .eqv. (.true.)).or.((OPT4).eqv.(.true.))) then
    ! ricker source time function (second derivative of a Gaussian)
    t = (it-1)*DELTAT
    source_term = ampl*exp(-a*(t-to)**2)/rho(ISOURCE)

    ! predictor step
    !g(:)=0._dp
    !g(ISOURCE)=(source_term*DELTAT**2)/rho(ISOURCE)

   ! wavefield extrapolation with conventional operators
   do i=3,NX-2
    value_du_dxx = unm1(i-2)*K_conv4(2,1,i)+unm1(i-1)*K_conv4(2,2,i)+unm1(i)*K_conv4(2,3,i) &
      +unm1(i+1)*K_conv4(2,4,i)+unm1(i+2)*K_conv4(2,5,i)
       u(i) = (value_du_dxx-unm1(i)*A_conv4(2,3,i)-unm2(i)*A_conv4(3,3,i))/A_conv4(1,3,i)
   end do

    if (CONV4 .eqv. (.true.)) then
    u(ISOURCE)=source_term*DELTAT
    end if
 end if
   !

 if ((OPT4).eqv.(.true.)) then
   ! corrector step
     g(:)=0._dp
     UI4(:,:,:) = 0._dp
     u1(:) = 0._dp

   do i = 3,NX-1
        UI4(1:3,1:5,i-2)= reshape((/real(kind=dp):: u(i-2), unm1(i-2), unm2(i-2), u(i-1), unm1(i-1), unm2(i-1), &
          u(i), unm1(i), unm2(i), u(i+1), unm1(i+1), unm2(i+1), u(i+2), unm1(i+2), unm2(i+2)/),(/3,5/))
       
        g(i)=-sum(sum((dA4(:,:,i-1)-dK4(:,:,i-1))*UI4(:,:,i-2),dim=1),dim=1)*dt2/rho(i)
   end do
   !

   !
   do i = 3, NX-1
    u(i) = u(i)+g(i)/A_conv4(1,3,i)
   end do
  
  u(ISOURCE) = source_term*DELTAT
 end if


    if ((OPT4).eqv.(.true.)) then
   !corrected displacement field
  ! utotal = u+u1

   !update of u_conv
   unm2 = unm1
   unm1 = u

   !update of du
   !u1nm2 = u1nm1
   !u1nm1 = u1

   !u=utotal
   dobs(:, it) = u 
   end if 

   if ((CONV4).eqv.(.true.)) then
   unm2 = unm1 ! update of time t-1 to t
   unm1 = u  ! update of time t 
   u_field(:,it)=u
   end if

   ! Dirichlet BC
   !u(1)=0
   !u(NX+1)=0
   
   end do

 end subroutine Optimally_accurate_24

end module mod_FW
!--------------------------------------------------------------------------------------------------------------------------------------------


