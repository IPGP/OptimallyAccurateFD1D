!
module mod_data
!
contains
!-----------------------------------------------------------------------
!
subroutine write_data (u_field, dobs, NSTEP)
  !
  implicit none
  !
  ! external variables
  integer, parameter :: dp=kind(1.0d0)
  integer :: NX, NSTEP, i, j
  real(dp), dimension(:,:) :: u_field, dobs
  !
  NX=size(u_field, 1)  

  !
  open(unit=1, file='wavefield.txt', status='unknown', &
       form='formatted', access='sequential', action='write')
  !
     do j = 1,NSTEP
       do i = 1,NX+1
            write (unit=1, fmt='(E30.20,1x)') (u_field(i,j))
       end do
     end do
  !
  close(unit=1)
  !  
  !
  !
  open(unit=10, file='dobs.txt', status='unknown', &
       form='formatted', access='sequential', action='write')
  !
     do j = 1,NSTEP
       do i = 1,NX+1
            write (unit=10, fmt='(E30.20,1x)') (dobs(i,j))
       end do
     end do
  !
  close(unit=10) 

  !
end subroutine write_data


subroutine write_data_bin (u_field, dobs, NSTEP)
  !
  implicit none
  !
  ! external variables
  integer, parameter :: dp=kind(1.0d0)
  integer :: NX, NSTEP, i, j
  real(dp), dimension(:,:) :: u_field, dobs
  !
  NX = size(u_field, 1)
  

  ! internal variables
  !

  open(unit=12, file='wavefield.dat', status='unknown', &
       form='unformatted', access='direct', recl=4*NX*NSTEP)
  !
     do j = 1,NSTEP-30000
            write (unit=12, rec=1) (u_field(i,j), i=1, NX)
     end do
  !
  close(unit=12)
  !  
 
  !
  open(unit=11, file='dobs.dat', status='unknown', &
       form='unformatted', access='direct', recl=4*NX*NSTEP)
  !
     do j = 1,NSTEP-30000
            write (unit=11, rec=1) (dobs(i,j), i=1, NX)
       end do
  !
  close(unit=11) 

  !

end subroutine write_data_bin

end module mod_data
