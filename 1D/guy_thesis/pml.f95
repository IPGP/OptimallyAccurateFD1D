module mod_pml
 !
 contains
 !
   subroutine pml(dx, vm, npml, sig, nz)
    !
    implicit none
    !
    !variables
    integer, parameter :: dp=kind(1.0d0)
    integer :: nz, npml, i
    real(dp), dimension(:) :: sig
    real(dp) :: vm, L, R, vpml, val, dx
    !
    L = real (npml-1)*dx
    R = 5.0_dp
    vpml =3.0_dp*vm/2._dp/L*R
    !
    sig = 0._dp
    ! 
    do i = 1, npml
       val = vpml * (real(i)/real(npml))**2
       sig(npml-i+1)=val
       sig(nz+npml+i)=val 
    end do
  end subroutine pml
  
  
  subroutine extend(cp, npml, cpe)
   !
   implicit none
   !
   integer, parameter :: dp=kind(1.0d0)
   integer :: npml, iz, NX
   !
   real(dp), dimension(:) :: cp, cpe
   !
   NX = size(cp,1)   
   !
   !----------- extended velocity model---------------------
   ! Inside the model
   cpe = 0._dp 
   !
   do iz = 1, NX
      cpe(iz+npml)=cp(iz)
   end do
   !
   ! Edges of the model
   !
   do iz = 1, npml
      cpe(iz) = cp(1)
      cpe(iz+NX+npml) = cp(NX)
   end do
   !
  end subroutine extend
end module mod_pml  
