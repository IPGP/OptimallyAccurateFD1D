	program normal2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computation of the synthetic seismograms in the time domain
c using the normal operators.
c
c						1995.10  N.Takeuchi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz
	real*8 pi
	parameter ( maxnz = 10000 )
	parameter ( pi=3.1415926535897932d0 )
c
c parameters for the gridding
	integer nt,nz,it,ist,isz
	real*8 dt,dz
c parameter for the structure
	real*8 rrho(4),kkappa(4)
	real*8 rho(maxnz+1),kappa(maxnz+1)
c parameter for the wavefield
	real*8 t
	real*8 u(maxnz+1),u1(maxnz+1),u2(maxnz+1),f(maxnz+1)
	real*8 e1(maxnz+1),e2(maxnz+1),e3(maxnz+1)
	real*8 f1(maxnz+1),f2(maxnz+1),f3(maxnz+1)
c parameter for the source
	real*8 tp,ts
c parameter for the receiver
	integer nr
c
c reading the parameter files
	call pinput( maxnz,nt,nz,dt,dz,rrho,kkappa,tp,ts,nr )
c Initializing the data
	call datainit( maxnz,u )
	call datainit( maxnz,u1 )
	call datainit( maxnz,u2 )
	call datainit( maxnz,f )
	call datainit( maxnz,rho )
	call datainit( maxnz,kappa )
c computing the intermediate parameters
	call calstruct( rrho,dz,nz,rho )
	call calstruct( kkappa,dz,nz,kappa )
	call cales( nz,rho,kappa,dt,dz,e1,e2,e3,f1,f2,f3 )
	ist = dnint( 2 * tp / dt )
	isz = nz / 2 + 1
c
	t = 0.d0
	write(6,*) real(t),real(u(nr))
	do 100 it=0,nt
	  call calf( it,t,ist,isz,dt,dz,rho(isz),tp,ts,f )
c evaluating the next step
	  call calstep( nz,e1,e2,e3,f1,f2,f3,u,u1,u2,f )
c increment of t
	  t = t + dt
	  if ( mod(it,10).eq.9 ) write(6,*) real(t),real(u(nr))
c	  write(6,*) real(t),real(u(nr))
  100	continue
c
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pinput( maxnz,nt,nz,dt,dz,rho,kappa,tp,ts,nr )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz,nt,nz,nr
	real*8 dt,dz,rho(4),kappa(4),tp,ts
	character*80 tmpfile,dummy
c
	data tmpfile / '/tmp/work' /
c
c temporary file open
	open( unit=11, file=tmpfile, status='unknown' )
c writing to the temporary file
  100	continue
	  read(5,110) dummy
  110	  format(a80)
	  if ( dummy(1:1).eq.'c' ) goto 100
	  if ( dummy(1:3).eq.'end' ) goto 120
	  write(11,*) dummy
	  goto 100
  120	continue
c temporary file close
	close(11)
c 
c temporary file open
	open( unit=11, file=tmpfile, status='unknown' )
c reading the parameter
	read(11,*) nt,nz
	if ( nz.gt.maxnz ) pause 'nz is too large (pinput).'
	read(11,*) dt,dz
	read(11,*) rho(1),rho(2),rho(3),rho(4)
	read(11,*) kappa(1),kappa(2),kappa(3),kappa(4)
	read(11,*) tp,ts
	read(11,*) nr
c temporary file close
	close(11)
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine datainit( maxnz,u )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz
	real*8 u(*)
	integer i
c
	do 100 i=1,maxnz+1
	  u(i) = 0.d0
  100	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calstruct( rrho,dz,nz,rho )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nz
	real*8 rrho(4),dz,rho(*)
	integer i,j
	real*8 r,rmax,coef,trho
c
	rmax = dz * nz
	do 100 i=1,nz+1
	  r = dble(i-1) * dz
          trho = 0.d0
	  do 110 j=1,4
	    if ( j.eq.1 ) then
	      coef = 1.d0
	    else
	      coef = coef * ( r / rmax )
	    endif
	    trho = trho + rrho(j) * coef
  110	  continue
          rho(i) = trho
  100	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cales( nz,rho,kappa,dt,dz,e1,e2,e3,f1,f2,f3 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nz
	real*8 rho(*),kappa(*),dt,dz,e1(*),e2(*),e3(*),f1(*),f2(*),f3(*)
	integer iz
	real*8 dt2,dz2
c
	dt2 = dt * dt
	dz2 = dz * dz
c
	e1(1) = 0.d0
	e2(1) = 2.d0 - ( kappa(1) + kappa(2) )
     &	               / rho(1) * dt2 / dz2
	e3(1) = ( kappa(1) + kappa(2) ) / rho(1) * dt2 / dz2
	f1(1) = 0.d0
	f2(1) = e2(1) / 12.d0
	f3(1) = ( e3(1) - 2.d0 ) / 12.d0
	do 100 iz=2,nz
	  e1(iz) = ( kappa(iz-1) + kappa(iz) )
     &	           / ( 2.d0 * rho(iz) ) * dt2 / dz2
	  e2(iz) = 2.d0
     &	           - ( kappa(iz-1) + 2.d0 * kappa(iz) + kappa(iz+1) )
     &	             / ( 2.d0 * rho(iz) ) * dt2 / dz2
	  e3(iz) = ( kappa(iz) + kappa(iz+1) )
     &	           / ( 2.d0 * rho(iz) ) * dt2 / dz2
	  f1(iz) = ( e1(iz) - 1.d0 ) / 12.d0
	  f2(iz) = e2(iz) / 12.d0
	  f3(iz) = ( e3(iz) - 1.d0 ) / 12.d0
  100	continue
	e1(nz+1) = ( kappa(nz) + kappa(nz+1) )
     &	           / rho(nz+1) * dt2 / dz2
	e2(nz+1) = 2.d0 - ( kappa(nz) + kappa(nz+1) )
     &	                  / rho(nz+1) * dt2 / dz2
	e3(nz+1) = 0.d0
	f1(nz+1) = ( e1(nz+1) - 2.d0 ) / 12.d0
	f2(nz+1) = e2(nz+1) / 12.d0
	f3(nz+1) = 0.d0
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calf( it,t,ist,isz,dt,dz,rho,tp,ts,f )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
	integer it,ist,isz
	real*8 t,dt,dz,rho,tp,ts,f(*)
	real*8 b
c
	if ( it.le.ist ) then
	  b = pi * ( t - ts ) / tp
	  f(isz) = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b) / dz
	  if ( (it.eq.0).or.(it.eq.ist) ) f(isz) = f(isz) / 2.d0
	  f(isz) = f(isz) * dt * dt / rho
	else
	  f(isz) = 0.d0
	endif
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calstep( nz,e1,e2,e3,f1,f2,f3,u,u1,u2,f )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nz
	real*8 e1(*),e2(*),e3(*),f1(*),f2(*),f3(*)
	real*8 u(*),u1(*),u2(*),f(*)
	integer iz
	real*8 tmp1,tmp2,tmp3
c
c evaluating the u using the unmodified operators
	u(1) = - u2(1)
     &	       + e2(1) * u1(1)
     &	       + e3(1) * u1(2)
     &	       + f(1)
	do 100 iz=2,nz
	  u(iz) = - u2(iz)
     &	          + e1(iz) * u1(iz-1)
     &	          + e2(iz) * u1(iz)
     &	          + e3(iz) * u1(iz+1)
     &	          + f(iz)
  100	continue
	u(nz+1) = - u2(nz+1)
     &	          + e1(nz+1) * u1(nz)
     &	          + e2(nz+1) * u1(nz+1)
     &	          + f(nz+1)
c computing du using 1-D Born approximation
	tmp2 = ( u(1) - u1(1) - u1(1) + u2(1) )
	do 110 iz=1,nz
	  tmp3 = ( u(iz+1) - u1(iz+1) - u1(iz+1) + u2(iz+1) )
	  u(iz) = u(iz)
     &	          + f1(iz) * tmp1
     &	          + f2(iz) * tmp2
     &	          + f3(iz) * tmp3
	  tmp1 = tmp2
	  tmp2 = tmp3
  110	continue
	u(nz+1) = u(nz+1)
     &	          + f1(nz+1) * tmp1
     &	          + f2(nz+1) * tmp2
c swapping u1 & u2
	do 120 iz=1,nz+1
	  u2(iz) = u1(iz)
	  u1(iz) = u(iz)
  120	continue
c
	return
	end
