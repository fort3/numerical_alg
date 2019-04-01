module newcot

	contains
	function fnc(x) result(y)
	use types
		real(kind=rkind), intent(in) :: x
		real(kind=rkind) 												:: y
				
		y = (x)/(sin(x) + 2)
		
	end function fnc
		
	subroutine trap(integ,a,b,fnc,nsteps)
	use types
		integer(kind=ikind) 												:: i,j
		integer(kind=ikind)													:: nsteps
		real(kind=rkind), intent(in) 			:: a, b
		real(kind=rkind), intent(out) 		:: integ
		real(kind=rkind) 															:: h,r4,x1,x2,dx,sum_
		integer																									::	datafile
		
	interface
	function fnc(x) result(y)
		use types
		real(kind=rkind), intent(in)		:: x
		real(kind=rkind)														:: y
	end function
	end interface	
	
	open(newunit=datafile, file="trap_output.mtx")
	
	nsteps = 1
	r4 = 0.5
	
	do
		integ = 0.0
		dx = (b-a)/nsteps
		
		do j = 1,nsteps
		x1= (j-1)*dx
		x2= j*dx
		if(x1==r4) then
			integ = 0.0
		else
			sum_ = (fnc(x1) + fnc(x2)) /2.0 * dx
		end if
			integ = integ + sum_
		end do
		print *, "Nsteps: ", nsteps, "Integ: ", integ, "Error: ", abs((1.0e-10_rkind)-integ)
		
		write(unit=datafile, fmt=*) nsteps, integ, abs((1.0e-10_rkind)-integ)
		write(*,*)
		if(nsteps >= 1000) then
			print*, "End"
			exit
		else
			nsteps=nsteps+1
		end if
		
	end do
	close(datafile)
 end subroutine trap
	
	subroutine rom (romMatrix,a,b,errf,fnc)
	use types
		real(kind=rkind), intent(in) 																	:: errf
		integer(kind=ikind)												 														:: nsteps,rombergdim,itmax=9
		real(kind=rkind)																														:: a,b,s,va
		integer(kind=ikind) 																										:: i,j,k,m
		real(kind=rkind), dimension(:,:), allocatable :: romMatrix
		integer																																							::	datfile
	interface
	function fnc(x) result(val)
		use types
		real(kind=rkind), intent(in) :: x
		real(kind=rkind) 												:: val
	end function
	end interface
	
	rombergdim = 10
	
	allocate(romMatrix(1:rombergdim,0:rombergdim-1))
	
	do i=lbound(romMatrix,1), ubound(romMatrix,1)
		romMatrix(i,0)=(b-a)/2.0_rkind*(fnc(a)+fnc(b))
	end do
	
	do i=1,itmax
		s=0.0_rkind
		do j=1,2**(i-1)
			s=s+fnc(a+(2.0_rkind*j-1.0_rkind)*(b-a)/2.0_rkind**i)
		end do
		romMatrix(1,i)=1.0_rkind/2*romMatrix(1,i-1)+(b-a)/2**i*s
		do m=1,i
			do k=i,1,-1
				va=4.0_rkind**m*romMatrix(m,k)-romMatrix(m,k-1)
				romMatrix(m,k-1)=va/(4.0_rkind**m-1.0_rkind)
				write(*,'(f11.5,t3)', advance="no"), romMatrix(m,k)
			end do
			write(*,*)
		end do
		if(abs(romMatrix(m,0)-romMatrix(m-1,0))< errf) then
			exit
		end if
	end do
	end subroutine rom
	
	subroutine lapoly(val,a,b,fnc,nst)
		use types
		integer(kind=ikind) 												:: i,j
		integer(kind=ikind)													:: nst
		real(kind=rkind), intent(in) 			:: a, b
		real(kind=rkind), intent(out) 		:: val
		real(kind=rkind) 															:: h,r4,x1,x2,x3,x4,dx,sum_
		integer																									::	dfile
		
		interface
		function fnc(a) result(d)
		use types
		real(kind=rkind), intent(in)				:: a
		real(kind=rkind)																:: d
		end function
		end interface
		
		open(newunit=dfile, file="newtonCotes_output.mtx")
	
		nst = 1
		r4 = 0.5
	
		do
		val = 0.0
		dx = (b-a)/nst
			do j = 1,nst
			x1= (j-3)*dx
			x2= (j-1)*dx
			x3= (j-1)*dx
			x4= (j)*dx
			if(x1==r4) then
				val = 0.0
			else
				sum_ = (fnc(x1) + 3*fnc(x2) + 3*fnc(x3) + fnc(x4)) /8.0 * (dx*3)
			end if
				val = val + sum_
			end do
		print *, "Nsteps: ", nst, "Integ: ", val, "Error: ", abs((1.0e-10_rkind)-val)
		
		write(unit=dfile, fmt=*) nst, val, abs((1.0e-10_rkind)-val)
		write(*,*)
		if(nst >= 10) then
			print*, "End"
			exit
		else
			nst=nst+1
		end if
		
	end do
	
	close(dfile)
	end subroutine lapoly
	
end module newcot



