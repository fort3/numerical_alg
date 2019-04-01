program main
	use types
	use newcot
	real(kind=rkind), dimension(:,:), allocatable:: rMatrix
	real(kind=rkind)																													::	inte,val
	integer(kind=ikind)																										::	nsteps,n
	
	print*,"This is the Trapezoid Approximation: "	
	call trap(inte,0.0_rkind,20.0_rkind,fnc,n)
	
	print*,"This is the Romberg Algorithm Matrix: "
 call rom(rMatrix,0.0_rkind,20.0_rkind,1.0e-10_rkind,fnc)
 
 print*,"This is the Newton Cotes Formular: "
 call lapoly(val,0.0_rkind,20.0_rkind,fnc,n)
 
end program
