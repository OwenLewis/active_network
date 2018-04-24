! ============================================================================
! Name        : im_bnd_opers_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing all subroutines pertaining to immersed
!				boundary functionality (spreading, interpolating, etc)
! ============================================================================
module im_bnd_opers_mod
USE global_param_mod
USE useful_functions_mod
implicit none
PRIVATE
PUBLIC :: grid_spread, boundary_interpolate, boundary_time_integrate

contains

!========================================================================================================
!THIS IS THE MAIN ROUTINE TO SPREAD FORCES ONTO THE EULERIAN GRID

	subroutine grid_spread(loc_bnd,force_bnd,force_horiz_grid,force_vert_grid)

	!FIRST TWO INPUTS SHOULD BE ARRAYS DEFINING LOCATION OF LAGRANGIAN POINTS
	!AND FORCES ACTING AT THOS LAGRANGIAN POINTS
	REAL, DIMENSION(1:Mb,1:2) :: loc_bnd, force_bnd
	!LAST TWO INPUTS ARE ARRAYS WHERE THE SPREAD FORCES WILL BE WRITTEN
	!NOTICE THAT THEY ARE SAME SIZE BECAUSE OF THE CENTERED GRID STRUCTURE
	REAL, DIMENSION(1:Nx,1:Ny) :: force_horiz_grid, force_vert_grid
	!VARIABLES FOR USE IN THE PROGRAM
	INTEGER :: I, left, bottom, innerI, innerK, xin, yin
	INTEGER, DIMENSION(1:4) :: xindeces, yindeces, modxin, modyin
	REAL, DIMENSION(1:4) :: coords, distance, xweights, yweights

	! begin by zeroing out the array we will spread the force into
	! POSSIBLY CHANGE THIS PIECE OF CODE, I'M NOT SURE
	force_horiz_grid = 0.0
	force_vert_grid = 0.0

	!loop through the length of the IB array
	do I = 1,Mb
	!FIRST WE WILL SPREAD THE HORIZONTAL FORCES TO THE HORIZONTAL GRID
		!Determine the index corresponding to X-coordinate immediately below the IB point
		left = floor(loc_bnd(I,1)/hx + 1./2.)
		!Form indeces defining  points to spread to
		xindeces = (/left-1,left,left+1,left+2/)
        !Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modxin = modulo(xindeces-1,Nx)+1
		!Calculate the x-coordinate of those points
		coords = hx*(xindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - loc_bnd(I,1)

		!Calculate weights via the approximate delta function
		xweights = better_phi(distance/hx)/hx



		!Determine the index corresponding to Y-coordinate immediately below the IB point
		bottom = floor(loc_bnd(I,2)/hy + 1./2.)
		!Form indeces defining  points to spread to
		yindeces = (/bottom-1, bottom, bottom+1, bottom+2/)
        !Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modyin = modulo(yindeces-1,Ny)+1
		!Calculate the y-coordinate of those points
		coords = hy*(yindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - loc_bnd(I,2)

		!Calculate weights via the approximate delta function
		yweights = better_phi(distance/hy)/hy


		!Now we loop through the appropiate indeces and spread the weighted values
		do innerI = 1,4
			do innerK = 1,4
				!Add the force in the X-direction times the weights from X and Y delta functions
				!and scale by the lagrangian arc length difference
				force_horiz_grid(modxin(innerI),modyin(innerK)) = force_horiz_grid(modxin(innerI),modyin(innerK)) &
				                            &+ force_bnd(I,1)*xweights(innerI)*yweights(innerK)*ref_bnd%ds0(I)
                force_vert_grid(modxin(innerI),modyin(innerK)) = force_vert_grid(modxin(innerI),modyin(innerK)) &
                                            &+ force_bnd(I,2)*xweights(innerI)*yweights(innerK)*ref_bnd%ds0(I)

			end do
		end do


!	!NOW WE WILL SPREAD THE VERTICAL FORCES TO THE VERTICAL GRID
!		!Determine the index corresponding to Y-coordinate immediately below the IB point
!		bottom = floor(loc_bnd(I,2)/hy + 1./2.)
!		!Form indeces defining  points to spread to
!		yindeces = (/bottom-1,bottom,bottom+1,bottom+2/)
!		!Calculate the y-coordinate of those points
!		coords = hy*(yindeces - 1./2.)
!		!Calculate distance from these Eulerian points to IM point
!		distance = coords - loc_bnd(I,2)
!
!		!Calculate weights via the approximate delta function
!		yweights = better_phi(distance/hy)/hy
!
!
!
!		!Determine the index corresponding to X-coordinate immediately below the IB point
!		left = floor(loc_bnd(I,1)/hx + 1./2.)
!		!Form indeces defining  points to spread to
!		xindeces = (/left-1, left, left+1, left+2/)
!		!Calculate the X-coordinate of those points
!		coords = hx*(xindeces - 1./2.)
!		!Calculate distance from these Eulerian points to IM point
!		distance = coords - loc_bnd(I,1)
!
!		!Calculate weights via the approximate delta function
!		xweights = better_phi(distance/hx)/hx
!
!
!		!Now we loop through the appropiate indeces and spread the weighted values
!		do innerI = 1,4
!			do innerK = 1,4
!				!Determine the indeces in the Eulerian grid that we will spread to.
!				xin = mymod(xindeces(innerI),Nx)
!				yin = mymod(yindeces(innerK),Ny)
!				!Add the force in the Y-direction times the weights from X and Y delta functions
!				!and scale by the lagrangian arc length difference
!				force_vert_grid(xin,yin) = force_vert_grid(xin,yin) + force_bnd(I,2)*xweights(innerI)*yweights(innerK)*ref_bnd%ds0(I)
!
!			end do
!		end do

	end do
	end subroutine grid_spread

!========================================================================================================
!THIS IS THE MAIN ROUTINE TO INTERPOLATE VELOCITIES ONTO THE LAGRANGIAN GRID

	subroutine boundary_interpolate(vel_horiz_grid,vel_vert_grid,loc_bnd,vel_bnd)

	!FIRST TWO INPUTS ARE ARRAYS WHERE THE GRID VELOCITIES WILL ARE STORED
	!NOTICE THAT THEY ARE DIFFERENT SIZES BECAUSE OF THE MAC GRID STRUCTURE
	REAL, DIMENSION(1:Nx,1:Ny) :: vel_horiz_grid, vel_vert_grid
	!LAST TWO INPUTS SHOULD BE ARRAYS DEFINING LOCATION OF LAGRANGIAN POINTS
	!AND VELOCITY ACTING AT THOSE LAGRANGIAN POINTS
	REAL, DIMENSION(Mb,2) :: loc_bnd, vel_bnd

	!VARIABLES FOR USE IN THE PROGRAM
	INTEGER :: I, left, bottom, innerI, innerK, xin, yin
	INTEGER, DIMENSION(1:4) :: xindeces, yindeces, modxin, modyin
	REAL, DIMENSION(1:4) :: coords, distance, xweights, yweights


	! begin by zeroing out the array we will spread the velocity into
	! POSSIBLY CHANGE THIS PIECE OF CODE, I'M NOT SURE
	vel_bnd = 0.0

	!loop through the length of the IB array
	do I = 1,Mb
	!FIRST WE WILL INTERPOLATE THE HORIZONTAL VELOCITIES FROM THE HORIZONTAL GRID
		!Determine the index corresponding to X-coordinate immediately below the IB point
		left = floor(loc_bnd(I,1)/hx + 1./2.)
		!Form indeces defining  points to spread to
		xindeces = (/left-1,left,left+1,left+2/)
		!Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modxin = modulo(xindeces-1,Nx)+1
		!Calculate the x-coordinate of those points
		coords = hx*(xindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - loc_bnd(I,1)

		!Calculate weights via the approximate delta function
		xweights = better_phi(distance/hx)/hx



		!Determine the index corresponding to Y-coordinate immediately below the IB point
		bottom = floor(loc_bnd(I,2)/hy + 1./2.)
		!Form indeces defining  points to spread to
		yindeces = (/bottom-1, bottom, bottom+1, bottom+2/)
		!Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modyin = modulo(yindeces-1,Ny)+1
		!Calculate the y-coordinate of those points
		coords = hy*(yindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - loc_bnd(I,2)

		!Calculate weights via the approximate delta function
		yweights = better_phi(distance/hy)/hy


		!Now we loop through the appropiate indeces and interpolate the weighted values
		do innerI = 1,4
			do innerK = 1,4
				!Add the velocity in the X-direction times the weights from X and Y
				!delta functions and scale by the Eulerian area difference
				vel_bnd(I,1) = vel_bnd(I,1) + vel_horiz_grid(modxin(innerI),modyin(innerK))*xweights(innerI)*yweights(innerK)*hx*hy
				vel_bnd(I,2) = vel_bnd(I,2) + vel_vert_grid(modxin(innerI),modyin(innerK))*xweights(innerI)*yweights(innerK)*hx*hy

			end do
		end do


!	!NOW WE WILL INTERPOLATE THE VERTICAL VELOCITIES FROM THE VERTICAL GRID
!		!Determine the index corresponding to Y-coordinate immediately below the IB point
!		bottom = floor(loc_bnd(I,2)/hy + 1./2.)
!		!Form indeces defining  points to spread to
!		yindeces = (/bottom-1,bottom,bottom+1,bottom+2/)
!		!Calculate the y-coordinate of those points
!		coords = hy*(yindeces - 1./2.)
!		!Calculate distance from these Eulerian points to IM point
!		distance = coords - loc_bnd(I,2)
!
!		!Calculate weights via the approximate delta function
!		yweights = better_phi(distance/hy)/hy
!
!
!
!		!Determine the index corresponding to X-coordinate immediately below the IB point
!		left = floor(loc_bnd(I,1)/hx + 1./2.)
!		!Form indeces defining  points to spread to
!		xindeces = (/left-1, left, left+1, left+2/)
!		!Calculate the X-coordinate of those points
!		coords = hx*(xindeces - 1./2.)
!		!Calculate distance from these Eulerian points to IM point
!		distance = coords - loc_bnd(I,1)
!
!		!Calculate weights via the approximate delta function
!		xweights = better_phi(distance/hx)/hx
!
!
!		!Now we loop through the appropiate indeces and interpolate the weighted values
!		do innerI = 1,4
!			do innerK = 1,4
!				!Determine the indeces in the Eulerian grid that we will spread to.
!				xin = mymod(xindeces(innerI),Nx)
!				yin = mymod(yindeces(innerK),Ny)
!				!Add the force in the Y-direction times the weights from X and Y delta functions
!				!and scale by the lagrangian arc length difference
!				vel_bnd(I,2) = vel_bnd(I,2) + vel_vert_grid(xin,yin)*xweights(innerI)*yweights(innerK)*hx*hy
!
!			end do
!		end do

	end do

	end subroutine boundary_interpolate


!========================================================================================================
!THIS IS THE MAIN ROUTINE TO UPDATE THE LAGRANGIAN POINTS BY PERFORMING TIME INTEGRATION

	subroutine boundary_time_integrate(loc_bnd_old,vel_bnd,loc_bnd_new)

	!FIRST INPUT IS M BY 2 ARRAY DEFINING THE 'OLD' BOUNDARY LOCATION
	!SECOND INPUT IS M BY 2 ARRAY DEFINING THE VELOCITY WITH WHICH BOUNDARY IS MOVING
	REAL, DIMENSION(1:Mb,1:2) :: loc_bnd_old, vel_bnd, loc_bnd_new
	!THIRD INPUT IS M BY 2 ARRAY IN WHICH THE 'NEW' BOUNDARY LOCATION WILL BE WRITTEN


	!We time step boundary using simple Foward Euler temporal discretization
	loc_bnd_new = loc_bnd_old + dt*vel_bnd

	end subroutine boundary_time_integrate


end module im_bnd_opers_mod
