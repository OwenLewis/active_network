! ============================================================================
! Name        : cortex_opers_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing operations that can be performed with the cortex (immersed network)
! ============================================================================

module cortex_opers_mod
USE global_param_mod
USE useful_functions_mod
implicit none
PRIVATE
PUBLIC :: network_grid_spread, network_interpolate, network_time_integrate
contains

!========================================================================================================
!THIS IS THE MAIN ROUTINE TO SPREAD FORCES ONTO THE EULERIAN GRID

	subroutine network_grid_spread(network,network_forces,force_horiz_grid,force_vert_grid)

	!FIRST TWO INPUTS SHOULD BE AN 'IM_NETWORK' TYPE DEFINING MY IMMERSED CORTEX
	!AND FORCES ACTING AT THE NODES OF SAID NETWORK
	type (im_network) :: network
	REAL, DIMENSION(:,:) :: network_forces
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



	!loop through the nodes of the im_network
	do I = 1,network%M_nodes
	!FIRST WE WILL SPREAD THE HORIZONTAL FORCES TO THE HORIZONTAL GRID
		!Determine the index corresponding to X-coordinate immediately below the CORTEX point
		left = floor(network%nodes(I)%location(1,1)/hx + 1./2)
		!Form indeces defining  points to spread to
		xindeces = (/left-1,left,left+1,left+2/)
		!Determine the indeces in the Eulerian grid that we will spread to.
		!Taking into account periodicity of the domain
        modxin = modulo(xindeces-1,Nx)+1
		!Calculate the x-coordinate of those points
		coords = hx*(xindeces - 1./2.)
		!Calculate distance from these Eulerian points to CORTEX point
		distance = coords - network%nodes(I)%location(1,1)
		!Calculate weights via the approximate delta function
		xweights = better_phi(distance/hx)/hx
!		print*, 'Sum of weights in x direction', sum(xweights)


		!Determine the index corresponding to Y-coordinate immediately below the CORTEX point
		bottom = floor(network%nodes(I)%location(1,2)/hy + 1./2.)
		!Form indeces defining  points to spread to
		yindeces = (/bottom-1, bottom, bottom+1, bottom+2/)
        !Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modyin = modulo(yindeces-1,Ny)+1
		!Calculate the y-coordinate of those points
		coords = hy*(yindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - network%nodes(I)%location(1,2)

		!Calculate weights via the approximate delta function
		yweights = better_phi(distance/hy)/hy
!		print*, 'Sum of weights in y direction', sum(yweights)


		!Now we loop through the appropiate indeces and spread the weighted values
		do innerI = 1,4
			do innerK = 1,4
				!Add the force in the X-direction times the weights from X and Y delta functions
				!and scale by the lagrangian area
				force_horiz_grid(modxin(innerI),modyin(innerK)) = &
                    &force_horiz_grid(modxin(innerI),modyin(innerK)) &
                    &+ network_forces(I,1)*xweights(innerI)*yweights(innerK)*network%nodes(I)%ref_area

                force_vert_grid(modxin(innerI),modyin(innerK)) = &
                    &force_vert_grid(modxin(innerI),modyin(innerK)) &
                    &+ network_forces(I,2)*xweights(innerI)*yweights(innerK)*network%nodes(I)%ref_area
			end do
		end do
!		print*, 'Total 2-D weighted force', foo
!		print*, 'Original force  (horiz) ', network_forces(I,1)

	!NOW WE WILL SPREAD THE VERTICAL FORCES TO THE VERTICAL GRID
		!Determine the index corresponding to Y-coordinate immediately below the CORTEX point
!		bottom = floor(network%nodes(I)%location(1,2)/hy + 1./2)
!		!Form indeces defining  points to spread to
!		yindeces = (/bottom-1,bottom,bottom+1,bottom+2/)
!		!Calculate the y-coordinate of those points
!		coords = hy*(yindeces - 1./2.)
!		!Calculate distance from these Eulerian points to CORTEX point
!		distance = coords - network%nodes(I)%location(1,2)
!
!		!Calculate weights via the approximate delta function
!		yweights = better_phi(distance/hy)/hy
!		print*, 'Sum of weights in y direction', sum(yweights)
!
!
!		!Determine the index corresponding to X-coordinate immediately below the CORTEX point
!		left = floor(network%nodes(I)%location(1,1)/hx + 1./2.)
!		!Form indeces defining  points to spread to
!		xindeces = (/left-1, left, left+1, left+2/)
!		!Calculate the X-coordinate of those points
!		coords = hx*(xindeces - 1./2.)
!		!Calculate distance from these Eulerian points to CORTEX point
!		distance = coords - network%nodes(I)%location(1,1)
!
!		!Calculate weights via the approximate delta function
!		xweights = better_phi(distance/hx)/hx
!		print*, 'Sum of weights in x direction', sum(xweights)
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
!				force_vert_grid(xin,yin) = force_vert_grid(xin,yin) &
!											&+ network_forces(I,2)*xweights(innerI)*yweights(innerK)*network%nodes(I)%ref_area
!			end do
!		end do
!		print*, 'Total 2-D weighted force', foo
!		print*, 'Original force  (vert)  ', network_forces(I,2)

	end do
	end subroutine network_grid_spread


!========================================================================================================
!THIS IS THE MAIN ROUTINE TO INTERPOLATE VELOCITIES ONTO THE IMMERSED CORTEX POINTS

	subroutine network_interpolate(vel_horiz_grid,vel_vert_grid,network,vel_cortex)

	!FIRST TWO INPUTS ARE ARRAYS WHERE THE GRID VELOCITIES WILL ARE STORED
	!NOTICE THAT THEY ARE DIFFERENT SIZES BECAUSE OF THE MAC GRID STRUCTURE
	REAL, DIMENSION(1:Nx,1:Ny) :: vel_horiz_grid, vel_vert_grid
	!THIRD INPUT SHOULD BE AN 'IM_NETWORK' STRUCTURE DEFINING THE CORTEX IMMERSED IN FLUID
	type(im_network) :: network
	REAL, DIMENSION(:,:) :: vel_cortex

	!VARIABLES FOR USE IN THE PROGRAM
	INTEGER :: I, left, bottom, innerI, innerK, xin, yin
	INTEGER, DIMENSION(1:4) :: xindeces, yindeces, modxin, modyin
	REAL, DIMENSION(1:4) :: coords, distance, xweights, yweights


	! begin by zeroing out the array we will spread the velocity into
	! POSSIBLY CHANGE THIS PIECE OF CODE, I'M NOT SURE
	vel_cortex = 0.0

	!loop through the length of the IB array
	do I = 1,network%M_nodes
	!FIRST WE WILL INTERPOLATE THE HORIZONTAL VELOCITIES FROM THE HORIZONTAL GRID
		!Determine the index corresponding to X-coordinate immediately below the IB point
		left = floor(network%nodes(I)%location(1,1)/hx + 1./2.)
		!Form indeces defining  points to spread to
		xindeces = (/left-1,left,left+1,left+2/)
		!Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modxin = modulo(xindeces-1,Nx)+1
		!Calculate the x-coordinate of those points
		coords = hx*(xindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - network%nodes(I)%location(1,1)

		!Calculate weights via the approximate delta function
		xweights = better_phi(distance/hx)/hx



		!Determine the index corresponding to Y-coordinate immediately below the IB point
		bottom = floor(network%nodes(I)%location(1,2)/hy + 1./2.)
		!Form indeces defining  points to spread to
		yindeces = (/bottom-1, bottom, bottom+1, bottom+2/)
		!Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modyin = modulo(yindeces-1,Ny)+1
		!Calculate the y-coordinate of those points
		coords = hy*(yindeces - 1./2.)
		!Calculate distance from these Eulerian points to IM point
		distance = coords - network%nodes(I)%location(1,2)

		!Calculate weights via the approximate delta function
		yweights = better_phi(distance/hy)/hy


		!Now we loop through the appropiate indeces and interpolate the weighted values
		do innerI = 1,4
			do innerK = 1,4
				!Add the velocity in the X-direction times the weights from X and Y
				!delta functions and scale by the Eulerian area difference
				vel_cortex(I,1) = vel_cortex(I,1) + &
				    &vel_horiz_grid(modxin(innerI),modyin(innerK))*xweights(innerI)*yweights(innerK)*hx*hy
                vel_cortex(I,2) = vel_cortex(I,2) + &
                    &vel_vert_grid(modxin(innerI),modyin(innerK))*xweights(innerI)*yweights(innerK)*hx*hy

			end do
		end do


!	!NOW WE WILL INTERPOLATE THE VERTICAL VELOCITIES FROM THE VERTICAL GRID
!		!Determine the index corresponding to Y-coordinate immediately below the IB point
!		bottom = floor(network%nodes(I)%location(1,2)/hy + 1./2.)
!		!Form indeces defining  points to spread to
!		yindeces = (/bottom-1,bottom,bottom+1,bottom+2/)
!		print*, 'Y indices are: ', yindeces
!		!Calculate the y-coordinate of those points
!		coords = hy*(yindeces - 1./2.)
!		print*, 'Which translates to y-coords:', coords
!		!Calculate distance from these Eulerian points to IM point
!		distance = coords - network%nodes(I)%location(1,2)
!		print*, 'That means distances to cortex point: ', distance
!
!		!Calculate weights via the approximate delta function
!		yweights = better_phi(distance/hy)/hy
!		print*, 'Weights from delta: ', yweights
!
!
!
!		!Determine the index corresponding to X-coordinate immediately below the IB point
!		left = floor(network%nodes(I)%location(1,1)/hx + 1./2.)
!		!Form indeces defining  points to spread to
!		xindeces = (/left-1, left, left+1, left+2/)
!		!Calculate the X-coordinate of those points
!		coords = hx*(xindeces - 1./2.)
!		!Calculate distance from these Eulerian points to IM point
!		distance = coords - network%nodes(I)%location(1,1)
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
!				vel_cortex(I,2) = vel_cortex(I,2) + vel_vert_grid(xin,yin)*xweights(innerI)*yweights(innerK)*hx*hy
!
!			end do
!		end do
!
	end do

	end subroutine network_interpolate

!========================================================================================================
!THIS IS THE MAIN ROUTINE TO TIME INTEGRATE THE POINTS OF THE IMMERSED NETWORK (CORTEX)
	subroutine network_time_integrate(network,velocity)
	!INPUTS ARE SIMPLY THE NETWORK BEING MOVED
	type (im_network) :: network
	!AS WELL AS THE VELOCITY FIELD THAT IT IS MOVING IN
	REAL, DIMENSION(:,:) :: velocity

	!IMPORTANT: velocity should be a 'M_nodes' x 2 array. OTHERWISE THIS SHIT WILL BE BROKEN!

	!This is a very simple forward euler discretization in time
	network%nodes(:)%location(1,1) = network%nodes(:)%location(1,1) + dt*velocity(:,1)
	network%nodes(:)%location(1,2) = network%nodes(:)%location(1,2) + dt*velocity(:,2)


	end subroutine network_time_integrate

!========================================================================================================


end module cortex_opers_mod
