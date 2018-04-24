! ============================================================================
! Name        : adhesion_opers_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing the internal functions of the cortex (immersed network)
! ============================================================================

module adhesion_opers_mod
USE global_param_mod
USE useful_functions_mod
implicit none
PRIVATE
PUBLIC :: adh_populate, adh_initialize, adh_time_integrate, adh_grid_spread
contains


!========================================================================================================

	!COMMENT THIS FUCKING SUBROUTINE YOU DUMBASS.
	subroutine adh_initialize(network,adhesions)
		!'network' is an im_network type that adhesions are connected to
		type (im_network) :: network
		!'adhesions' is an 'adh_complex' type that we are initializing
		type (adh_complex) :: adhesions
		!Just a dummy integer used for looping & such
		INTEGER :: I, K
		!Some dummy reals used for calculating the initialization adhesion strength
		real :: left, right

		!First, we determine how big the network (and associated adhesion complex) is

		K = network%M_nodes


		!Now allocate space within adhesion complex for all the necessary variables
		allocate(adhesions%location(1:K,1:2),adhesions%adh_dir(1:K,1:2),adhesions%adh_length(1:K),&
				&adhesions%adh_stiff(1:K),adhesions%adh_slip(1:K),adhesions%ref_location(1:K,1:2))

		!also allocate space for the indeces of connected points
		allocate(adhesions%network_buddy(1:K))
		!Adn set the total number of adhesions
		adhesions%M_adh = K


		!Set the parameters equal to their global inpus
		!And set the direction and length of every link to zero
		adhesions%adh_dir = 0
		adhesions%adh_length = 0
		adhesions%adh_stiff = adhesion_stiff
		adhesions%adh_slip = adhesion_slip

		left = minval(network%nodes(:)%location(1,1))
		right = maxval(network%nodes(:)%location(1,1))

!		! We loop through the number of network nodes (and thus adhesions)
		do I = 1,K
!			print*, 'We are dealing with node', I
!			print*, 'The network point is at      ', network%nodes(I)%location(1,:)
			!We will initially set all the adhesions to the same location as their nodes
			adhesions%location(I,:) = network%nodes(I)%location(1,:)
!			print*, 'And I placed the adhesion at ', adhesions%location(I,:)
			!And their 'buddy' has the same index that they do
			adhesions%network_buddy(I) = I


            !I'm not entirely sure why this code still exists. I really need to think about how my shit scales
            !I'm about 99% sure these lines are unnecessary. Forces are in units of force per unit
            !Lagrangian (network) area. When I spread to the substrate, I multiply by network area and
            !divide by lab (eulerian) area. Everything seems consistent & it means that my parameters
            !Dont need to be adjusted as I refine the grid. I'm pretty sure

!			adhesions%adh_slip(I) = adhesion_slip*((adhesions%location(I,1) - left)/(right - left) + 1)
!			adhesions%adh_slip(I) = adhesion_slip*network%nodes(I)%ref_area
!			adhesions%adh_stiff(I) = adhesion_stiff/network%nodes(I)%ref_area
		end do

		adhesions%ref_location = adhesions%location

	end subroutine adh_initialize

!========================================================================================================

	!COMMENT THIS ROUTINE YOU FUCKING DUMBASS
	subroutine adh_populate(network,adhesions)
		!'network' is an im_network type that adhesions are connected to
		type (im_network) :: network
		!'adhesions' is an 'adh_complex' type that we are initializing
		type (adh_complex) :: adhesions
		!Just a dummy integer used for looping & such
		INTEGER :: I, K
		!Vector used for calculations
		REAL, DIMENSION(1,1:2) :: ray
		REAL :: length


		!We will loop through all points in the adhesion complex

		do I = 1,adhesions%M_adh
			!First get the index that determines which network point this adhesion is connected to
			K = adhesions%network_buddy(I)

			!Here we just perform a check to make sure nothing went wront during the linking

			if (.false.) then

			end if

			!Calculate the ray from the cortex to the adhesion
			ray(1,:) = network%nodes(K)%location(1,:) - adhesions%location(I,:)
			!and its length
			length = sqrt(ray(1,1)**2 + ray(1,2)**2)

			!Normalize that ray
			ray(1,1) = ray(1,1)/length
			ray(1,2) = ray(1,2)/length

			!If that ray is zero (or close), i just created some NaN's
			!This is just a check statment to deal with them
			if (isNan(ray(1,1)).or.isNan(ray(1,2))) then
				ray = 0
			end if

			!Now write these quantities into the adhesion

			adhesions%adh_dir(I,:) = ray(1,:)
			adhesions%adh_length(I) = length


		end do

        !Now just calculate the force associated with links between network and adhesion complex
		call adh_force_calc(network,adhesions)

	end subroutine adh_populate


!========================================================================================================

	!YOU DIDN'T PROPERLY COMMENT ANY OF THESE ROUTINES YOU GODDAMN FOOL

	subroutine adh_force_calc(network,adhesions)
		!'network' is an im_network type that adhesions are connected to
		type (im_network) :: network
		!'adhesions' is an 'adh_complex' type that we are initializing
		type (adh_complex) :: adhesions
		!Just a dummy integer used for looping & such
		INTEGER :: I, K
		!Vector used for calculations
		REAL, DIMENSION(1,1:2) :: force


		!First, zero out all the adhesion forces, just in case
		network%nodes(:)%adh_force(1,1) = 0.
		network%nodes(:)%adh_force(1,2) = 0.

		!And lets just zero out that strain energy too
		adhesions%strain_energy = 0.

		!Now loop through each of the adhesions
		do I = 1,adhesions%M_adh

			K = adhesions%network_buddy(I)

			!Now calculate the force in the link
			force(1,1) = adhesions%adh_stiff(I)*adhesions%adh_length(I)*adhesions%adh_dir(I,1)
			force(1,2) = adhesions%adh_stiff(I)*adhesions%adh_length(I)*adhesions%adh_dir(I,2)

            !And while we're at it, lets update the strain energy
			adhesions%strain_energy = adhesions%strain_energy + &
                                               & network%nodes(K)%ref_area*adhesions%adh_stiff(I)*(adhesions%adh_length(I)**2)/2


			!NOTICE THAT HERE I DO NOT DIVIDE THE FORCE BY REFERENCE AREA. THIS HAS IMPLICATIONS FOR
			!THE UNITS OF ADH_STIFF, BUT MORE IMPORTANTLY IT HAS IMPLICATIONS FOR HOW I CALCULATE U_SUB
			!I SHOULD REVISIT THIS IN THE FUTURE AFTER MY CONFERENCE.
			network%nodes(K)%adh_force = -force!/network%nodes(K)%ref_area
		end do
	end subroutine adh_force_calc


!========================================================================================================

	subroutine adh_time_integrate(adhesions,velocity)
		!INPUTS ARE SIMPLY THE NETWORK BEING MOVED
	type (adh_complex) :: adhesions
	!AS WELL AS THE VELOCITY FIELD THAT IT IS MOVING IN
	REAL, DIMENSION(:,:) :: velocity

	!IMPORTANT: velocity should be a 'M_adh' x 2 array. OTHERWISE THIS SHIT WILL BE BROKEN!

	!This is a very simple forward euler discretization in time
	adhesions%location(:,1) = adhesions%location(:,1) + dt*velocity(:,1)
	adhesions%location(:,2) = adhesions%location(:,2) + dt*velocity(:,2)

	end subroutine adh_time_integrate

!========================================================================================================

!========================================================================================================
!THIS IS THE MAIN ROUTINE TO SPREAD FORCES ONTO THE EULERIAN GRID

    subroutine adh_grid_spread(adhesions,network,adh_forces,force_horiz_subs,force_vert_subs)

    !FIRST INPUT SHOULD BE AN 'ADH_COMPLEX' DEFINING WHERE THE ADHESION FORCES TO BE SPREAD ARE LOCATED
    type (adh_complex) :: adhesions
    !NEXT TWO INPUTS SHOULD BE AN 'IM_NETWORK' TYPE DEFINING MY IMMERSED CORTEX
    !AND FORCES ACTING AT THE NODES OF PREVIOUSLY MENTIONED 'ADH_COMPLEX'
    type (im_network) :: network
    REAL, DIMENSION(:,:) :: adh_forces
    !LAST TWO INPUTS ARE ARRAYS WHERE THE SPREAD FORCES WILL BE WRITTEN
    !NOTICE THAT THEY ARE DIFFERENT SIZES BECAUSE OF THE MAC GRID STRUCTURE
    REAL, DIMENSION(1:Nx,1:Ny) :: force_horiz_subs
    REAL, DIMENSION(1:Nx,1:Ny) :: force_vert_subs
    !VARIABLES FOR USE IN THE PROGRAM
    INTEGER :: I, left, bottom, innerI, innerK
    INTEGER, DIMENSION(1:4) :: xindeces, yindeces, modxin, modyin
    REAL, DIMENSION(1:4) :: coords, distance, xweights, yweights

    ! begin by zeroing out the array we will spread the force into
    ! POSSIBLY CHANGE THIS PIECE OF CODE, I'M NOT SURE
    force_horiz_subs = 0.0
    force_vert_subs = 0.0



    !loop through the nodes of the adh_complex
    do I = 1,adhesions%M_adh
    !FIRST WE WILL SPREAD THE HORIZONTAL FORCES TO THE HORIZONTAL GRID
        !Determine the index corresponding to X-coordinate immediately below the CORTEX point
        left = floor(adhesions%location(I,1)/hx + 1./2.)
        !Form indeces defining  points to spread to
        xindeces = (/left-1,left,left+1,left+2/)
        !Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modxin = modulo(xindeces-1,Nx)+1
        !Calculate the x-coordinate of those points
        coords = hx*(xindeces - 1./2.)
        !Calculate distance from these Eulerian points to CORTEX point
        distance = coords - adhesions%location(I,1)
        !Calculate weights via the approximate delta function
        xweights = better_phi(distance/hx)/hx
!       print*, 'Sum of weights in x direction', sum(xweights)


        !Determine the index corresponding to Y-coordinate immediately below the CORTEX point
        bottom = floor(adhesions%location(I,2)/hy + 1./2.)
        !Form indeces defining  points to spread to
        yindeces = (/bottom-1, bottom, bottom+1, bottom+2/)
        !Determine the indeces in the Eulerian grid that we will spread to.
        !Taking into account periodicity of the domain
        modyin = modulo(yindeces-1,Ny)+1
        !Calculate the y-coordinate of those points
        coords = hy*(yindeces - 1./2.)
        !Calculate distance from these Eulerian points to IM point
        distance = coords - adhesions%location(I,2)

        !Calculate weights via the approximate delta function
        yweights = better_phi(distance/hy)/hy
!       print*, 'Sum of weights in y direction', sum(yweights)


        !Now we loop through the appropiate indeces and spread the weighted values
        do innerI = 1,4
            do innerK = 1,4
                !Add the force in the X-direction times the weights from X and Y delta functions
                !and scale by the lagrangian area
                force_horiz_subs(modxin(innerI),modyin(innerK)) = force_horiz_subs(modxin(innerI),modyin(innerK)) &
                    &+ adh_forces(I,1)*xweights(innerI)*&
                    &yweights(innerK)*network%nodes(adhesions%network_buddy(I))%ref_area

                force_vert_subs(modxin(innerI),modyin(innerK)) = force_vert_subs(modxin(innerI),modyin(innerK)) &
                    &+ adh_forces(I,2)*xweights(innerI)*&
                    &yweights(innerK)*network%nodes(adhesions%network_buddy(I))%ref_area
            end do
        end do!End the loop over the local patch defined by the adhesive point in question


    end do!End the loop through all adhesive points
    end subroutine adh_grid_spread


end module adhesion_opers_mod
