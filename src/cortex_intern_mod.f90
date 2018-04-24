! ============================================================================
! Name        : cortex_intern_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing the internal functions of the cortex (immersed network)
! ============================================================================

module cortex_intern_mod
USE global_param_mod
implicit none
REAL :: strain_energy
PRIVATE
PUBLIC :: network_populate, network_initialize, linkup, link_force_calc, vec_network_integrate, scalar_network_integrate
contains


!========================================================================================================

	!We will now update all of the various quantities associated with a network
	!Assuming that the node points have been updated through a fluid step or some other means
	subroutine network_initialize(network)
		!'network' is an im_network type that we will be updating
		type (im_network) :: network
		!Just a dummy integer used for looping & such
		INTEGER :: I, point1, point2, point3, K
		REAL :: area, length, avg_area = 0, left, right
		REAL, DIMENSION(1:1,2) :: centroid, center, ray

        !These 2 lines of code are used when setting the spatially dependent spring constants
        !of the various links. To express everything in body coordinates, we have to calculate
        !the right- and left-most points in the network
        right = maxval(network%nodes(:)%location(1,1))
        left = minval(network%nodes(:)%location(1,1))

		!We will loop through all the faces, updating as we go

		do I = 1,network%M_faces
			!Grab the integers that define the nodes of this face
			point1 = network%faces(I)%my_nodes(1)
			point2 = network%faces(I)%my_nodes(2)
			point3 = network%faces(I)%my_nodes(3)

			!Calculate the area of this face by taking a cross product of the coordinates
			!of the bounding nodes
			area = abs((network%nodes(point1)%location(1,1)*network%nodes(point2)%location(1,2)&
				&-network%nodes(point2)%location(1,1)*network%nodes(point1)%location(1,2)&
				&+network%nodes(point2)%location(1,1)*network%nodes(point3)%location(1,2)&
				&-network%nodes(point3)%location(1,1)*network%nodes(point2)%location(1,2)&
				&+network%nodes(point3)%location(1,1)*network%nodes(point1)%location(1,2)&
				&-network%nodes(point1)%location(1,1)*network%nodes(point3)%location(1,2))/2.)

			!Now calculate the centroid (its just the average of the points)
			centroid = (network%nodes(point1)%location + network%nodes(point2)%location&
						&+ network%nodes(point3)%location)/3

			!Now write these two quantities into this face
			network%faces%ref_area = area
			network%faces(I)%ref_center = centroid
			network%nodes(point1)%ref_area = network%nodes(point1)%ref_area + area/3
			network%nodes(point2)%ref_area = network%nodes(point2)%ref_area + area/3
			network%nodes(point3)%ref_area = network%nodes(point3)%ref_area + area/3

		end do

		!Now we will loop through all of the links in the network
		do I = 1,network%M_links
			!Grab the integers that point to the nodes at either end of the link
			point1 = network%links(I)%my_nodes(1)
			point2 = network%links(I)%my_nodes(2)

			!Calculate the area associated with this link
			avg_area = (network%nodes(point1)%ref_area + network%nodes(point2)%ref_area)*0.5

			!Calculate the length of this link
			length = sqrt((network%nodes(point1)%location(1,1) - network%nodes(point2)%location(1,1))**2&
					&+(network%nodes(point1)%location(1,2) - network%nodes(point2)%location(1,2))**2)

			!Now calculate the center of this link (ITS EASY)
			center = (network%nodes(point1)%location(:,:)+network%nodes(point2)%location(:,:))/2

			!Now write these quantities into this link
			network%links(I)%ref_length = length
			network%links(I)%ref_center = center

			!We also set the LOCAL spring constant in order to preserve scaling
			!We have quite a few choices here
            select case(homogeneous)
                case(1)
                    !This is the standard (SPATIALLY UNIFORM) method where spring constants are set
                    !to reproduce the desired lame constant
                    network%links(I)%stiff = 8*lame*avg_area/(3*length)

    		    case(0)
                    !This is the modified (SPATIALLY DEPENDENT) method where spring constants are set like a step function
                    !to reproduce the desired lame constant
                    network%links(I)%stiff = merge(1.,head_multiple,center(1,1).lt.0.8)*8*lame*avg_area/(3*length)

                case(2)
			        !THIS IS A MODIFIED BLOCK OF CODE WHERE WE ARE SETTING THE SPRING CONSTANTS TO BE SPATIALLY
			        !DEPENDENT YOU KNOW, BECAUSE SHIT IS AWESOME LIKE THAT
                    network%links(I)%stiff = network%links(I)%stiff*&
                        &(1. - (1.-head_multiple)*(network%links(I)%ref_center(1,1)-left)/(right-left))
            end select

			!Finally, we attempt to identify if we are at the boundary of the network
			!K identifies if this link has '0' as one of its faces ('0' is the exterior face)
			K = count(mask=network%links(I)%my_faces.eq.0)
			!If this is true
			if (K.gt.0) then
				!We identify this link and both its nodes as boundary elements
				network%links(I)%boundary = .true.
				network%nodes(point1)%boundary = .true.
				network%nodes(point2)%boundary = .true.
			end if

		end do

		!Here i'm going to initialize some node based quantities
                do I = 1,network%M_nodes
                   !Here I am scaling the spring constant with which a boundary node is connected to the membrane
                   !Based on how many connections it has to the interior network. This was done to prevent boundary kinking as the 
                   !Grid is refined.
                   network%nodes(I)%connect_k = connect_stiff*network%nodes(I)%link_count

                   !Here I am initializing the drag parameter and 'center' value at immersed network nodes.
                   network%nodes(I)%node_drag = network_slip
                   center = network%nodes(I)%location

                   !here are some optional statements to build the FLOW CHANNEL (2 different geometries) into the cell.
                   !if ((9*(center(1,1)-0.5)**2 + 144*(center(1,2)-0.5)**2).lt.1) network%nodes(I)%node_drag = network_slip/2.0
                   !if (((center(1,1).gt.0.25).and.(center(1,2).lt.9./16.).and.(center(1,2).gt.7./16.)).or.(center(1,1).gt.0.75)) network%nodes(I)%node_drag = network_slip/10.
                end do


	end subroutine network_initialize

!========================================================================================================

	!We will now update all of the various quantities associated with a network
	!Assuming that the node points have been updated through a fluid step or some other means
	subroutine network_populate(network)
		!'network' is an im_network type that we will be updating
		type (im_network) :: network
		!Just a dummy integer used for looping & such
		INTEGER :: I, point1, point2, point3
		REAL :: area, length
		REAL, DIMENSION(1:1,2) :: centroid, center, ray

		!We will loop through all the faces, updating as we go

		!First, zero out all the nodes' areas. These will be updated in the loop
		network%nodes(:)%area = 0

		!Also, zero out the strain energy
		network%strain_energy = 0.
		
		do I = 1,network%M_faces
			!Grab the integers that define the nodes of this face
			point1 = network%faces(I)%my_nodes(1)
			point2 = network%faces(I)%my_nodes(2)
			point3 = network%faces(I)%my_nodes(3)

			!Calculate the area of this face by taking a cross product of the coordinates
			!of the bounding nodes
			area = abs((network%nodes(point1)%location(1,1)*network%nodes(point2)%location(1,2)&
				&-network%nodes(point2)%location(1,1)*network%nodes(point1)%location(1,2)&
				&+network%nodes(point2)%location(1,1)*network%nodes(point3)%location(1,2)&
				&-network%nodes(point3)%location(1,1)*network%nodes(point2)%location(1,2)&
				&+network%nodes(point3)%location(1,1)*network%nodes(point1)%location(1,2)&
				&-network%nodes(point1)%location(1,1)*network%nodes(point3)%location(1,2))/2.)

			!Now calculate the centroid (its just the average of the points)
			centroid = (network%nodes(point1)%location + network%nodes(point2)%location&
						&+ network%nodes(point3)%location)/3

			!Now write these two quantities into this face
			network%faces(I)%area = area
			network%faces(I)%face_center = centroid

			!Also distribute this area to the bounding nodes
			network%nodes(point1)%area = network%nodes(point1)%area + area/3
			network%nodes(point2)%area = network%nodes(point2)%area + area/3
			network%nodes(point3)%area = network%nodes(point3)%area + area/3
		end do


		!Now we will loop through all of the links in the network
		do I = 1,network%M_links
			!Grab the integers that point to the nodes at either end of the link
			point1 = network%links(I)%my_nodes(1)
			point2 = network%links(I)%my_nodes(2)
			!Calculate the vector that defines this link
			ray = network%nodes(point2)%location - network%nodes(point1)%location
			!Calculate the length of this link
			length = sqrt((network%nodes(point1)%location(1,1) - network%nodes(point2)%location(1,1))**2&
					&+(network%nodes(point1)%location(1,2) - network%nodes(point2)%location(1,2))**2)


			!Now calculate the center of this link (ITS EASY)
			center = (network%nodes(point1)%location(:,:)+network%nodes(point2)%location(:,:))/2

			!Now write these quantities into this link
			network%links(I)%rayup = ray
			network%links(I)%raydown = -ray
			network%links(I)%link_length = length
			network%links(I)%link_center = center
		end do

		call elastic_force_calc(network)

	end subroutine network_populate


!========================================================================================================

	!We will now update all of the various quantities associated with a network
	!Assuming that the node points have been updated through a fluid step or some other means
	subroutine elastic_force_calc(network)
		!'network' is an im_network type that we will be updating
		type (im_network) :: network


		!DUMMY VARIABLES USED FOR CALCULATING FORCE IN EACH LINK
		REAL :: strain, tens
		REAL, DIMENSION(1,1:2) :: direction, scaled_force
		!These are dummy variables for looping and such
		INTEGER :: I, point1, point2

		!First we will make sure that forces are zerod out. KEEP THIS IN MIND WHEN WRITING CALING ROUTINE
		network%nodes(:)%elastic_force(1,1) = 0
		network%nodes(:)%elastic_force(1,2) = 0
!		print*, 'forces has size', size(forces,1), size(forces,2)

		!Now loop through each of the links
		do I = 1,network%M_links
			!Grab the indices of the two points this link connects
			point1 = network%links(I)%my_nodes(1)
			point2 = network%links(I)%my_nodes(2)
!			print*, 'This link connects nodes', point1, 'and', point2

			!Now calculate the strain in the link
			strain = network%links(I)%link_length - network%links(I)%ref_length
!			print*, 'It has this much dimensional strain', strain
			!Normalize by reference length to make sure 'strain' is actually a strain
			strain = strain/network%links(I)%ref_length
!			print*, 'It has this much non-dimensional strain', strain


			!Calculate the strain energy in this link and add it to the total for the network
			network%strain_energy = network%strain_energy + (strain**2)*network%links(I)%stiff*network%links(I)%ref_length/2

			!'tens' is just the 'spring constant' (scaling taken care of perviously)
			!times the strain which we just calculate
			tens = network%links(I)%stiff*strain
			!tens = network%links(I)%stiff*(strain + network%links(I)%active_tens)

!			print*, 'Which means this much tension', tens
			!Now we need to put the properly scaled forces in the correct places
				!Lets deal with node number 1 first
			!First, we make a UNIT vector defining the direction of force
			direction = network%links(I)%rayup
!			print*, 'vector connecting the two nodes:', direction
			direction = direction/network%links(I)%link_length
!			print*, 'Now its normalized:', direction
!			print*, 'Has size', size(direction,1), size(direction,2)
			!Now we multiply this unit vector by the 'tens' and scale it by the area of node 1
			scaled_force = (tens+network%links(I)%active_tens)*direction/network%nodes(point1)%ref_area
!			scaled_force(1,2) = direction(1,2)*network%links(I)%active_tens/network%nodes(point1)%ref_area

			!Now add this force (per unit area) to the entry corresponding to node 1
			network%nodes(point1)%elastic_force = network%nodes(point1)%elastic_force + scaled_force

			!Now we deal with node number 2

			!Simply reverse the direction
			direction = -direction

			!Scale it by 'tens' and the area of NODE 2 THIS TIME
			scaled_force = (tens+network%links(I)%active_tens)*direction/network%nodes(point2)%ref_area
!			scaled_force(1,2) = direction(1,2)*network%links(I)%active_tens/network%nodes(point2)%ref_area

			!Now add this force (per unit area) to the entry corresponding to node 2
			network%nodes(point2)%elastic_force = network%nodes(point2)%elastic_force + scaled_force

		end do
	end subroutine elastic_force_calc


!========================================================================================================

	subroutine linkup(network,boundary)
		!subroutine takes as input a network and a boundary
		!THIS IS IMPORTANT: THESE SHOULD BE IN A CONFIGURATION WHERE IT IS OBVIOUS WHO GETS LINKED TO WHOM
		!THAT IS TO SAY, POINTS TO BE LINKED SHOULD HAVE THE SAME (OR VERY CLOSE) LOCATION
		type(im_network) :: network
		type(im_bnd) :: boundary

		!Arrays to construct distance from cortex point to immersed boundary point
		REAL, DIMENSION(1:Mb,1:2) :: difference
		REAL, DIMENSION(1:Mb) :: distance
		!This is a real for 'tolerance'. i haven't quite figured out how to use this just yet.
		REAL :: tol = 10**(-4)

		!An integer for looping and an integer to store the index of linking that I find
		INTEGER :: I
		INTEGER, DIMENSION(1) :: A

		!We are going to loop through all of the nodes of the immersed network
		do I = 1,network%M_nodes
			!We only wish to 'linkup' the boundary nodes
			if (network%nodes(I)%boundary) then
!				print*, 'Node number ', I, 'is a boundary node'
				!We construct an array of vectors representing the ray from the
				!cortex node to the immersed boundary points
				difference(:,1) = boundary%im_bnd_loc(:,1) - network%nodes(I)%location(1,1)
				difference(:,2) = boundary%im_bnd_loc(:,2) - network%nodes(I)%location(1,2)

				!Turn those vectors into distances
				distance = sqrt(difference(:,1)**2 + difference(:,2)**2)
!				print*, 'Distances', distance

				!Now we pick out the index of the immersed boundary point which is closes to the
				!cortex point under examination (i.e. the "I'th" cortex point

				A = minloc(distance)
!				print*, 'Minloc returns int of rank', kind(minloc(distance))!,mask=distance.lt.tol)
!				print*, 'This node is connected to boundary point', A
!				print*, 'They are ', minval(distance), 'apart'

				!Now we set this cortex point's "buddy" to be A and mark it as linked
				network%nodes(I)%connect_buddy = A(1)
				network%nodes(I)%linked = .true.
				!We also set the A'th boundary points buddy to be this cortex point
				boundary%connect_buddy(A(1)) = I
			end if
		end do

	end subroutine linkup

	!========================================================================================================

	subroutine link_force_calc(network,boundary)
		!subroutine takes as input a network and a boundary
		!THIS IS IMPORTANT: THESE SHOULD BE IN A CONFIGURATION WHERE IT IS OBVIOUS WHO GETS LINKED TO WHOM
		!THAT IS TO SAY, POINTS TO BE LINKED SHOULD HAVE THE SAME (OR VERY CLOSE) LOCATION
		type(im_network) :: network
		type(im_bnd) :: boundary

		!Arrays to construct distance from cortex point to immersed boundary point
		REAL, DIMENSION(1,1:2) :: difference, force, direction
		REAL :: distance


		!An integer for looping and an integer to store the index of linking
		INTEGER :: I, K

		network%nodes(:)%connect_force(1,1) = 0
		network%nodes(:)%connect_force(1,2) = 0
		boundary%connect_force = 0
		!We are going to loop through all of the nodes of the immersed network
		do I = 1,network%M_nodes
			!We only wish to 'linkup' the boundary nodes
			if (network%nodes(I)%linked) then
!				print*, 'Node number ', I, 'is a boundary node'
				!We construct an array of vectors representing the ray from the
				!cortex node to the immersed boundary points
				K = network%nodes(I)%connect_buddy
				difference(1,1) = boundary%im_bnd_loc(K,1) - network%nodes(I)%location(1,1)
				difference(1,2) = boundary%im_bnd_loc(K,2) - network%nodes(I)%location(1,2)

				!Turn those vectors into distances
				distance = sqrt(difference(1,1)**2 + difference(1,2)**2)

				direction = difference/distance
				if (distance.eq.0) direction = 0

				force = network%nodes(I)%connect_k*distance*direction

				network%nodes(I)%connect_force = force*ref_bnd%ds0(K)

				boundary%connect_force(K,:) = -force(1,:)*network%nodes(I)%ref_area!*hx*hy

			end if
		end do

	end subroutine link_force_calc
!========================================================================================================
!=============================================================================================
!THIS SUBROUTINE INTEGRATES A VECTOR QUANTITY ON THE IMMERSED NETWORK. THE NETWORK
!MUST BE SPECIFIED SO THAT IT KNOWS WHAT AREA TO INTEGRATE WITH RESPECT TO
subroutine vec_network_integrate(out,in,network)
!'NETWORK' is a of type 'im_network' and supplies the area needed to integrate over
type(im_network) :: network
!THE 'IN'-PUT IS A VECTOR QUANTITY DEFINED ON SOME IMMERSED BOUNDARY
REAL, DIMENSION(:,:) :: in
!THE 'OUT'-PUT IS A LENGTH 1 VECTOR QUANTITY THAT WILL BE CALCULATED
REAL, DIMENSION(1,1:2) :: out

!we simply multiply the input by the network's reference area & sum them
out(1,1) = sum(in(:,1)*network%nodes(:)%ref_area)
out(1,2) = sum(in(:,2)*network%nodes(:)%ref_area)
!AND WE'RE DONE

end subroutine vec_network_integrate


!=============================================================================================
!THIS SUBROUTINE INTEGRATES A SCALAR QUANTITY ON THE IMMERSED NETWORK. THE NETWORK
!MUST BE SPECIFIED SO THAT IT KNOWS WHAT AREA TO INTEGRATE WITH RESPECT TO
subroutine scalar_network_integrate(out,in,network)
!'NETWORK' is a of type 'im_network' and supplies the area needed to integrate over
type(im_network) :: network
!THE 'IN'-PUT IS A SCALAR QUANTITY DEFINED ON SOME IMMERSED BOUNDARY
REAL, DIMENSION(:) :: in
!THE 'OUT'-PUT IS A LENGTH 1 SCALAR QUANTITY THAT WILL BE CALCULATED
REAL :: out

!we simply multiply the input by the network's reference area & sum them
out = sum(in(:)*network%nodes(:)%ref_area)
!AND WE'RE DONE

end subroutine scalar_network_integrate


end module cortex_intern_mod
