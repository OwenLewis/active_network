! ============================================================================
! Name        : io_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing the input/output interface of the cortex (immersed network)
! ============================================================================

module io_mod
USE global_param_mod
implicit none
PRIVATE

PUBLIC :: network_read_in, network_write_out, bnd_read_in, bnd_write_out, fluid_write_out, adh_write_out,&
		  &substrate_write_out, lag_force_write_out


contains
	!THESE ARE THE ROUTINES TO WRITE OUT FLUID VARIABLES (2 VELOCITIES AND PRESSURE) AND
	!IMMERSED BOUNDARY VARIABLES (CURRENTLY JUST THE POSITION OF THE IMMERSED BOUNDARY, LOCATION
	!OF THE TETHER POINTS, AND THE 2 FORCES ASSOCIATED WITH THAT CONFIGURATION).


	!========================================================================================================
    !This subroutine is the one that writes out fluid variables that are passed to it.
    subroutine substrate_write_out(force_horiz,force_vert,plotcount)
        !This integer defines what plot we are on
        INTEGER :: plotcount
        !Forces are defined on a collocated eulerian grid
        REAL, DIMENSION(1:Nx,1:Ny) :: force_horiz,force_vert


        INTEGER :: A, B
        CHARACTER(len = 30) :: filename
            write(filename,'(2a,i4.4)') runname,'.substrate.', plotcount
            print*, filename
            !THIS IS THE CODE THAT WRITES TO A BINARY FILE. THIS ONLY WORKS ON AN
!           !IFORT COMPILE.
!           open(12,file=filename,form='binary')
!           do A=1,Ny
!               do B =1,Nx
!                   write(12) force_horiz(B,A), force_vert(B,A)
!               end do
!           end do
!           close(12)

            !THIS IS THE CODE THAT WRITES TO A BINARY FILE. THIS IS A SHITTY
            !WORKAROUND FOR GFORTRAN COMPILED EXECUTABLES
            open(12,file=filename,form='formatted')
            do A=1,Nx
                do B=1,Ny
                    write(12, '(E25.12,E25.12)'), &
                    & force_horiz(A,B),force_vert(A,B)
                end do
            end do
            close(12)

    end subroutine substrate_write_out

!========================================================================================================
	!This subroutine is the one that writes out fluid variables that are passed to it.
	subroutine fluid_write_out(vel_horiz,vel_vert,pressure,plotcount)
		!This integer defines what plot we are on
		INTEGER :: plotcount
		!All quantities are defined on collocated eulerian grid
		REAL, DIMENSION(1:Nx,1:Ny) :: vel_horiz, vel_vert, pressure

		INTEGER :: A, B
		CHARACTER(len = 30) :: filename
		    write(filename,'(2a,i4.4)') runname,'.fluid.', plotcount
			print*, filename
		    !THIS IS THE CODE THAT WRITES TO A BINARY FILE. THIS ONLY WORKS ON AN
!		    !IFORT COMPILE.
!  			open(20,file=filename,form='binary')
!  			do A=1,Ny
!  				do B =1,Nx
!        			write(20) vel_horiz(B,A), vel_vert(B,A), pressure(B,A)
!  				end do
!  			end do
!  			close(20)

  			!THIS IS THE CODE THAT WRITES TO A BINARY FILE. THIS IS A SHITTY
  			!WORKAROUND FOR GFORTRAN COMPILED EXECUTABLES
			open(20,file=filename,form='formatted')
  			do A=1,Nx
  				do B=1,Ny
	        		write(20, '(E25.12,E25.12,E25.12)'), &
    	    		& vel_horiz(A,B),vel_vert(A,B),pressure(A,B)
       			end do
       		end do
  			close(20)

	end subroutine fluid_write_out


!========================================================================================================
	!THIS SUBROUTINE IS SIMPLY PASSED THE INT 'PLOTCOUNT'
	!IT MUST HAVE ACCESS TO THE STRING 'RUNNAME' IN ORDER TO CONSTRUCT A FILE
	!TO WRITE OUT THE CURRENT IMMERSED BOUNDARY CONFIGURATION
	subroutine bnd_write_out(boundary,plotcount)
		!This integer that is passed in defines which output were on
		INTEGER :: plotcount
		type(im_bnd) :: boundary

		INTEGER :: A
		CHARACTER(len = 30) :: filename

		write(filename,'(2a,i4.4)') runname,'.bnd.', plotcount
		print*, filename

		!THIS IS THE VERSION WHERE WE WRITE OUT AS TEXT
  			open(30,file=filename,form='formatted')
  			do A=1,Mb
        		write(30, '(E25.12,E25.12)')&!,E25.12,E25.12,E25.12,E25.12,E25.12,E25.12,I4)'), &
        		& boundary%im_bnd_loc(A,1),boundary%im_bnd_loc(A,2)!,&
!        		& boundary%tether_pt_loc(A,1),boundary%tether_pt_loc(A,2),&
!        		& boundary%stretch_force(A,1),boundary%stretch_force(A,2),&
!        		& boundary%tether_force(A,1),boundary%tether_force(A,2),&
!        		& boundary%flag(A)
  			end do

  			close(30)

		!THIS IS THE VERSION WHERE WE WRITE OUT AS BINARY
!		open(30,file=filename,form='binary')
!		do A=1,M
!			write(30) boundary%im_bnd_loc(A,1), boundary%im_bnd_loc(A,2),&
!					& boundary%tether_pt_loc(A,1), boundary%tether_pt_loc(A,2),&
!					& boundary%stretch_force(A,1), boundary%stretch_force(A,2),&
!					& boundary%tether_force(A,1), boundary%tether_force(A,2), &
!                    & boundary%stretch_kp(A), boundary%tether_k(A), &
!					& real(boundary%flag(A))
!		end do
!		close(30)
	end subroutine bnd_write_out


!========================================================================================================

!This is the routine to write an immersed network type out to a file for visualization in MATLAB
	subroutine network_write_out(network,plotcount)
		!The inputs are an 'im_network' type called network (i'm clever, I know)
		type(im_network) :: network
		INTEGER :: plotcount
		!and a string that defines the file to be written out to
		CHARACTER(len = 30) :: filename

		!and here is just a nice little integer for looping purposes
		INTEGER :: I


		write(filename,'(2a,i4.4)') runname,'.network.', plotcount
		print*, filename

		open(80,file=filename,form='formatted')

	  	write(80, '(i8)'), network%M_nodes
	  	write(80, '(i8)'), network%M_faces

	  	!First we loop through all the nodes of the network, writing their coordinates
	  	do I = 1,network%M_nodes
			write(80,'(E25.12,E25.12)'), network%nodes(I)%location(1,1), network%nodes(I)%location(1,2)
	  	end do

	  	!Now we loop through all the faces of the network, writing their vertices
	  	do I = 1,network%M_faces
			write(80,'(i5,i5,i5)'), network%faces(I)%my_nodes(1), network%faces(I)%my_nodes(2), network%faces(I)%my_nodes(3)
	  	end do
	  	!close that file and be done with it!
		close(80)

        !This is an adendum to write out the elastic strain energy associated with the network
		write(filename,'(2a)') runname,'.networkenergy'
		open(01,file=filename,form='formatted',position='append')
		write(01,'(E25.12)'),network%strain_energy

	end subroutine network_write_out


!========================================================================================================

!This is the routine to write forces associated with a network type out to a file for visualization in MATLAB
    subroutine lag_force_write_out(network,plotcount)
        !The inputs are an 'im_network' type called network (i'm clever, I know)
        type(im_network) :: network
        INTEGER :: plotcount
        !and a string that defines the file to be written out to
        CHARACTER(len = 30) :: filename

        !and here is just a nice little integer for looping purposes
        INTEGER :: I


        write(filename,'(2a,i4.4)') runname,'.lag_force.', plotcount
        print*, filename

        open(40,file=filename,form='formatted')

        write(40, '(i8)'), network%M_nodes

        !First we loop through all the nodes of the network, writing the elastic forces, followed by connection force
        !followed by adhesive forces
        do I = 1,network%M_nodes
            write(40,'(E25.12,E25.12,E25.12,E25.12,E25.12,E25.12)'), network%nodes(I)%elastic_force(1,1), &
            &network%nodes(I)%elastic_force(1,2), network%nodes(I)%adh_force(1,1), network%nodes(I)%adh_force(1,2),&
            &network%nodes(I)%connect_force(1,1), network%nodes(I)%connect_force(1,2)
        end do

        !close that file and be done with it!
        close(40)

    end subroutine lag_force_write_out


!========================================================================================================

!This is the routine to write an adhesion complex type out to a file for visualization in MATLAB
	subroutine adh_write_out(adhesions,plotcount)
		!The inputs are an 'adh_complex' type called adhesions (i'm clever, I know)
		type(adh_complex) :: adhesions
		INTEGER :: plotcount
		!and a string that defines the file to be written out to
		CHARACTER(len = 30) :: filename

		!and here is just a nice little integer for looping purposes
		INTEGER :: I

		!construct the name of teh file to be written to
		write(filename,'(2a,i4.4)') runname,'.adhesion.', plotcount
		print*, filename
		!Open that file
		open(90,file=filename,form='formatted')

		!first, write out the number of adhesions in the complex
	  	write(90, '(i8)'), adhesions%M_adh

	  	!First we loop through all the nodes of the network, writing their coordinates and the network point they link to
	  	do I = 1,adhesions%M_adh
			write(90,'(E25.12,E25.12,i8)'), adhesions%location(I,1), adhesions%location(I,2), adhesions%network_buddy(I)
	  	end do

	  	!close that file and be done with it!
		close(90)

		!This is an adendum to write out the elastic strain energy associated with the adhesions
        write(filename,'(2a)') runname,'.adhesionenergy'
        open(02,file=filename,form='formatted',position='append')
        write(02,'(E25.12)'),adhesions%strain_energy

	end subroutine adh_write_out

	!========================================================================================================
	!THIS SUBROUTINE IS PASSED AN IMMERSED BOUNDARY TYPE. THEREFORE THERE ARE NO SCOPING ISSUES.
	!HOWEVER, IT MUST HAVE ACCESS TO THE CHAR 'RUNNAME' AND THE INT 'PLOTCOUNT'
	subroutine bnd_read_in(boundary)
		TYPE (im_bnd) :: boundary
		!CHARACTER(len = 30) :: filename
		INTEGER :: A
		REAL :: holder

		!write(filename,'(2a,i4.4)') runname,'.bnd.', 9999
		print*, 'incoming_boundary.txt'
		!THIS IS THE VERSION WHERE WE READ IN AS TEXT
  		open(40,file='incoming_boundary.txt',form='formatted')

  		read(40, '(i8)'), Mb !Now we know how many points are in the boundary


  		!Now we can begin allocating the component parts of the immersed boundary structure
  		!Locations and forces
  		allocate(boundary%im_bnd_loc(1:Mb,1:2),boundary%stretch_force(1:Mb,1:2),boundary%connect_force(1:Mb,1:2))
  		!Elastic parameter scalar arrays
  		allocate(boundary%stretch_kp(1:Mb),boundary%stretch_km(1:Mb),boundary%gammap(1:Mb),boundary%gammam(1:Mb))
        !Geometric vectors
        allocate(boundary%taup(1:Mb,1:2),boundary%taum(1:Mb,1:2),boundary%tau0(1:Mb,1:2),boundary%norm0(1:Mb,1:2))
        !Arc length measurements
        allocate(boundary%dsp(1:Mb),boundary%dsm(1:Mb),boundary%ds0(1:Mb))
        !Miscelaneous
        allocate(boundary%connect_buddy(1:Mb),boundary%flag(1:Mb))

  		do A=1,Mb
    		read(40, '(E25.12,E25.12)'), &
    		& boundary%im_bnd_loc(A,1),boundary%im_bnd_loc(A,2)

  		end do

  		close(40)

		!THIS IS THE VERSION WHERE WE WRITE OUT AS BINARY
!		open(40,file=filename,form='binary')
!		do A=1,Mb
!			read(30) boundary%im_bnd_loc(A,1), boundary%im_bnd_loc(A,2),&
!					& boundary%tether_pt_loc(A,1), boundary%tether_pt_loc(A,2),&
!					& boundary%stretch_force(A,1), boundary%stretch_force(A,2),&
!					& boundary%tether_force(A,1), boundary%tether_force(A,2), &
!                    & boundary%stretch_kp(A), boundary%tether_k(A), &
!					& holder
!					boundary%flag(A) = int(holder)
!		end do
!		close(40)
	end subroutine bnd_read_in


!========================================================================================================

!This is the routine to construct an immersed network type by reading in a file dumped out of MATLAB
!and DISTMESH. This is the tricky part. Test this fucker extremely well
!there are many parts to this routine
!IN THE FUTURE I MAY WANT TO BREAK THIS UP INTO MULTIPLE SUB-SUBROUTINES
	subroutine network_read_in(network)
		!The arguments passed in are an im_network type, and a string declaring the name of file to read from
		type(im_network) :: network
		!CHARACTER(len = 30) :: readfile
		!These are integers used to store the number of nodes, faces, unique links
		!and some dummy integers used for looping and counting
		INTEGER :: inc_nodes_M, inc_faces_M, inc_links_M
		!Several allocatable arrays of integers
		!inc_faces, inc_links will hold the actual faces and links when done
		INTEGER, DIMENSION(:,:), ALLOCATABLE :: inc_faces, inc_links
		!An allocatable array to read nodes into
		REAL, DIMENSION(:,:), ALLOCATABLE :: inc_nodes
		!we have a logical used to search for unique entries in all these lists
		LOGICAL :: found = .false.
		!These are integers used to store the number of nodes, faces, unique links
		!and some dummy integers used for looping and counting
		INTEGER :: cur_count, I, K
		!Several allocatable arrays of integers
		!holding _links and store_links will be used for the sorting algorithm
		INTEGER, DIMENSION(:,:), ALLOCATABLE :: holding_links, store_links




		!THIS IS THE VERSION WHERE WE READ IN FROM TEXT
	  	open(60,file='incoming_network.txt',form='formatted')
	  	!first we grab the 2 integers that tell us how many mesh nodes and how many faces
	  	!the immersed network will have
	  	read(60, '(i8)'), inc_nodes_M
	  	read(60, '(i8)'), inc_faces_M

!       print*, 'I NOW KNOW THE SIZE TO IMPORT!'
!       print*, inc_nodes_M, inc_faces_M

        !now we allocate arrays for the nodes and faces. we do not yet know how many links there
        !will be
  		allocate(inc_nodes(1:inc_nodes_M,1:2),inc_faces(1:inc_faces_M,1:3))

  		!Now we read in the locations of the nodes from the file
  		do I = 1,inc_nodes_M
  			read(60,'(E25.12,E25.12)'), inc_nodes(I,1), inc_nodes(I,2)
  		end do
!		print*, 'I have now read in the node locations'

		!Now we read in the array of integers that identifies which node triples
		!define each face of the immersed structure
		do I = 1,inc_faces_M
			read(60,'(i5,i5,i5)'), inc_faces(I,1), inc_faces(I,2), inc_faces(I,3)
		end do
		!Close the file we're reading out of
		close(60)
!		print*, 'I have now read in the face identifications!'
		!this is just a set of print statements for debugging purposes
!  		print*, inc_nodes_M, inc_faces_M
!  		do I = 1,inc_nodes_M
!			print*, inc_nodes(I,:)
!  		end do


		!NOW WE DO A SHIT LOAD OF WORK TO IDENTIFY AND ORGANIZE ALL OF THE
		!LINKS WITHIN THE IMMERSED NETWORK

		!Allocate a holding variable story all the NON-UNIQUE links of the network
		allocate(holding_links(1:3*inc_faces_M,1:2))

		!Write the integers defining all these NON-UNIQUE links into the holding array
		holding_links(1:inc_faces_M,:) = inc_faces(:,1:2)
		holding_links(inc_faces_M+1:2*inc_faces_M,:) = inc_faces(:,2:3)
		holding_links(2*inc_faces_M+1:3*inc_faces_M,1) = inc_faces(:,3)
		holding_links(2*inc_faces_M+1:3*inc_faces_M,2) = inc_faces(:,1)

		!Now loop through the holding array and sort each link so that its nodes are in
		!ascending order
		do I = 1,3*inc_faces_M
			if (holding_links(I,1).gt.holding_links(I,2))  then
				K = holding_links(I,1)
				holding_links(I,1) = holding_links(I,2)
				holding_links(I,2) = k
			end if
		end do

		!THIS IS JUST A PRINT STATEMENT FOR DEBUGGING PURPOSES
!		print*, 'NON-UNIQUE LINKS IN THE NETWORK'
!		do I = 1,3*inc_faces_M
!			print*, holding_links(I,:)
!		end do

		!I'm now allocating an array where I will place unique links one by one
		allocate(store_links(1:3*inc_faces_M,1:2))
		store_links = 0 !and setting it to zero

		!I assume that I have one unique link in the network (i.e. the first one)
		inc_links_M = 1
		!So I place it in my storage array
		store_links(1,:) = holding_links(1,:)

		!Now, I'm going to loop through my remaining NON-UNIQUE links and compare them to
		!All links with index K < I
		do I = 2,3*inc_faces_M
			!OK, I'M NOW EXAMINING LINK 'I'

			K = 1 !I start back from the beginning of non-unique links
			found = .false. !Assume that I have not found this link yet
			!And loop until I catch up to link 'I', or I find a repeat
			do while ((K.lt.I).and.(.not.found))
				!Check statement to see if 2 links under examination are the same
				if ((store_links(K,1).eq.holding_links(I,1)).and.&
					&(store_links(K,2).eq.holding_links(I,2))) then
					!If so, I've found myself
					found = .true.
				end if
				!And now I move on to examine the next link 'K'
				K = K+1
			end do

			!At this point, I have compared link 'I' to all links 'K' where K < I
			!If 'found' was never set to true, then link I is unique (so far)
			if (.not.found) then
				!So i increment the number of unique links
				inc_links_M = inc_links_M+1
				!PRINT STATMENTS FOR DEBUGGING
!				print*, 'FOUND UNIQUE LINK NUMBER:', inc_links_M
!				print*, holding_links(I,:)
				!And now I write this unique link into my storage array
				store_links(inc_links_M,:) = holding_links(I,:)
			end if
			!NOW I'M DONE WITH LINK 'I' AND I MOVE ON TO THE NEXT ONE
		end do

		!Having checked all links, i now know how many unique ones there are
		!so i allocate space for that
		allocate(inc_links(1:inc_links_M,1:2))
		!and write the array of unique links
		inc_links = store_links(1:inc_links_M,:)

		!THESE ARE PRINT STATEMENTS FOR DEBUGGING
!		print*, 'We found this many inc_links_M links', inc_links_M
!		do I = 1, inc_links_M
!			print*, inc_links(I,:)
!		end do

	!THIS REPRESENTS A LARGE MILESTONE IN THE SUBROUTINE. WE'VE NOW ORGANIZED THE INCOMING DATA
	!NOW WE MUST GO THROUGH BUILDING THE IM_NETWORK TYPE AND LINKING EVERYTHING PROPERLY
	!Now we allocate the arrays of nodes, faces and faces within the im_network type
		allocate(network%nodes(1:inc_nodes_M),network%faces(1:inc_faces_M),network%links(1:inc_links_M))

		!Assign the number of nodes and faces within the network type
		network%M_nodes = inc_nodes_M
		network%M_faces = inc_faces_M
		network%M_links = inc_links_M

!		print*, 'I HAVE NOW ALLOCATED SPACE FOR LINKS'


		!We will now loop through every NODE and fill in all its pertinent connectivity
		!as well as its location (this is needed to get everything started)

		do I = 1,inc_nodes_M !'I' is the index of which node we're taking care of

			!write the node location in (this is the easy part)
			network%nodes(I)%location(1,:) = inc_nodes(I,:)
!			network%nodes(I)%my_links = -1

			!NOW WE WILL MAKE THE NODE REFERENCE THE PROPER LINK
			!Count the number of links this node is in
			K = count(mask=inc_links.eq.I)
			!print statement for debugging
			!print*, 'Node', I, 'has', K, 'links coming out of it!'

			!Now allocate that many links within this network_node type
			allocate(network%nodes(I)%my_links(1:K))
			network%nodes(I)%link_count = K
			network%nodes(I)%my_links = 0

			!Now we will loop through ALL links and determine if node 'I' is contained in them
			cur_count = 1 !initilize cur_count at 1 because we're searching for this node's first link
			do K = 1,inc_links_M !'K' is the index that tracks which link we're examining

				!WE CHECK TO SEE IF THE K'TH LINK CONTAINS THE I'TH NODE
				if (count(mask=inc_links(K,:).eq.I).gt.0) then
!					print*, 'Node', I, 'is in the', K, 'link!'
!					print*, 'SEE:', inc_links(K,:), 'contains', I

					!if the answer is yes, i want to make the i'th node point to the k'th link
					network%nodes(I)%my_links(cur_count) = K
!					print*, 'Just altered node ', I
!					print*, network%nodes(I)%my_links
!					print*, 'Writing link', K, 'into spot ', cur_count, 'for node', I

					!and now we incriment cur_count to make sure we put the next reference in
					!the place in the 'my_links' array
					cur_count = cur_count + 1
!				    print*, 'AT STEP', I, 'THE FIRST NODE IS IN THESE LINKS'
!			        print*, network%nodes(1)%my_links
				end if
			end do



			!NOW WE WILL MAKE THE NODE REFERENCE THE PROPER FACES
			!Count the number of faces this node is in
			K = count(mask=inc_faces.eq.I)
			!print statement for debugging
!			print*, 'Node', I, 'is in', K, 'faces!!'
			!Now allocate that many faces within this network_node type
			allocate(network%nodes(I)%my_faces(1:K))
			network%nodes(I)%face_count = K
			!Now we will loop through ALL faces and determine if node 'I' is contained in them
			cur_count = 1 !initilize cur_count at 1 because we're searching for this node's first face
			do K = 1,inc_faces_M !'K' is the index that tracks which face we're examining
				!WE CHECK TO SEE IF THE K'TH LINK CONTAINS THE I'TH NODE
				if (count(mask=inc_faces(K,:).eq.I).gt.0) then
					!if the answer is yes, i want to make the i'th node point to the k'th link
					network%nodes(I)%my_faces(cur_count) = K
					!and now we incriment cur_count to make sure we put the next reference in
					!the place in the 'my_faces' array
					cur_count = cur_count + 1
				end if
			end do
		end do
		!THAT COMPLETES THE CREATION OF ALL THE NODES IN OUR NETWORK!


		!We will now loop through every FACE and fill in all its pertinent connectivity
		do I = 1,inc_faces_M !'I' is the index keeping track of which face we're dealing with
			!Write integer triples defining each face into the im_network type
			network%faces(I)%my_nodes(:) = inc_faces(I,:)
			!That takes care of all the node connectivity
			!now we need to link the face to its links

			!we will loop through ever LINK and find the ones that correspond to this FACE
			cur_count = 1 !we set this because we're looking for the first 'correct' link
			do K = 1,inc_links_M

			!this is the check statement we want the first AND second entries in link 'k'
			!to be in face 'I'
				if ((count(mask=inc_faces(I,:).eq.inc_links(K,1)) +&
					&count(mask=inc_faces(I,:).eq.inc_links(K,2))).eq.2) then
					!if that's true, then this link is in face 'I'
					network%faces(I)%my_links(cur_count) = K
					!now we proceed to look for the next link in this face
					cur_count = cur_count + 1
				end if
			end do
			!At this point, we've examined all possible links
		end do
		!THAT COMPLETES THE CREATION OF ALL THE FACES IN OUR NETWORK!


		!We will now loop through every link and fill in all its pertinent connectivity
		do I = 1,inc_links_M !'I' is the index keeping track of which link we're dealing with
			!Write integer doubles defining each link into the im_network type
			network%links(I)%my_nodes(:) = inc_links(I,:)
			!That takes care of all the node connectivity
			!now we need to link the link to its faces

			!INITIALIZE THIS VALUE TO ZERO. IF IT REMAINS ZERO AFTER EVERTYHING IS DONE
			!THAT INDICATES THAT THIS LINK IS AN 'EXTERIOR' LINK
			network%links(I)%my_faces = 0

			!we will loop through ever FACE and find the ones that correspond to this LINK
			cur_count = 1 !we set this because we're looking for the first 'correct' FACE
			do K = 1,inc_faces_M

			!this is the check statement we want the first AND second entries in link 'I'
			!to be in face 'K'
				if ((count(mask=inc_faces(K,:).eq.inc_links(I,1)) +&
					&count(mask=inc_faces(K,:).eq.inc_links(I,2))).eq.2) then
					!if that's true, then this face is bounded by link 'I'
					network%links(I)%my_faces(cur_count) = K
					!now we proceed to look for the next link in this face
					cur_count = cur_count + 1
				end if
			end do
			!At this point, we've examined all possible faces
		end do

		!Lets deallocate this shit just for good measure
		deallocate(holding_links,store_links)
	!THAT COMPLETES THE CREATION OF ALL THE LINKS IN OUR NETWORK!
	end subroutine network_read_in






end module io_mod
