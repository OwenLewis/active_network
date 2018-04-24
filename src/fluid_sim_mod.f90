! ============================================================================
! Name        : fluid_sim_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing the variables defining a current state of
!				an immersed boundary simulation as well as the routines to
!				evolve the simulation
! ============================================================================

module fluid_sim_mod

	USE global_param_mod
	USE io_mod
	USE im_bnd_intern_mod
	USE im_bnd_opers_mod
	USE cortex_intern_mod
	USE cortex_opers_mod
	USE adhesion_opers_mod
	USE stokes_2D_fft_mod
	USE geometry_mod

	implicit none

	PRIVATE

	REAL, DIMENSION(:,:), ALLOCATABLE :: pressure, Xc, Yc, eul_torque1, eul_torque2
	REAL, DIMENSION(:,:), ALLOCATABLE :: vel_horiz, Fh, Fhfromnet, Fhfrombnd, FhSubstrate
	REAL, DIMENSION(:,:), ALLOCATABLE :: vel_vert, Fv, Fvfromnet, Fvfrombnd, FvSubstrate
	REAL, DIMENSION(:,:), ALLOCATABLE :: bnd_force, vel_bnd, loc_bnd_new
	REAL, DIMENSION(:), ALLOCATABLE :: mem_torque1, mem_torque2
	REAL, DIMENSION(:,:), ALLOCATABLE :: network_force, network_vel, network_interp_vel, adh_correction, adh_vel
	REAL, DIMENSION(:), ALLOCATABLE :: adh_weight, net_torque1, net_torque2
	REAL :: scalar, foo
	REAL :: tol
	INTEGER :: tot_markers, tot_nodes, tot_links, tot_faces, tot_adh, plots = 0
	TYPE (im_bnd) :: current_bnd
	TYPE (im_network), SAVE :: immersed_network
	TYPE (adh_complex) :: network_adhesions


	public initialize_fluid, im_bnd_network_sim_step, write_everything_out
	public pressure, vel_horiz, vel_vert, Fh, Fv

!=========================================================================================
	contains

!THIS IS THE BIG ONE THIS IS THE ROUTINE THAT TAKES A FULL IMMERSED BOUNDARY STEP
	subroutine im_bnd_network_sim_step(time)
		REAL :: TIME
		!THIS IS JUST A DUMBY INT USED FOR INDEXING LOOPS
		INTEGER :: I
		REAL, DIMENSION(1,1:2):: update_vec, test_vec, u_sub = 0.
        REAL, DIMENSION(1:Mb,1:2) :: force_vec

		!AND HERE WE GO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!WE RECONFIGURE THE INTERNAL NETWORK ACCORDING TO CURRENT SIMULATION TIME
		call configure_network(immersed_network,network_adhesions,time)

		!NOW WE UPDATE THE QUANTITIES ASSOCIATED WITH THE NEW
		!CONFIGURATION
		call im_bnd_populate(current_bnd)
		call network_populate(immersed_network)
		call link_force_calc(immersed_network,current_bnd)
		call adh_populate(immersed_network,network_adhesions)



		!WE SUM ALL FORCES ON THE NETWORK SO THAT THEY CAN BE SPREAD TO GRID
		network_force(:,1) = immersed_network%nodes(:)%elastic_force(1,1) + immersed_network%nodes(:)%connect_force(1,1) &
						&+ immersed_network%nodes(:)%adh_force(1,1)
		network_force(:,2) = immersed_network%nodes(:)%elastic_force(1,2) + immersed_network%nodes(:)%connect_force(1,2) &
							&+ immersed_network%nodes(:)%adh_force(1,2)
!		print*, 'Elastic force on network', &
!		&sum(immersed_network%nodes(:)%elastic_force(1,1)*immersed_network%nodes(:)%ref_area), &
!		&sum(immersed_network%nodes(:)%elastic_force(1,2)*immersed_network%nodes(:)%ref_area)
!		print*, 'Maximum Elastic force on network', &
!		&maxval(abs(immersed_network%nodes(:)%elastic_force(1,1))), &
!		&maxval(abs(immersed_network%nodes(:)%elastic_force(1,2)))
!
!		print*, 'Connection force on network', &
!		&sum(immersed_network%nodes(:)%connect_force(1,1)*immersed_network%nodes(:)%ref_area), &
!		&sum(immersed_network%nodes(:)%connect_force(1,2)*immersed_network%nodes(:)%ref_area)

!		print*, 'Maximum Connection force on network', &
!		&maxval(abs(immersed_network%nodes(:)%connect_force(1,1))), &
!		&maxval(abs(immersed_network%nodes(:)%connect_force(1,2)))
!
!		print*, 'Adhesion force on network', &
!		&sum(immersed_network%nodes(:)%adh_force(1,1)*immersed_network%nodes(:)%ref_area), &
!		&sum(immersed_network%nodes(:)%adh_force(1,2)*immersed_network%nodes(:)%ref_area)

!		print*, 'Maximum Adhesion force on network', &
!		&maxval(abs(immersed_network%nodes(:)%adh_force(1,1))), &
!		&maxval(abs(immersed_network%nodes(:)%adh_force(1,2)))


!		print*, 'Total force on network', sum(network_force(:,1)*immersed_network%nodes(:)%ref_area), &
!										& sum(network_force(:,1)*immersed_network%nodes(:)%ref_area)
!
!		print*, 'Connection force on boundary', &
!		&sum(current_bnd%connect_force(:,1)*ref_bnd%ds0), sum(current_bnd%connect_force(:,2)*ref_bnd%ds0)




		!Now we sum all the forces on the boundary so that they can be spread to the grid
		bnd_force = current_bnd%stretch_force + current_bnd%connect_force



		!Now we spread both the boundary and network forces to the fluid domain
		call grid_spread(current_bnd%im_bnd_loc,bnd_force,Fhfrombnd,Fvfrombnd)
		call network_grid_spread(immersed_network,network_force,Fhfromnet,Fvfromnet)

		!Add those togeter for the total force acting on the fluid
		Fh = Fhfrombnd + Fhfromnet
		Fv = Fvfrombnd + Fvfromnet


		if ( ( abs(sum(Fh(1:Nx,1:Ny))*hx*hy).GT.tol ) .or. ( abs(sum(Fv(1:Nx,1:Ny))*hx*hy).GT.tol) ) then
        print*, 'SHIT WE HAVE A PROBLEM!!'
		print*, 'Total force in Eulerian Domain', sum(Fh(1:Nx,1:Ny))*hx*hy, sum(Fv(1:Nx,1:Ny))*hx*hy
		print*, 'At time: ', time
        end if

        !Now we are going to calculate the torques being applied to the fluid.
        !This is for debugging purposes.
!		eul_torque1 = Fh*Yc
!
!		eul_torque2 = Fv*Xc

        !PRINT AN ERROR STATEMENT IF THAT SHIT IS TOO BIG
!        if (abs((sum(eul_torque1-eul_torque2)*hx*hy)).gt.tol) then
!        print*, 'SHIT WE HAVE A PROBLEM'
!		print*, 'Total torque on Eulerian Domain', sum(eul_torque1-eul_torque2)*hx*hy
!        PRINT*, 'At time: ', time
!        end if


!        !Now we are going to calculate the torques on the membrane from internal forces.
!        !This is for debugging purposes.
!        mem_torque1 = current_bnd%stretch_force(:,1)*current_bnd%im_bnd_loc(:,2)
!        mem_torque2 = current_bnd%stretch_force(:,2)*current_bnd%im_bnd_loc(:,1)
!        mem_torque1 = mem_torque1 - mem_torque2
!
!        !PRINT AN ERROR STATEMENT IF THAT SHIT IS TOO BIG
!        if (abs(sum(mem_torque1*ref_bnd%ds0)).gt.tol) then
!        print*, 'SHIT WE HAVE A PROBLEM'
!        print*, 'Total torque on Lagrangian Membrain (due to internal forces)', sum(mem_torque1*ref_bnd%ds0)
!        print*, 'At time: ', time
!        end if
!
!        !Now we are going to calculate the torques on the membrane from connection to network.
!        !This is for debugging purposes.
!
!        mem_torque1 = current_bnd%connect_force(:,1)*current_bnd%im_bnd_loc(:,2)
!        mem_torque2 = current_bnd%connect_force(:,2)*current_bnd%im_bnd_loc(:,1)
!        mem_torque1 = mem_torque1 - mem_torque2
!
!
!        !PRINT AN ERROR STATEMENT IF THAT SHIT IS TOO BIG
!        if (sum(mem_torque1*ref_bnd%ds0).gt.tol) then
!        print*, 'SHIT WE HAVE A PROBLEM'
!        print*, 'Total torque on Lagrangian Membrain (due to connection force)', sum(mem_torque1*ref_bnd%ds0)
!        print*, 'At time: ', time
!        end if

        !Now we are going to calculate the torques on the network from connection forces.
        !This is for debugging purposes.

!        net_torque1 = immersed_network%nodes(:)%connect_force(1,1)*immersed_network%nodes(:)%location(1,2)
!        net_torque2 = immersed_network%nodes(:)%connect_force(1,2)*immersed_network%nodes(:)%location(1,1)
!        net_torque1 = net_torque1 - net_torque2
!
!        call scalar_network_integrate(foo,net_torque1,immersed_network)
!
!
!        !PRINT AN ERROR STATEMENT IF THAT SHIT IS TOO BIG
!        if (abs(foo).gt.tol) then
!        print*, 'SHIT WE HAVE A PROBLEM'
!        print*, 'Total torque on Lagrangian Network (due to connection forces)', foo
!        print*, 'At time :', time
!        end if
!
!
!        !Now we are going to calculate the torques on the entwork from internal forces.
!        !This is for debugging purposes.
!
!        net_torque1 = immersed_network%nodes(:)%elastic_force(1,1)*immersed_network%nodes(:)%location(1,2)
!        net_torque2 = immersed_network%nodes(:)%elastic_force(1,2)*immersed_network%nodes(:)%location(1,1)
!        net_torque1 = net_torque1 - net_torque2
!
!        call scalar_network_integrate(foo,net_torque1,immersed_network)
!
!        !PRINT AN ERROR STATEMENT IF THAT SHIT IS TOO BIG
!        if (abs(foo) .gt. tol) then
!        print*, 'SHIT WE HAVE A PROBLEM'
!        print*, 'Total torque on Lagrangian Network (due to internal elastic forces)', foo
!        print*, 'At time :', time
!        end if
!
!
!        !Now we are going to calculate the torques on the entwork from adhesive forces.
!        !This is for debugging purposes.
!        net_torque1 = immersed_network%nodes(:)%adh_force(1,1)*immersed_network%nodes(:)%location(1,2)
!        net_torque2 = immersed_network%nodes(:)%adh_force(1,2)*immersed_network%nodes(:)%location(1,1)
!        net_torque1 = net_torque1 - net_torque2
!
!        call scalar_network_integrate(foo,net_torque1,immersed_network)

        !PRINT AN ERROR STATEMENT IF THAT SHIT IS TOO BIG
!        if (abs(foo) .gt.tol) then
!        print*, 'SHIT WE HAVE A PROBLEM'
!        print*, 'Total torque on Lagrangian Network (due to adhesive forces)', foo
!        print*, 'At time: ', time
!        end if


!AND SO ENDS THE HUGE ASS STRING OF DEBUGGING PRINT STATEMENTS

		!NOW SOLVE THE FLUID EQUATION
		call stokes_solve_2d_fft(Fh,Fv,vel_horiz,vel_vert,pressure)

		!We interpolate the fluid velocity to the immersed boundary and network structures
		call boundary_interpolate(vel_horiz,vel_vert,current_bnd%im_bnd_loc,vel_bnd)
		call network_interpolate(vel_horiz,vel_vert,immersed_network,network_interp_vel)

!		print*, 'Maximum boundary velocity:', maxval(vel_bnd(:,1)), maxval(vel_bnd(:,2))
!		print*, 'Maximum velocity interpolated to network', maxval(network_interp_vel(:,1)),maxval(network_interp_vel(:,2))

		!Calculate the actual translocation speed of the immersed network
		network_vel(:,1) = network_force(:,1)/immersed_network%nodes(:)%node_drag + network_interp_vel(:,1)
                network_vel(:,2) = network_force(:,2)/immersed_network%nodes(:)%node_drag + network_interp_vel(:,2)

!		print*, 'Maximum network velocity', maxval(network_vel(:,1)),maxval(network_vel(:,2))



		!NOW HERE IS THE TRICKY PART. THIS IS WHERE WE CALULATE U_SUB TO ENSURE THAT FORCES REMAINED BALANCED AT THE NEXT TIME STEP
        if (translate.eq.1) then

!           print*, 'WE ARE DOING THAT WEIRD SUBSTRATE TRANSLATION THING' 
            !1 by 1 we will add the terms to the correction term used to calculate usub
            !first is ahesion force times (k^2/xsi + k^2/eta)
            adh_correction(:,1) = immersed_network%nodes(:)%adh_force(1,1)*&
                &(network_adhesions%adh_stiff(:)/immersed_network%nodes(:)%node_drag + &
                &network_adhesions%adh_stiff(:)/network_adhesions%adh_slip(:))

            adh_correction(:,2) = immersed_network%nodes(:)%adh_force(1,2)*&
                &(network_adhesions%adh_stiff(:)/immersed_network%nodes(:)%node_drag + &
                &network_adhesions%adh_stiff(:)/network_adhesions%adh_slip(:))

            !The second term looks like k times the interplated fluid velocity
            adh_correction(:,1) = adh_correction(:,1) + network_adhesions%adh_stiff(:)*network_interp_vel(:,1)
            adh_correction(:,2) = adh_correction(:,2) + network_adhesions%adh_stiff(:)*network_interp_vel(:,2)

            !The third term looks like k times the elastic and connection forces already on the network
         adh_correction(:,1) = adh_correction(:,1) + network_adhesions%adh_stiff(:)*(immersed_network%nodes(:)%elastic_force(1,1)&
                                &+ immersed_network%nodes(:)%connect_force(1,1))/immersed_network%nodes(:)%node_drag
         adh_correction(:,2) = adh_correction(:,2) + network_adhesions%adh_stiff(:)*(immersed_network%nodes(:)%elastic_force(1,2)&
                                &+ immersed_network%nodes(:)%connect_force(1,2))/immersed_network%nodes(:)%node_drag

            !We also have this weighting term which looks like the stiffnesses of the adhesive links
            adh_weight = network_adhesions%adh_stiff

            !Integrate the weighting term and the correction term around the immersed network
            call vec_network_integrate(u_sub,adh_correction,immersed_network)
            call scalar_network_integrate(scalar,adh_weight,immersed_network)

            !the speed of the substrate (u_sub) is the ratio of these two terms (from above integrals)
            u_sub = u_sub/scalar

    !		print*, 'Substrate velocity is', u_sub

            !We integrate the markers with the substrate velocity
    !		do I = 1,tot_markers
    !			markers(I,:) = markers(I,:) + dt*u_sub(1,:)
    !		end do

    !		print*, 'Center of mass is:', markers

            !Now we adjust the fluid velocity, membrane velocity, and network velocity according to "U_sub"
            vel_horiz = vel_horiz - u_sub(1,1)
            vel_vert = vel_vert - u_sub(1,2)

            network_vel(:,1) = network_vel(:,1) - u_sub(1,1)
            network_vel(:,2) = network_vel(:,2) - u_sub(1,2)

            vel_bnd(:,1) = vel_bnd(:,1) - u_sub(1,1)
            vel_bnd(:,2) = vel_bnd(:,2) - u_sub(1,2)
        end if


		!Now we calculate the velocity of the adhesions
		adh_vel(:,1) = -immersed_network%nodes(:)%adh_force(1,1)/network_adhesions%adh_slip(:)
		adh_vel(:,2) = -immersed_network%nodes(:)%adh_force(1,2)/network_adhesions%adh_slip(:)

		!We integrate in time the boundary (with fluid velocity)
		call boundary_time_integrate(current_bnd%im_bnd_loc,vel_bnd,loc_bnd_new)
		current_bnd%im_bnd_loc = loc_bnd_new

		!We also integrate the network and the adhesions with their own velocities
		call network_time_integrate(immersed_network,network_vel)
		call adh_time_integrate(network_adhesions,adh_vel)

		!THAT'S IT, WE'RE ALL GOOD
		!READY FOR NEXT STEP

	end subroutine im_bnd_network_sim_step


!========================================================================================================

!========================================================================================================
	!This is the subroutine that initializes the fluid simulation
	subroutine initialize_fluid
		INTEGER :: I
		character(len = 20) :: filename = 'input_param.txt'

		!READ THE SOFT PARAMETERS IN FROM EXTERNAL FILE
		OPEN(15,file=filename,form='formatted')
		READ(15, nml=soft_params)
		CLOSE(15)

        hx = Xmax/Nx
        hy = Ymax/Ny

        allocate(Xc(1:Nx,1:Ny),Yc(1:Nx,1:Ny))
        allocate(pressure(1:Nx,1:Ny),vel_horiz(1:Nx,1:Ny),vel_vert(1:Nx,1:Ny))
        allocate(Fh(1:Nx,1:Ny),Fv(1:Nx,1:Ny),Fhfrombnd(1:Nx,1:Ny),Fvfrombnd(1:Nx,1:Ny))
        allocate(Fhfromnet(1:Nx,1:Ny),Fvfromnet(1:Nx,1:Ny),FhSubstrate(1:Nx,1:Ny),FvSubstrate(1:Nx,1:Ny))

		!construct the various grids associated with the Eulerian domain
		do I = 1,Nx
			Xc(I,:) = (I - 0.5)*hx
		end do

		do I = 1,Ny
			Yc(:,I) = (I - 0.5)*hy
		end do


		!zero out the velocity fields
		vel_horiz = 0.
		vel_vert = 0.
		pressure = 0.

		!Zero out the substrate stresses
		FhSubstrate = 0.
		FvSubstrate = 0.

		!This sets the reference boundary configuration
		call bnd_read_in(ref_bnd)
		call im_bnd_initialize(ref_bnd)
		call im_bnd_populate(ref_bnd)

		allocate(bnd_force(1:Mb,1:2),vel_bnd(1:Mb,1:2),loc_bnd_new(1:Mb,1:2))
		allocate(mem_torque1(1:Mb),mem_torque2(1:Mb))

		!NOW WE SET THE CURRENT BOUNDARY TO THE STARTING POSITION
		current_bnd = ref_bnd

        !NOW WE RIGHT OUT THE REFERENCE CONFIGURATION
        call bnd_write_out(current_bnd,9999)

		call network_read_in(immersed_network)
		call network_initialize(immersed_network)
		call network_populate(immersed_network)
		call linkup(immersed_network,current_bnd)

		call adh_initialize(immersed_network,network_adhesions)
		call adh_populate(immersed_network,network_adhesions)


        !WRITE ALL THE PARAMETERS OUT TO EXTERNAL FILE
        write(filename,'(2a)') runname,'.param'
        print*, filename
        open(16,file=filename)
        WRITE(16,nml=soft_params)
        CLOSE(16)

		I = immersed_network%M_nodes

		allocate(network_force(1:I,1:2),network_vel(1:I,1:2),network_interp_vel(1:I,1:2))

		!These torque variables are used for debugging purposes
		allocate(net_torque1(1:I),net_torque2(1:I))

		I = network_adhesions%M_adh
		allocate(adh_correction(1:I,1:2),adh_vel(1:I,1:2),adh_weight(1:I))

        tol = 10.**(-11.)
        print*, 'Tolerance is set to', tol

	end subroutine initialize_fluid

!===============================================================================================================

	subroutine write_everything_out
	    INTEGER :: I


		call fluid_write_out(vel_horiz,vel_vert,pressure,plots)
		call bnd_write_out(current_bnd,plots)
		call network_write_out(immersed_network,plots)
		call adh_write_out(network_adhesions,plots)
		call lag_force_write_out(immersed_network,plots)

        !Now we're going to spread the spread the adhesion forces to the 'substrate' and write them out
        do I = 1,network_adhesions%M_adh
            adh_vel(I,:) = -immersed_network%nodes(network_adhesions%network_buddy(I))%adh_force(1,:)
        end do

        call adh_grid_spread(network_adhesions,immersed_network,adh_vel,FhSubstrate,FvSubstrate)

        call substrate_write_out(FhSubstrate,FvSubstrate,plots)
		plots = plots+1

	end subroutine write_everything_out



end module fluid_sim_mod
