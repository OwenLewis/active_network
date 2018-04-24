! ============================================================================
! Name        : geometry_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Initialization and update module for all goemetric information
!               This will set and update the boundary, tethers, and material
!				properties. This includes the reference boundary config. Starting
!				boundary config, etc.
!THIS VERSION OF THE MODULE IS FOR THE CAPPED TUBE GEOMETRY
! ============================================================================

module geometry_mod
USE global_param_mod
USE useful_functions_mod
implicit none
private


	REAL, DIMENSION(:), ALLOCATABLE :: tens_wave, link_body_coord, node_body_coord, adh_wave, scaled_amp

public configure_network

!========================================================================================================
contains



!========================================================================================================
!This subroutine will configure the immersed network according to simulation time. THIS ROUTINE IS ENTIRELY
!SIMULATION AND CONTEXT DEPENDENT. I CAN CHANGE IT WHENEVER I WANT.


subroutine configure_network(network,adhesions,time)
	!We input an 'im_network' type and the current time of the simulation
	TYPE (im_network) :: network
	TYPE (adh_complex) :: adhesions
	REAL :: time, left, right, pulse, modular
	!We also have a couple of reals used for calculation
	!and an integer used for some shit
	INTEGER :: I
	!'tens_wave' and 'body_coord' are stored at the module level.
	!We see if they are allocated yet. If not, do that and populate the body coordinate
	if (.not.allocated(tens_wave)) then
		!The wave should be the same size as number of links
		I = network%M_links
		allocate(tens_wave(1:I),link_body_coord(1:I),scaled_amp(1:I))
		I = adhesions%M_adh
		allocate(node_body_coord(1:I),adh_wave(1:I))

		!Find the left & right bounds of the network (@ reference)
 		left = minval(network%links(:)%ref_center(1,1))
 		right = maxval(network%links(:)%ref_center(1,1))
		!Now we make a body coordinate for each link that goes from zero to one
		link_body_coord = (network%links(:)%ref_center(1,1) - left)/(right - left)
!		print*, 'Body coordinate used to input wave of contraction'
!		print*, body_coord
!		print*,'They run from ', maxval(body_coord), 'to', minval(body_coord)

        left = minval(adhesions%ref_location(:,1))
        right = maxval(adhesions%ref_location(:,1))
        node_body_coord = (adhesions%ref_location(:,1) - left)/(right - left)


        !Here is where I change the amplitude of the active contractive wave locally for each link
        select case(isotropic)

            !This is the original version where I do not adjust the amplitude of contraction according to
            !the orientation of the link
            case(1)
                print*, 'We are not scaling contraction: Isotropic'
                scaled_amp = Amp

            !This is the version where I simply adjusted the amplitude according to the orientation of the link
            !I.E. vertical links will contract (actively) more than horizontal ones
            case(0)
                print*, 'We are scaling contraction: Anisotropic'
                scaled_amp = Amp*abs(network%links(:)%rayup(1,2)/network%links(:)%link_length)
        end select

        !Now we will adjust the contraction wave to be spacially heterogeneous or not

        select case(homogeneous)

            !This is the version where the amplitude is indipendent of spacial coordinate
            case(1)
                print*, 'We are not spacialy adjusting contraction strength: Homogeneous'

            !This is the version where we scale contraction by a heavyside function
            case(0)
                print*, 'We are spacially adjusting contraction strength: Heavyside heterogeneity'
                scaled_amp = scaled_amp*(merge(1.,head_multiple,link_body_coord.lt.0.8))

            !This is the version where we scale contraction by a linearly decreasing function
            case(2)
                print*, 'We are spacially adjusting contraction strength: Linear heterogeneity'
                scaled_amp = scaled_amp*(1. - link_body_coord*(1.-head_multiple))
        end select


        !This is just some print statements for the log file so that I can see what I've done
        do I = 1,network%M_links
            if (isnan(scaled_amp(I))) then
                print*, 'Something went wrong. Divided by zero'
            end if
            print*, 'the Amplitude of contraction has been scaled by', scaled_amp(I)/Amp
        end do
	end if




	!The tension wave can be chose as one of several wave forms

	select case(contractwave)
	    !This is the case where we have a spacially uniform contraction
	    case(0)
	        tens_wave = scaled_amp

        !This is the case where we have our standard traveling wave of contraction
        case(1)
            tens_wave = scaled_amp*(cos(con_wave_num*link_body_coord - freq*time) + 1)/2.

        !This is the case where we have a uniform, oscillatory contraction
        case(2)
            tens_wave = scaled_amp*(sin(freq*time) + 1)/2.
    end select


    !Here we can do an increasing exponential to eliminate the "jerk" associated with instantly turning on contraction.
    if (time.le.1.0) tens_wave = tens_wave*(1.-exp(-time/contau))


    !Finally, we construct the adhesion wave
    select case(adhwave)
        case(1)
            !The adhesion wave a phase wave plus a constant 'base_adh_slip' to avoid zero adhesion
            adh_wave = adhesion_slip*(cos(adh_wave_num*node_body_coord - freq*time + adh_offset) + 1)*0.5 + base_adh_slip

        case(0)
            !The adhesion wave is a slightly more complex standing wave fo the body coordinate
            adh_wave = adhesion_slip*(cos(Pi*node_body_coord)*cos(freq*time + adh_offset) + 1)*0.5 + base_adh_slip
    end select


    !Put that shit in the proper place
	network%links(:)%active_tens = tens_wave
	adhesions%adh_slip = adh_wave

end subroutine configure_network



end module geometry_mod
