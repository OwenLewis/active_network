! ============================================================================
! Name        : im_bnd_intern.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing immersed boundary data and internal functions to
!				the immersed boundary (such as force calculations).
! ============================================================================
module im_bnd_intern_mod
USE global_param_mod
implicit none
private
public im_bnd_populate, scalar_bnd_integrate, vec_bnd_integrate, im_bnd_initialize


contains

!=========================================================================================================


    subroutine im_bnd_initialize(boundary)
    type(im_bnd) :: boundary
    INTEGER :: I
    REAL , DIMENSION(1:Mb) :: linkmid
    REAL :: left, right
    LOGICAL, DIMENSION(1:Mb) :: rightleftlog


    !We are going to calculate the middle of each link in the immersed boundary
    linkmid(1:Mb-1) = (boundary%im_bnd_loc(1:Mb-1,1) + boundary%im_bnd_loc(2:Mb,1))/2.
    linkmid(Mb) =  (boundary%im_bnd_loc(Mb,1) + boundary%im_bnd_loc(1,1))/2.

    !This is all code that I use to make various different spacially dependent elastic parameters for the boundary

    !THIS ONE IS A STEP FUNCTION
    !Now create a logical array that says if each point is in the left or right region of cell
    rightleftlog = linkmid.lt.0.8

    !turn that logical into values by which the membrane elasticity will be scaled
    linkmid = merge(1.,head_multiple,rightleftlog)

    !Now actually calculate the link-by-link elastic constant
    !we have a few cases

    select case (homogeneous)
        !The first case is where our cell is homogeneous
        case(1)
            !we simply make the boundary elastic parameter a constant
            boundary%stretch_kp = base_stretch_k
            boundary%stretch_km = base_stretch_k


        !The second case is a heterogeneous cell
        case(0)
            !we make the boundary elastic parameter a step funciton
            boundary%stretch_kp(1:Mb) = base_stretch_k*linkmid(1:Mb)
            boundary%stretch_km(2:Mb) = base_stretch_k*linkmid(1:Mb-1)
            boundary%stretch_km(1) = base_stretch_k*linkmid(Mb)

        !The third case is a linearly hetereogeneous cell
        case(2)
            !Find the right and left edges of the cell
            right = maxval(boundary%im_bnd_loc(:,1))
            left = minval(boundary%im_bnd_loc(:,1))

            !Scale the link coordinate to make it run from zero to one
            linkmid = (linkmid - left)/(right - left)

            !Now set the linearly decreasing elastic parameters
            boundary%stretch_kp = base_stretch_k*(1 - linkmid/2.)
            boundary%stretch_km(2:Mb) = boundary%stretch_kp(1:Mb-1)
            boundary%stretch_km(1) = boundary%stretch_kp(Mb)
    end select



    !Loading pressure
    boundary%gammap = load_gamma
    boundary%gammam = load_gamma

end subroutine im_bnd_initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!THIS INITIALIZATION ROUTINE POUPLATES MOST OF THE VARIOUS QUANTITIES ASSOCIATED WITH THE
!IMMERSED BOUNDARY TYPE. IN ORDER TO RUN PROPERLY, THE LOCATION OF THE BOUNDARY, THE LOCATION
!OF THE TETHER POINTS, CONSTANTS ASSOCIATED WITH THE VARIOUS FORCES, AND TARGET 'LENGTHS' AND 'CURVATURES'
!MUST BE POPULATED ALREADY
!=========================================================================================================
subroutine im_bnd_populate(boundary)
!The only input is the immersed boundary which is of type 'im_bnd'
!IMPORTANT: THIS ROUTINE ASSUMES THAT THE BOUNDARY POINTS AND TETHER POINTS HAVE ALREADY BEEN UPDATED
!IT SIMPLY UPDATES ALL OF THE OTHER QUANTITIES ASSOCIATED WITH THE IMMERSED BOUNDARY
type (im_bnd) :: boundary

!Calculate forward difference vector
boundary%taup(1:Mb-1,:) = boundary%im_bnd_loc(2:Mb,:) - boundary%im_bnd_loc(1:Mb-1,:)
boundary%taup(Mb,:) = boundary%im_bnd_loc(1,:) - boundary%im_bnd_loc(Mb,:)
!Calculate backward difference vector
boundary%taum(2:Mb,:) = boundary%im_bnd_loc(2:Mb,:) - boundary%im_bnd_loc(1:Mb-1,:)
boundary%taum(1,:) = boundary%im_bnd_loc(1,:) - boundary%im_bnd_loc(Mb,:)
!Average them to approximate centered difference vector
boundary%tau0 = (boundary%taup + boundary%taum)/2.

!Now we calculate the normal by rotating the tangent
!This code assumes that the boundary is parametrized counter clockwise
!And generates the inward pointing normal vector
boundary%norm0(:,1) = -boundary%tau0(:,2)
boundary%norm0(:,2) = boundary%tau0(:,1)

!Now we calculate the difference lengths used in calculating deriviatives
boundary%dsp = sqrt(boundary%taup(:,1)**2. + boundary%taup(:,2)**2.)
boundary%dsm = sqrt(boundary%taum(:,1)**2. + boundary%taum(:,2)**2.)
!boundary%ds0 = sqrt(boundary%tau0(:,1)**2. + boundary%tau0(:,2)**2)
boundary%ds0 = (boundary%dsp + boundary%dsm)/2.


!print*, 'Before calculating, stretch force in boundary is: '
!print*, boundary%stretch_force(:,1)
!print*, boundary%stretch_force(:,2)

call calc_stretch_force(boundary)

end subroutine im_bnd_populate




!=============================================================================================
!THIS SUBROUTINE UPDATES THE INTERNAL FORCES DUE TO STRETCHING ASSOCIATED WITH A BOUNDARY
!CONFIGURATION. IN ORDER FOR THIS TO WORK, THE IMMERSED BOUNDARY POINTS MUST ALREADY BE UPDATED
!SIMILARLY, THE PARAMETER 'STRETCH_K' SHOULD BE UPDATED (IF DESIRED) BEFORE CALLING THIS
subroutine calc_stretch_force(boundary)
	!SIMPLY PASS THE IMMERSED BOUNDARY TYPE THAT THIS WILL WORK ON
	TYPE (im_bnd) :: boundary
	!INTERNAL VARIABLES FOR CALCULATES
	REAL, DIMENSION(1:Mb,1:2) :: deriv_p, deriv_m, oriented_p, oriented_m
	REAL, DIMENSION(1:Mb) :: tens_p, tens_m

	!First we calculate the first derivative of the current IM
    !with respect to REFERENCE ARC LENGTH
    deriv_p(:,1) = boundary%taup(:,1)/ref_bnd%dsp
	deriv_p(:,2) = boundary%taup(:,2)/ref_bnd%dsp
	deriv_m(:,1) = boundary%taum(:,1)/ref_bnd%dsm
	deriv_m(:,2) = boundary%taum(:,2)/ref_bnd%dsm

	!We now calculate the stretching tension associated
	!with the deformed state of the boundary
	tens_p = boundary%stretch_kp*(sqrt(deriv_p(:,1)**2. + deriv_p(:,2)**2.) - 1.) + boundary%gammap
	tens_m = boundary%stretch_km*(sqrt(deriv_m(:,1)**2. + deriv_m(:,2)**2.) - 1.) + boundary%gammam


	!Now we multiply the tension by the unit tangent vector to the
	!current boundary configuration
	oriented_p(:,1) = tens_p*deriv_p(:,1)/sqrt(deriv_p(:,1)**2. + deriv_p(:,2)**2.)
	oriented_p(:,2) = tens_p*deriv_p(:,2)/sqrt(deriv_p(:,1)**2. + deriv_p(:,2)**2.)
	oriented_m(:,1) = tens_m*deriv_m(:,1)/sqrt(deriv_m(:,1)**2. + deriv_m(:,2)**2.)
	oriented_m(:,2) = tens_m*deriv_m(:,2)/sqrt(deriv_m(:,1)**2. + deriv_m(:,2)**2.)


	!Finally, we take the derivative of oriented tension WITH RESPECT TO
	!REFERENCE ARC LENGTH to calculate the associated force

	boundary%stretch_force(:,1) = (oriented_p(:,1) - oriented_m(:,1))/ref_bnd%ds0
	boundary%stretch_force(:,2) = (oriented_p(:,2) - oriented_m(:,2))/ref_bnd%ds0
        
end subroutine calc_stretch_force


!=============================================================================================
!THIS SUBROUTINE INTEGRATES A VECTOR QUANTITY ON THE IMMERSED BOUNDARY. THE DIFFERENTIAL FORM
!MUST BE SPECIFIED SO THAT IT KNOWS WHAT ARC LENGTH TO INTEGRATE WITH RESPECT TO
subroutine vec_bnd_integrate(out,in,form)
!THE 'IN'-PUT IS A LENGTH M VECTOR QUANTITY DEFINED ON SOME IMMERSED BOUNDARY
REAL, DIMENSION(1:Mb,1:2) :: in
!THE 'OUT'-PUT IS A LENGTH 1 VECTOR QUANTITY THAT WILL BE CALCULATED
REAL, DIMENSION(1,1:2) :: out
!THE DIFFERENTIAL FORM IS A LENGTH M SCALAR QUANTITY DEFINED ON THE IMMERSED BOUNDARY
REAL, DIMENSION(1:Mb,1) :: form

!we simply multiply the input by the form & sum them
out(1,1) = sum(in(:,1)*form(:,1))
out(1,2) = sum(in(:,2)*form(:,1))
!AND WE'RE DONE

end subroutine vec_bnd_integrate


!=============================================================================================
!THIS SUBROUTINE INTEGRATES A SCALAR QUANTITY ON THE IMMERSED BOUNDARY. THE DIFFERENTIAL FORM
!MUST BE SPECIFIED SO THAT IT KNOWS WHAT ARC LENGTH TO INTEGRATE WITH RESPECT TO
subroutine scalar_bnd_integrate(out,in,form)
!THE 'IN'-PUT IS A LENGTH M VECTOR QUANTITY DEFINED ON SOME IMMERSED BOUNDARY
REAL, DIMENSION(1:Mb,1) :: in
!THE 'OUT'-PUT IS A LENGTH 1 VECTOR QUANTITY THAT WILL BE CALCULATED
REAL :: out
!THE DIFFERENTIAL FORM IS A LENGTH M SCALAR QUANTITY DEFINED ON THE IMMERSED BOUNDARY
REAL, DIMENSION(1:Mb,1) :: form

!we simply multiply the input by the form & sum them
out = sum(in*form)
!AND WE'RE DONE

end subroutine scalar_bnd_integrate

end module im_bnd_intern_mod
