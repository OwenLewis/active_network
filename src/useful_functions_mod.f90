! ============================================================================
! Name        : useful_functions_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing an assortment of useful functions
! ============================================================================

module useful_functions_mod
USE global_param_mod
implicit none
PRIVATE
PUBLIC :: mymod, better_phi, four_part_wave
contains

!========================================================================================================

	!Function to apply 'mod' in a way that returns values 1 through 'base' and
	!treats negative integers in the manner that works best for my spreading and
	!interpolating operators
	function mymod(argint,base)
		!'argint' is the integer that will be mod'ed
		! base defines the base of the modular arithmetic
		INTEGER :: argint, base
		!RETURNED VALUE WILL BE BETWEEN 1 AND 'BASE'
		INTEGER :: mymod
		!initialize the returned variable
		mymod = argint
		!IF THE ARGUMENT IS 'NEGATIVE' ADD MULTIPLES OF
		!THE 'BASE' UNTIL THIS IS NO LONGER TRUE
		do
			if (mymod > 1) exit
			mymod = mymod + base
		end do
		!subract one, mod by base, then add one
		!THIS IS WHAT GETS RETURNED VALUES IN 1 TO BASE
		!INSTEAD OF 0 TO BASE-1 (WHICH IS DEFAULT)
		mymod = mymod - 1
		mymod = mod(mymod,base)
		mymod = mymod + 1
	end function mymod


!========================================================================================================

	!Function to apply Peskin's simple phi function to a small array of distances
	function simple_phi(distance)
		!Return array is of size 4 (reflecting the footprint of the phi function)
		REAL, DIMENSION(1:4) :: simple_phi
		!Take in an array of size 4 for same reason
		REAL, INTENT(IN), DIMENSION(1:4) :: distance
		!Here I use Peskin's original phi function (the cosine)
		simple_phi = (1.0 + cos(PI*distance/2.0) )/4.0

	end function simple_phi

!========================================================================================================

	!Function to apply Peskin's well derived approximate phi function
	function better_phi(distance)
		!Return array is of size 4 (reflecting the footprint of the phi function)
		REAL, DIMENSION(1:4) :: better_phi
		!Take in an array of size 4 for the same reasons
		REAL, INTENT(IN), DIMENSION(1:4) :: distance
		!Here I use Peskin's original delta function (the cosine)
		better_phi(1) = (5. + 2.*distance(1) - (-7. - 12.*distance(1) - 4.*distance(1)**2.)**(1./2.) )/8.
		better_phi(2) = (3. + 2.*distance(2) + (1. - 4.*distance(2) - 4.*distance(2)**2.)**(1./2.) )/8.
		better_phi(3) = (3. - 2.*distance(3) + (1. + 4.*distance(3) - 4.*distance(3)**2.)**(1./2.) )/8.
		better_phi(4) = (5. - 2.*distance(4) - (-7. + 12.*distance(4) - 4.*distance(4)**2.)**(1./2.) )/8.

	end function better_phi

!========================================================================================================


    !This is a function to apply my piecewise defined 4 part wave. It should take in body coordinate array
    function four_part_wave(body_coord,time)
        !Input array is real, but of unknown size
        REAL, INTENT(IN), DIMENSION(:) :: body_coord
        REAL, INTENT(IN) :: time
        !Output array is a real of unknown size (it will be the same as input)
        REAL, DIMENSION(lbound(body_coord,1):ubound(body_coord,1)) :: four_part_wave, shifted_coord

        shifted_coord = con_wave_num*body_coord - freq*time
        shifted_coord = modulo(shifted_coord,2*PI)

        where(shifted_coord.lt.(PI/2)) four_part_wave = 0.
        where((shifted_coord.gt.(PI/2)).and.(shifted_coord.lt.PI)) four_part_wave = 2*shifted_coord/PI - 1.
        where((shifted_coord.gt.PI).and.(shifted_coord.lt.(3*PI/2))) four_part_wave = 1.
        where(shifted_coord.gt.(3.*PI/2.)) four_part_wave = 4-shifted_coord*(2./Pi)


    end function four_part_wave


end module useful_functions_mod
