! ============================================================================
! Name        : stokes_2D_fft_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module to perform a single Stokes Solve using a 2D fft poisson
!				solver
! ============================================================================
module stokes_2D_fft_mod
USE global_param_mod
implicit none
INCLUDE 'fftw3.f'
PRIVATE
PUBLIC :: stokes_solve_2D_fft
!EIGENVALUES FOR LAPLACE AND DIFFERENTIATION
REAL, DIMENSION(:,:), ALLOCATABLE:: laplace_eigs
COMPLEX, DIMENSION(:,:), ALLOCATABLE :: diff_eigs_x, diff_eigs_y
!DUMBY VARIABLES TO PERFORM POISSON SOLVES
COMPLEX, DIMENSION(:,:), ALLOCATABLE :: data_real, data_hat
!DUMBY INTEGERS FOR ITERATIONS
INTEGER :: I,K
!LOGICAL VALUE FOR INITIALIZATION FLAG
INTEGER :: isInit = 0
!INTEGERS TO HOLD FFTW PLANS
INTEGER*8 :: fft, ifft


!========================================================================================================

contains
!Subroutine to perform a single stokes solve
!forces are input, velocity and pressure are output
	subroutine stokes_solve_2d_fft(input_force_horiz,input_force_vert,output_vel_horiz,output_vel_vert,output_pres)
    !All arrays of input and output are Nx by Ny and real
	REAL, DIMENSION(1:Nx,1:Ny) :: input_force_horiz, input_force_vert, output_vel_horiz, output_vel_vert, output_pres
	!All arrays used to do computations are of size Nx by Ny and complex
	COMPLEX, DIMENSION(1:Nx,1:Ny) :: output_pres_hat, output_vel_h_hat, output_vel_v_hat
	COMPLEX, DIMENSION(1:Nx,1:Ny) :: force_horiz_hat, grad_pres_horiz, force_vert_hat, grad_pres_vert, div_force_hat
	COMPLEX, DIMENSION(1:Nx,1:Ny) :: rhs_horiz, rhs_vert


	!Check to see if initialization has happened. If not, do it & flag
	if (isInit == 0) then
	    print*, 'We are initializing all of the fft plans and differentiation operators'
		call stokes_initialize
		isInit = 1
	end if

    !We take the FFT of the input forces, Horizontal first
        data_real = input_force_horiz!copy real data in complex array
        call dfftw_execute_dft(fft,data_real,data_hat) !take the fourier transform
        force_horiz_hat = data_hat !store that output


    !Now vertical
        data_real = input_force_vert!copy real data into complex array
        call dfftw_execute_dft(fft,data_real,data_hat)!take the fourier transform
        force_vert_hat = data_hat!store that output


	!Take the divergence of the input force (in frequency space)
	call divergence(force_horiz_hat,force_vert_hat,div_force_hat)

	!solve a cell centered poisson problem
	!This gives us the pressure, BUT IN FREQUENCY SPACE STILL
	call poisson_solve_2D_fft(div_force_hat,output_pres_hat)

	!Take the gradient of the pressure
	call gradient(output_pres_hat,grad_pres_horiz,grad_pres_vert)

	!form the right hand side of the poisson solve for velocity
	rhs_horiz = (grad_pres_horiz - force_horiz_hat)/mu
	rhs_vert = (grad_pres_vert - force_vert_hat)/mu

	! call a poisson solve on the horiz and vert components
	call poisson_solve_2D_fft(rhs_horiz,output_vel_h_hat)
	call poisson_solve_2D_fft(rhs_vert,output_vel_v_hat)

	!Now we need to get everything back to real space
	!dont forget to normalize by the grid size

	!Horizontal Velocity
	data_hat = output_vel_h_hat
	call dfftw_execute_dft(ifft,data_hat,data_real)
	output_vel_horiz = data_real/(Nx*Ny)

    !Vertical Velocity
    data_hat = output_vel_v_hat
    call dfftw_execute_dft(ifft,data_hat,data_real)
    output_vel_vert = data_real/(Nx*Ny)

	!Pressure
    data_hat = output_pres_hat
    call dfftw_execute_dft(ifft,data_hat,data_real)
    output_pres = data_real/(Nx*Ny)

    !Now we will normalize the 'baseline pressure' so that it is zero 'outside' the cell
    output_pres = output_pres - output_pres(1,1)


01	end subroutine stokes_solve_2d_fft
!========================================================================================================

!Initialization sub-routine to populate arrays full of eigenvalues of the discrete
!2D laplacian for the 3 fluid variables (sizes vary on a MAC grid)
!Also generates the plans for the forward & inverse transforms involved in poisson solve
	subroutine stokes_initialize
	INTEGER :: up_rangex, down_rangex, up_rangey, down_rangey


        allocate(laplace_eigs(1:Nx,1:Ny),diff_eigs_x(1:Nx,1:Ny),diff_eigs_y(1:Nx,1:Ny))
        allocate(data_real(1:Nx,1:Ny),data_hat(1:Nx,1:Ny))

		!To populate the differentiation operators we will need to do some clever indexing.

		!We check to see if we are on an even or odd sized grid
		if (mod(Nx,2).eq.0) then !If even size
		    !We don't quite get Nx/2 positive frequencies
		    up_rangex = Nx/2 - 1
		    !But we do get Nx/2 negative frequencies
		    down_rangex = -Nx/2
		else !If the grid is of odd size
		    !We get Nx/2 positive frequencies
		    up_rangex = floor(Nx/2.)
		    !And Nx/2 negative frequencies
		    down_rangex = -up_rangex
        endif


        !We check to see if we are on an even or odd sized grid
		if (mod(Ny,2).eq.0) then !If even size
            !We don't quite get Nx/2 positive frequencies
		    up_rangey = Ny/2 - 1
            !But we do get Nx/2 negative frequencies
		    down_rangey = -Ny/2
        else !If the grid is of odd size
            !We get Nx/2 positive frequencies
            up_rangey = floor(Ny/2.)
            !And Nx/2 negative frequencies
            down_rangey = -up_rangey
        endif


		!Now I will generate eigenvalues of the first derivative in both directions
        !Notice that each index is bumped up by one. This is to get the 'zeroth' frequency
        !in the first entry of the array
        do I = down_rangex,up_rangex
            do K = down_rangey,up_rangey
                !This is the derivative in the x direction
                diff_eigs_x(modulo(I,Nx)+1,modulo(K,Ny)+1) = (2.0*(I)*PI*J/Xmax)
                !This is the derivative in the y direction
                diff_eigs_y(modulo(I,Nx)+1,modulo(K,Ny)+1) = (2.0*(K)*PI*J/Ymax)
            end do
        end do

        !Generate eigenvalues of laplacian for periodic scalar field
        !In frequency space, shit sure is easy
        laplace_eigs = diff_eigs_x**2 + diff_eigs_y**2
        !Be sure to zero out the first eigenvalue
        laplace_eigs(1,1) = 0


		!MAKE FORWARD & BACKWARD PLANS FOR PERIODIC SCALARS
		call dfftw_plan_dft_2d(fft,Nx,Ny,data_real,data_hat,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_plan_dft_2d(ifft,Nx,Ny,data_hat,data_real,FFTW_BACKWARD,FFTW_ESTIMATE)


02	end subroutine stokes_initialize

!========================================================================================================

!Subroutine to take the gradient of the pressure (or any scalar field)
!NOTICE THAT THIS IS ALL DONE IN FREQUENCY SPACE. I AM ASSUMING THAT FOURIER
!TRANSFORMS HAVE ALREADY BEEN DONE IN THE CALLING FUNCTION. THIS IS FOR SPEED.
	subroutine gradient(input_scalar_hat,output_horiz_hat,output_vert_hat)

	!Remember that on a colocated grid, vector components are of the same
	!rank as cell-centered scalar fields
	COMPLEX, DIMENSION(1:Nx,1:Ny) :: input_scalar_hat, output_horiz_hat, output_vert_hat


	!I take the derivative in the x direction by multiplying by the corresponding
	!eigenvalues stored already
	output_horiz_hat = diff_eigs_x*input_scalar_hat

	!I take the derivative in the x direction by multiplying by the corresponding
    !eigenvalues stored already
    output_vert_hat = diff_eigs_y*input_scalar_hat

03	end subroutine gradient

!========================================================================================================

!Subroutine to take the divergence of the velocity (or any vector field)
!NOTICE THAT THIS IS ALL DONE IN FREQUENCY SPACE. I AM ASSUMING THAT FOURIER
!TRANSFORMS HAVE ALREADY BEEN DONE IN THE CALLING FUNCTION. THIS IS FOR SPEED.
	subroutine divergence(input_horiz_hat,input_vert_hat,output_scalar_hat)

	COMPLEX, DIMENSION(1:Nx,1:Ny) :: input_horiz_hat,input_vert_hat,output_scalar_hat

	!We take the divergence of the vector field by multiplying by eigenvalues of x and y
	!derivatives, then summing the result
	output_scalar_hat = input_horiz_hat*diff_eigs_x + input_vert_hat*diff_eigs_y

04	end subroutine divergence

!========================================================================================================

!Subroutine to solve scalar valued poisson problem on a periodic grid
!THIS IS USED TO SOLVE FOR PRESSURE AND BOTH COMPONENTS OF VELOCITY
	subroutine poisson_solve_2d_fft(input_hat,output_hat)
		COMPLEX, DIMENSION(1:Nx,1:Ny) :: input_hat, output_hat


		!Divide transformed data by eigenvalues of laplacian
		output_hat = input_hat/laplace_eigs
		!make sure to zero out the component in null-space
		output_hat(1,1) = 0


	end subroutine poisson_solve_2d_fft


end module stokes_2D_fft_mod
