! ============================================================================
! Name        : spectral_cortex.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : This is the main program which runs simulations.
! ============================================================================

program spectral_cortex

    USE fluid_sim_mod
    USE global_param_mod

implicit none
!Simple variables to track the number of plots and solves which have taken place
!Also a real to store the current time
INTEGER :: plot_count = 0, solvecount
REAL :: cur_time = 0.0


!Initialize the fluid simulation
call initialize_fluid
!Write out initial conditions to file.
call write_everything_out

!Calculate the total number of steps to take
steps = Tmax/dt

!Run the simulation for the number of steps calcualted
do solvecount = 0,steps
    !Calculate the current time
    cur_time = solvecount*dt
    !Perform an immersed boundary step
    call im_bnd_network_sim_step(cur_time)

    !If this step is a multiple of 'out_steps' write the current
    !state of the simulation to file.
    if (mod(solvecount+1,out_steps).eq.0) then
	    call write_everything_out
    end if

end do
!And we're done!

end program spectral_cortex
