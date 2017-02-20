MODULE time_advancement

USE variables
USE essentials

IMPLICIT NONE

TYPE(block), DIMENSION(:), ALLOCATABLE :: gradp    ! block storing pressure gradients
TYPE(block), DIMENSION(:), ALLOCATABLE :: totder   ! support block used to store sum of derivatives for RK method
! residuals (velocity differences)
REAL, DIMENSION(3)                     :: vel_diff_proc       ! max velocity components differences for each process
REAL                                   :: max_vel_diff_proc   ! max velocity difference among all components, but for each process
REAL                                   :: vel_res             ! max difference among all components and processes (true residual)


CONTAINS


  SUBROUTINE set_grad_p

    !USE variables
    USE MPI_module, ONLY: ndims
    USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

    IMPLICIT NONE

    INTEGER            :: i
    REAL               :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
    NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

    ALLOCATE(gradp(ndims))

    DO i = 1, ndims
      gradp(i)%b = uvwp(i)%b
      ALLOCATE(gradp(i)%values(gradp(i)%b(1, 1):gradp(i)%b(1, 2), &
                               gradp(i)%b(2, 1):gradp(i)%b(2, 2), &
                               gradp(i)%b(3, 1):gradp(i)%b(3, 2)))
      gradp(i)%values = NaN
    END DO

  END SUBROUTINE set_grad_p


  SUBROUTINE set_tot_der

    !USE variables
    USE MPI_module, ONLY: ndims
    USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

    IMPLICIT NONE

    INTEGER            :: i
    REAL               :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
    NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

    ALLOCATE(totder(ndims))

    DO i = 1, ndims
      totder(i)%b = uvwp(i)%b
      ALLOCATE(totder(i)%values(totder(i)%b(1, 1):totder(i)%b(1, 2), &
                               totder(i)%b(2, 1):totder(i)%b(2, 2), &
                               totder(i)%b(3, 1):totder(i)%b(3, 2)))
      totder(i)%values = NaN
    END DO

  END SUBROUTINE set_tot_der



  SUBROUTINE divergence_calc
    ! Routine calculating divergence of starred velocity field, that is to say the
    ! source term of the PPE. This quantity is stored inside the variable "b", whose
    ! name was chosen for the sake of clarity...cough! cough!

    USE set_pressure
    USE bandedmatrix

    IMPLICIT NONE

    INTEGER, PARAMETER              :: nn = n_flow_variables-1
    INTEGER                         :: i, j, k
    REAL, DIMENSION(:), ALLOCATABLE :: vct          ! just to meet the "ALLOCATABLE" requirement for CDS multiplication

    ! Initialization
    b = 0

    ! du/dx + ...
    i = 1
    ALLOCATE(vct(GraDiv(i, 2)%lb(2):GraDiv(i, 2)%ub(2)))
    DO j = uvwp(i)%b(2, 1), uvwp(i)%b(2, 2)
      DO k = uvwp(i)%b(3, 1), uvwp(i)%b(3, 2)

        vct = uvwp(i)%values(GraDiv(i, 2)%lb(2):GraDiv(i, 2)%ub(2), j, k)
        b(:, j, k) = b(:, j, k) + GraDiv(i, 2)*vct

      END DO
    END DO
    DEALLOCATE(vct)

    ! ... + dv/dy + ...
    i = 2
    ALLOCATE(vct(GraDiv(i, 2)%lb(2):GraDiv(i, 2)%ub(2)))
    DO j = uvwp(i)%b(1, 1), uvwp(i)%b(1, 2)
      DO k = uvwp(i)%b(3, 1), uvwp(i)%b(3, 2)

        vct = uvwp(i)%values(j, GraDiv(i, 2)%lb(2):GraDiv(i, 2)%ub(2), k)
        b(j, :, k) = b(j, :, k) + GraDiv(i, 2)*vct

      END DO
    END DO
    DEALLOCATE(vct)

    ! ... + dw/dz
    i = 3
    ALLOCATE(vct(GraDiv(i, 2)%lb(2):GraDiv(i, 2)%ub(2)))
    DO j = uvwp(i)%b(1, 1), uvwp(i)%b(1, 2)
      DO k = uvwp(i)%b(2, 1), uvwp(i)%b(2, 2)

        vct = uvwp(i)%values(j, k, GraDiv(i, 2)%lb(2):GraDiv(i, 2)%ub(2))
        b(j, k, :) = b(j, k, :) + GraDiv(i, 2)*vct

      END DO
    END DO
    DEALLOCATE(vct)

  END SUBROUTINE divergence_calc



  SUBROUTINE p_grad_calc
    ! Routine calculating pressure gradients used in the correction step.

    USE set_pressure
    USE bandedmatrix

    IMPLICIT NONE

    INTEGER, PARAMETER              :: nn = n_flow_variables-1
    INTEGER                         :: i, j, k
    REAL, DIMENSION(:), ALLOCATABLE :: vct          ! just to meet the "ALLOCATABLE" requirement for CDS multiplication
    REAL, DIMENSION(:), ALLOCATABLE :: vct2         ! this is due to the fact that GraDiv calculates gradients even at the boundaries
                                                    ! and they are not meaningful


    !!!!!!!!!! dp/dx !!!!!!!!!!
    i = 1
    ALLOCATE(vct(GraDiv(i, 1)%lb(2):GraDiv(i, 1)%ub(2)))
    ALLOCATE(vct2(GraDiv(i, 1)%lb(1):GraDiv(i, 1)%ub(1)))
    DO j = gradp(i)%b(2, 1), gradp(i)%b(2, 2)
      DO k = gradp(i)%b(3, 1), gradp(i)%b(3, 2)

        vct = uvwp(4)%values(GraDiv(i, 1)%lb(2):GraDiv(i, 1)%ub(2), j, k)
        vct2 = GraDiv(i, 1)*vct
        gradp(i)%values(:, j, k) = vct2(gradp(i)%b(i, 1):gradp(i)%b(i, 2))

      END DO
    END DO
    DEALLOCATE(vct)
    DEALLOCATE(vct2)

    !!!!!!!!!! dp/dy !!!!!!!!!!
    i = 2
    ALLOCATE(vct(GraDiv(i, 1)%lb(2):GraDiv(i, 1)%ub(2)))
    ALLOCATE(vct2(GraDiv(i, 1)%lb(1):GraDiv(i, 1)%ub(1)))
    DO j = gradp(i)%b(1, 1), gradp(i)%b(1, 2)
      DO k = gradp(i)%b(3, 1), gradp(i)%b(3, 2)

        vct = uvwp(4)%values(j, GraDiv(i, 1)%lb(2):GraDiv(i, 1)%ub(2), k)
        vct2 = GraDiv(i, 1)*vct
        gradp(i)%values(j, :, k) = vct2(gradp(i)%b(i, 1):gradp(i)%b(i, 2))

      END DO
    END DO
    DEALLOCATE(vct)
    DEALLOCATE(vct2)

    !!!!!!!!!! dp/dy !!!!!!!!!!
    i = 3
    ALLOCATE(vct(GraDiv(i, 1)%lb(2):GraDiv(i, 1)%ub(2)))
    ALLOCATE(vct2(GraDiv(i, 1)%lb(1):GraDiv(i, 1)%ub(1)))
    DO j = gradp(i)%b(1, 1), gradp(i)%b(1, 2)
      DO k = gradp(i)%b(2, 1), gradp(i)%b(2, 2)

        vct = uvwp(4)%values(j, k, GraDiv(i, 1)%lb(2):GraDiv(i, 1)%ub(2))
        vct2 = GraDiv(i, 1)*vct
        gradp(i)%values(j, k, :) = vct2(gradp(i)%b(i, 1):gradp(i)%b(i, 2))

      END DO
    END DO
    DEALLOCATE(vct)
    DEALLOCATE(vct2)


  END SUBROUTINE p_grad_calc



  SUBROUTINE residual_eval(dt)

    USE MPI

    INTEGER          :: i, ierr
    REAL, INTENT(IN) :: dt

    DO i = 1, 3

      vel_diff_proc(i) = MAXVAL(ABS(uvwp(i)%values(uvwp(i)%b(1, 1) : uvwp(i)%b(1, 2), &
                                                   uvwp(i)%b(2, 1) : uvwp(i)%b(2, 2), &
                                                   uvwp(i)%b(3, 1) : uvwp(i)%b(3, 2)) - uvwp_old(i)%values))

    END DO

    max_vel_diff_proc = MAXVAL(vel_diff_proc)
    CALL MPI_ALLREDUCE(max_vel_diff_proc, vel_res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    vel_res = vel_res/dt


  END SUBROUTINE residual_eval



  SUBROUTINE ExplEuler(dt, ni)
    ! Routine that calculates the velocity field at the next time step. This version
    ! supports explicit Euler formula only. More sophisticated time schemes will
    ! be implemented in the near future.

    USE diffusive_term
    USE convective_term
    USE solve_pressure
    USE SPIKE

    IMPLICIT NONE

    REAL, INTENT(IN) :: dt, ni
    INTEGER          :: ic

    ! saving old values of flow variables
    DO ic = 1, n_flow_variables-1
      uvwp_old(ic)%values = uvwp(ic)%values(uvwp(ic)%b(1, 1) : uvwp(ic)%b(1, 2), &
                                            uvwp(ic)%b(2, 1) : uvwp(ic)%b(2, 2), &
                                            uvwp(ic)%b(3, 1) : uvwp(ic)%b(3, 2))
    END DO

    CALL diff_calc
    CALL conv_calc

    !!!!!!!!!! Starred field !!!!!!!!!!
    DO ic = 1, ndims
      ! only internal faces must be updated
      uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                      uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                      uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) = uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                                                                           uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                                                                           uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) + &
                                                           dt*ni*diffvel(ic)%values                           + &
                                                           dt*( -convvel(ic)%values)
    END DO
    CALL SPIKE_exchange_uvw

    !!!!!!!!!! Divergence !!!!!!!!!!
    CALL divergence_calc

    !!!!!!!!!! Pressure term !!!!!!!!!!
    CALL Solve_p
    CALL p_grad_calc

    !!!!!!!!!! Pressure correction !!!!!!!!!!
    DO ic = 1, ndims
      ! only internal faces must be updated
      uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                      uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                      uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) = uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                                                                           uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                                                                           uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) + &
                                                           ( -gradp(ic)%values)
    END DO
    CALL SPIKE_exchange_uvw



  END SUBROUTINE ExplEuler



  SUBROUTINE RK4(dt, ni)
    ! Routine that calculates the velocity field at the next time step by means of
    ! fourth order Runge-Kutta method

    USE diffusive_term
    USE convective_term
    USE solve_pressure
    USE SPIKE

    IMPLICIT NONE

    REAL, INTENT(IN)   :: dt, ni
    INTEGER            :: ic, istg
    REAL, DIMENSION(4) :: aa, bb

    aa = [0.5, 0.5, 1., 1.]
    bb = [1./6, 1./3, 1./3, 1./6]

    ! saving old values of flow variables
    DO ic = 1, n_flow_variables-1
      uvwp_old(ic)%values = uvwp(ic)%values(uvwp(ic)%b(1, 1) : uvwp(ic)%b(1, 2), &
                                            uvwp(ic)%b(2, 1) : uvwp(ic)%b(2, 2), &
                                            uvwp(ic)%b(3, 1) : uvwp(ic)%b(3, 2))
      totder(ic)%values = 0
    END DO

    !!!!!!!!!! Stages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO istg = 1, 4
      CALL diff_calc
      CALL conv_calc

      DO ic = 1, ndims
        ! only internal faces must be updated
        uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                        uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                        uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) = uvwp_old(ic)%values        + &
                                                             aa(istg)*dt*ni*diffvel(ic)%values + &
                                                             aa(istg)*dt*( -convvel(ic)%values)
        ! accumulating derivative values
        totder(ic)%values = totder(ic)%values + bb(istg)*dt*(-convvel(ic)%values + ni*diffvel(ic)%values)
      END DO
      CALL SPIKE_exchange_uvw

      !!!!!!!!!! Divergence !!!!!!!!!!
      CALL divergence_calc

      !!!!!!!!!! Pressure term !!!!!!!!!!
      CALL Solve_p
      CALL p_grad_calc

      !!!!!!!!!! Pressure correction !!!!!!!!!!
      DO ic = 1, ndims
        ! only internal faces must be updated
        uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                        uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                        uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) = uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                                                                             uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                                                                             uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) + &
                                                             ( -gradp(ic)%values)
        ! accumulating derivative values
        totder(ic)%values = totder(ic)%values + bb(istg)*(-gradp(ic)%values)
        CALL SPIKE_exchange_uvw



      END DO
    END DO


    !!!!!!!!!! Final evaluation !!!!!!!!!!
    DO ic = 1, ndims
      ! only internal faces must be updated
      uvwp(ic)%values(uvwp(ic)%b(1, 1):uvwp(ic)%b(1, 2), &
                      uvwp(ic)%b(2, 1):uvwp(ic)%b(2, 2), &
                      uvwp(ic)%b(3, 1):uvwp(ic)%b(3, 2)) = uvwp_old(ic)%values + &
                                                           totder(ic)%values
    END DO
    CALL SPIKE_exchange_uvw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE RK4



END MODULE time_advancement
