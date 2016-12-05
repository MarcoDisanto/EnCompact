MODULE diffusive_term

USE variables,  !ONLY : n_flow_variables, uvwp, block
USE MPI_module, ONLY : ndims, myid
USE compact

IMPLICIT NONE

INTEGER, PARAMETER :: nn = n_flow_variables-1 !temporary upper bound

! For the time being the block type has been recycled to store diffusive terms as well
TYPE(block), DIMENSION(nn) :: diffvel ! single index, one for each component

CONTAINS

SUBROUTINE set_diff_bounds

  USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

  IMPLICIT NONE

  INTEGER            :: ic
  REAL               :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
  NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

  for_each_component: DO ic = 1, ndims

      diffvel(ic)%b = uvwp(ic)%b        ! only internal nodes
      ALLOCATE(diffvel(ic)%values(uvwp(ic)%b(1, 1) : uvwp(ic)%b(1, 2), &
                                  uvwp(ic)%b(2, 1) : uvwp(ic)%b(2, 2), &
                                  uvwp(ic)%b(3, 1) : uvwp(ic)%b(3, 2)))
      diffvel(ic)%values = NaN

  END DO for_each_component

END SUBROUTINE set_diff_bounds



SUBROUTINE diff_calc

  USE SPIKE
  USE essentials
  USE bandedmatrix
  USE Thomas_suite
  USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

  IMPLICIT NONE

  INTEGER :: ic, id, ider = 2, istag
  INTEGER :: j, k
  REAL, DIMENSION(:),       ALLOCATABLE        :: q
  REAL, DIMENSION(:, :, :), ALLOCATABLE        :: temp   ! used to temporarily store block-diagonal solution
  REAL                                         :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
  NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

  for_each_component: DO ic = 1, nn

    ! Debugging allocation
    diffvel(ic)%values = NaN

    ! I need this block to apply high level SPIKE logic. If I called SPIKE_solve2 at id = 2 without temp, I would have
    ! solved a reduced system using d/dx + d/dy as coupling terms, which is not what I want.
    ALLOCATE(temp(diffvel(ic)%b(1, 1):diffvel(ic)%b(1, 2), &
                  diffvel(ic)%b(2, 1):diffvel(ic)%b(2, 2), &
                  diffvel(ic)%b(3, 1):diffvel(ic)%b(3, 2)))

    !!!!!!!!!! x derivatives !!!!!!!!!!
    id = 1
    istag = logical2integer(id==ic)+1
    ALLOCATE(q(uvwp(ic)%b_bo(id, 1):uvwp(ic)%b_bo(id, 2)))
    DO j = diffvel(ic)%b(2, 1), diffvel(ic)%b(2, 2)
      DO k = diffvel(ic)%b(3, 1), diffvel(ic)%b(3, 2)

        q = uvwp(ic)%values(:, j, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    diffvel(ic)%values(:, j, k))

      END DO
    END DO
    DEALLOCATE(q)
    CALL SPIKE_solve2(id, ider, istag, diffvel(ic)%values)


    !!!!!!!!!! y derivatives !!!!!!!!!!
    id = 2
    istag = logical2integer(id==ic)+1
    ALLOCATE(q(uvwp(ic)%b_bo(id, 1):uvwp(ic)%b_bo(id, 2)))
    DO j = diffvel(ic)%b(1, 1), diffvel(ic)%b(1, 2)
      DO k = diffvel(ic)%b(3, 1), diffvel(ic)%b(3, 2)

        q = uvwp(ic)%values(j, :, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, :, k))

      END DO
    END DO
    DEALLOCATE(q)
    CALL SPIKE_solve2(id, ider, istag, temp)
    diffvel(ic)%values = diffvel(ic)%values + temp


    !!!!!!!!!! z derivatives !!!!!!!!!!
    id = 3
    istag = logical2integer(id==ic)+1
    ALLOCATE(q(uvwp(ic)%b_bo(id, 1):uvwp(ic)%b_bo(id, 2)))
    DO j = diffvel(ic)%b(1, 1), diffvel(ic)%b(1, 2)
      DO k = diffvel(ic)%b(2, 1), diffvel(ic)%b(2, 2)

        q = uvwp(ic)%values(j, k, :)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, k, :))

      END DO
    END DO
    DEALLOCATE(q)
    CALL SPIKE_solve2(id, ider, istag, temp)
    diffvel(ic)%values = diffvel(ic)%values + temp


    DEALLOCATE(temp)
  END DO for_each_component


END SUBROUTINE diff_calc

END MODULE diffusive_term
