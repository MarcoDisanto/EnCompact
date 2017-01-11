MODULE convective_term

USE variables,  !ONLY : n_flow_variables, uvwp, block
USE MPI_module, ONLY : ndims, mycoords, dims, myid, ndims, idm, idp, procs_grid
USE compact,    ONLY : cmp

IMPLICIT NONE

INTEGER, PARAMETER :: nn = n_flow_variables-1 !temporary upper bound

! For the time being the block type has been recycled to store convective terms as well
TYPE(block), DIMENSION(2, nn)  :: velvel    ! velvel(i, j) = uvwp(i)*uvwp(j) not the exact formal operation, but you get the idea;
                                            ! on the first row only the homogeneus products are stored (i.e. u^2, v^2 and w^2), while
                                            ! the cross-products (i.e. uv, vw and wu) are put on the second row.
TYPE(block), DIMENSION(nn)     :: intvel    ! Interpolated values of velocity; only used for cross-products uv, vw,and wu;
                                            ! values are calculated and immediately stored (or multiplied) in velvel so as to free this variable
TYPE(block), DIMENSION(nn)     :: convvel   ! Convective term for momentum equation

INTEGER, DIMENSION(:, :), ALLOCATABLE :: MPI_slab_m, MPI_slab_p          ! Almost identical to MPI_flats

CONTAINS

  SUBROUTINE set_conv_bounds

    USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines
    USE essentials
    USE MPI,        ONLY : MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_ORDER_FORTRAN, MPI_PROC_NULL, MPI_ADDRESS_KIND

    IMPLICIT NONE

    INTEGER                               :: lid, uid      ! number of lower and upper diagonals
    INTEGER                               :: i, ic, id, id2, j
    INTEGER                               :: ierr
    INTEGER, DIMENSION(ndims)             :: array_of_sizes
    INTEGER, DIMENSION(ndims)             :: array_of_subsizes_m, array_of_subsizes_p
    INTEGER, DIMENSION(ndims)             :: array_of_starts
    REAL               :: r = 0.0, NaN  ! r is a dummy real used to define a NaN of type real
    NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r


    DO i = 1, ndims

      !!!!!!!!!! Convective term !!!!!!!!!!
      ! This block contains data regarding internal nodes only (the ones being updated)
      convvel(i)%b = uvwp(i)%b
      ! Allocation
      ALLOCATE(convvel(i)%values(convvel(i)%b(1, 1):convvel(i)%b(1, 2), &
                                 convvel(i)%b(2, 1):convvel(i)%b(2, 2), &
                                 convvel(i)%b(3, 1):convvel(i)%b(3, 2)))
      convvel(i)%values = NaN

      !!!!!!!!!! velocity*velocity !!!!!!!!!!
      ! Homogeneus products
      ! NOTE: number of diagonals of B(i, 1, 1); overlap is determined by how many
      ! extra values are needed for derivation of u^2, v^2 and w^2 with staggering 1 (cell to face)
      lid = cmp(i, 1, 1)%B%lid
      uid = cmp(i, 1, 1)%B%uid
      ! Internal values
      velvel(1, i)%b = uvwp(i)%b
      velvel(1, i)%b(i, 1) = velvel(1, i)%b(i, 1) - logical2integer(mycoords(i)==0)
      ! Internal values + overlap (fictional overlaps included)
      velvel(1, i)%b_ol = velvel(1, i)%b
      velvel(1, i)%b_ol(i, 1) = velvel(1, i)%b_ol(i, 1) - lid
      velvel(1, i)%b_ol(i, 2) = velvel(1, i)%b_ol(i, 2) + uid
      ! Internal values + boundary
      velvel(1, i)%b_bc = velvel(1, i)%b
      velvel(1, i)%b_bc(i, 1) = velvel(1, i)%b_bc(i, 1) - logical2integer(mycoords(i)==0)
      velvel(1, i)%b_bc(i, 2) = velvel(1, i)%b_bc(i, 2) + logical2integer(mycoords(i)==dims(i)-1)
      ! Boundary and/or overlap
      velvel(1, i)%b_bo = velvel(1, i)%b
      velvel(1, i)%b_bo(i, 1) = velvel(1, i)%b_bo(i, 1) - lid*(logical2integer(mycoords(i)/=0)) &
                                                        - logical2integer(mycoords(i)==0)
      velvel(1, i)%b_bo(i, 2) = velvel(1, i)%b_bo(i, 2) + uid*(logical2integer(mycoords(i)/=dims(i)-1)) &
                                                        + logical2integer(mycoords(i)==dims(i)-1)
      ! Allocation
      ALLOCATE(velvel(1, i)%values(velvel(1, i)%b_bo(1, 1):velvel(1, i)%b_bo(1, 2), &
                                   velvel(1, i)%b_bo(2, 1):velvel(1, i)%b_bo(2, 2), &
                                   velvel(1, i)%b_bo(3, 1):velvel(1, i)%b_bo(3, 2)))
      velvel(1, i)%values = NaN

      ! Cross products
      j = i+1
      IF (j>ndims) j = 1
      ! NOTE: number of diagonals of B(i, 1, 2); overlap is determined by how many
      ! extra values are needed for derivation of uv, vw and wu with staggering 2 (face to cell)
      ! For the time being the same schemes are assumed to hold for both directions i and j
      lid = cmp(i, 1, 2)%B%lid
      uid = cmp(i, 1, 2)%B%uid
      ! Internal values
      velvel(2, i)%b = uvwp(i)%b
      velvel(2, i)%b(j, 1) = velvel(2, i)%b(j, 1) + logical2integer(mycoords(j)==0)
      velvel(2, i)%b(j, 2) = velvel(2, i)%b(j, 2)
      ! Internal values + overlap (fictional overlaps included)
      velvel(2, i)%b_ol = velvel(2, i)%b
      velvel(2, i)%b_ol(i, 1) = velvel(2, i)%b_ol(i, 1) - lid
      velvel(2, i)%b_ol(i, 2) = velvel(2, i)%b_ol(i, 2) + uid
      velvel(2, i)%b_ol(j, 1) = velvel(2, i)%b_ol(j, 1) - lid
      velvel(2, i)%b_ol(j, 2) = velvel(2, i)%b_ol(j, 2) + uid
      ! Internal values + boundary
      velvel(2, i)%b_bc = velvel(2, i)%b
      velvel(2, i)%b_bc(i, 1) = velvel(2, i)%b_bc(i, 1) - logical2integer(mycoords(i)==0)
      velvel(2, i)%b_bc(i, 2) = velvel(2, i)%b_bc(i, 2) + logical2integer(mycoords(i)==dims(i)-1)
      velvel(2, i)%b_bc(j, 1) = velvel(2, i)%b_bc(j, 1) - logical2integer(mycoords(j)==0)
      velvel(2, i)%b_bc(j, 2) = velvel(2, i)%b_bc(j, 2) + logical2integer(mycoords(j)==dims(j)-1)
      ! Boundary and/or overlap
      velvel(2, i)%b_bo = velvel(2, i)%b
      velvel(2, i)%b_bo(i, 1) = velvel(2, i)%b_bo(i, 1) - lid*(logical2integer(mycoords(i)/=0)) &
                                                        - logical2integer(mycoords(i)==0)
      velvel(2, i)%b_bo(i, 2) = velvel(2, i)%b_bo(i, 2) + uid*(logical2integer(mycoords(i)/=dims(i)-1)) &
                                                        + logical2integer(mycoords(i)==dims(i)-1)
      velvel(2, i)%b_bo(j, 1) = velvel(2, i)%b_bo(j, 1) - lid*(logical2integer(mycoords(j)/=0)) &
                                                        - logical2integer(mycoords(j)==0)
      velvel(2, i)%b_bo(j, 2) = velvel(2, i)%b_bo(j, 2) + uid*(logical2integer(mycoords(j)/=dims(j)-1)) &
                                                        + logical2integer(mycoords(j)==dims(j)-1)

      ! Allocation
      ALLOCATE(velvel(2, i)%values(velvel(2, i)%b_bo(1, 1):velvel(2, i)%b_bo(1, 2), &
                                   velvel(2, i)%b_bo(2, 1):velvel(2, i)%b_bo(2, 2), &
                                   velvel(2, i)%b_bo(3, 1):velvel(2, i)%b_bo(3, 2)))
      velvel(2, i)%values = NaN

      !!!!!!!!!! interpolated velocity !!!!!!!!!!
      ! NOTE: same dimensions as velvel(2, i)
      intvel(i)%b = velvel(2, i)%b
      intvel(i)%b_ol = velvel(2, i)%b_ol
      intvel(i)%b_bc = velvel(2, i)%b_bc
      intvel(i)%b_bo = velvel(2, i)%b_bo

      ! Allocation
      ALLOCATE(intvel(i)%values(intvel(i)%b_bo(1, 1):intvel(i)%b_bo(1, 2), &
                                intvel(i)%b_bo(2, 1):intvel(i)%b_bo(2, 2), &
                                intvel(i)%b_bo(3, 1):intvel(i)%b_bo(3, 2)))
      intvel(i)%values = NaN

    END DO


    !!!!!!!!!! Slabs definition !!!!!!!!!!
    ALLOCATE(MPI_slab_m(ndims, ndims))
    ALLOCATE(MPI_slab_p(ndims, ndims))
    for_each_component: DO ic = 1, ndims

          id = ic
          ! Slab type prerequisites
          array_of_sizes = SHAPE(velvel(1, ic)%values)
          DO id2 = 1, ndims
            array_of_subsizes_m(id2) = KronDelta(id, id2)*(velvel(1, ic)%b(id2, 1)-velvel(1, ic)%b_ol(id2, 1)) &
                                       + (1-KronDelta(id, id2))*(velvel(1, ic)%b(id2, 2)-velvel(1, ic)%b(id2, 1)+1)
            array_of_subsizes_p(id2) = KronDelta(id, id2)*(velvel(1, ic)%b_ol(id2, 2)-velvel(1, ic)%b(id2, 2)) &
                                       + (1-KronDelta(id, id2))*(velvel(1, ic)%b(id2, 2)-velvel(1, ic)%b(id2, 1)+1)
          END DO
          array_of_starts = [0, 0, 0]


          ! Creating slab types
          CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                      & array_of_sizes, array_of_subsizes_m, array_of_starts, &
                                      & MPI_ORDER_FORTRAN, &
                                      & MPI_DOUBLE_PRECISION, MPI_slab_m(ic,id), &
                                      & ierr)
          CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                      & array_of_sizes, array_of_subsizes_p, array_of_starts, &
                                      & MPI_ORDER_FORTRAN, &
                                      & MPI_DOUBLE_PRECISION, MPI_slab_p(ic,id), &
                                      & ierr)
          ! Committing types
          CALL MPI_TYPE_COMMIT(MPI_slab_m(ic,id), ierr)
          CALL MPI_TYPE_COMMIT(MPI_slab_p(ic,id), ierr)

          id = ic + 1
          IF (id>ndims) id = 1
          ! Slab type prerequisites
          array_of_sizes = SHAPE(velvel(2, ic)%values)
          DO id2 = 1, ndims
            array_of_subsizes_m(id2) = KronDelta(id, id2)*(velvel(2, ic)%b(id2, 1)-velvel(2, ic)%b_ol(id2, 1)) &
                                       + (1-KronDelta(id, id2))*(velvel(2, ic)%b(id2, 2)-velvel(2, ic)%b(id2, 1)+1)
            array_of_subsizes_p(id2) = KronDelta(id, id2)*(velvel(2, ic)%b_ol(id2, 2)-velvel(2, ic)%b(id2, 2)) &
                                       + (1-KronDelta(id, id2))*(velvel(2, ic)%b(id2, 2)-velvel(2, ic)%b(id2, 1)+1)
          END DO
          array_of_starts = [0, 0, 0]

          ! Creating slab types
          CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                      & array_of_sizes, array_of_subsizes_m, array_of_starts, &
                                      & MPI_ORDER_FORTRAN, &
                                      & MPI_DOUBLE_PRECISION, MPI_slab_m(ic,id), &
                                      & ierr)
          CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                      & array_of_sizes, array_of_subsizes_p, array_of_starts, &
                                      & MPI_ORDER_FORTRAN, &
                                      & MPI_DOUBLE_PRECISION, MPI_slab_p(ic,id), &
                                      & ierr)
          ! Committing types
          CALL MPI_TYPE_COMMIT(MPI_slab_m(ic,id), ierr)
          CALL MPI_TYPE_COMMIT(MPI_slab_p(ic,id), ierr)


          id = ic - 1
          IF (id<1) id = ndims
          ! Slab type prerequisites
          array_of_sizes = SHAPE(velvel(2, id)%values)
          DO id2 = 1, ndims
            array_of_subsizes_m(id2) = KronDelta(id, id2)*(velvel(2, id)%b(id2, 1)-velvel(2, id)%b_ol(id2, 1)) &
                                       + (1-KronDelta(id, id2))*(velvel(2, id)%b(id2, 2)-velvel(2, id)%b(id2, 1)+1)
            array_of_subsizes_p(id2) = KronDelta(id, id2)*(velvel(2, id)%b_ol(id2, 2)-velvel(2, id)%b(id2, 2)) &
                                       + (1-KronDelta(id, id2))*(velvel(2, id)%b(id2, 2)-velvel(2, id)%b(id2, 1)+1)
          END DO
          array_of_starts = [0, 0, 0]

          ! Creating slab types
          CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                      & array_of_sizes, array_of_subsizes_m, array_of_starts, &
                                      & MPI_ORDER_FORTRAN, &
                                      & MPI_DOUBLE_PRECISION, MPI_slab_m(ic,id), &
                                      & ierr)
          CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                      & array_of_sizes, array_of_subsizes_p, array_of_starts, &
                                      & MPI_ORDER_FORTRAN, &
                                      & MPI_DOUBLE_PRECISION, MPI_slab_p(ic,id), &
                                      & ierr)
          ! Committing types
          CALL MPI_TYPE_COMMIT(MPI_slab_m(ic,id), ierr)
          CALL MPI_TYPE_COMMIT(MPI_slab_p(ic,id), ierr)

    END DO for_each_component

  END SUBROUTINE set_conv_bounds



  SUBROUTINE conv_interp
    ! Subroutine that interpolates and multiply velocities to put into velvel block
    USE SPIKE
    USE essentials
    USE bandedmatrix
    USE Thomas_suite
    USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

    IMPLICIT NONE

    INTEGER :: ic, id, ider = 0, istag
    INTEGER :: j, k
    REAL, DIMENSION(:), ALLOCATABLE :: q
    REAL                            :: r = 0.0, NaN
    NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r


    !NOTE: homogeneus terms have additional layers only along the direction of their component
    ! i.e. u^2 only along x. This means that it is fine to give SPIKE_solve2 the whole 3d array
    ! using ":" for the remaining directions. The result is something like:
    !
    ! velvel(1, 1)%values(velvel(1, 1)%b(1, 1) : velvel(1, 1)%b(1, 2), :, :)
    !
    ! which translates into:
    !
    ! u^2(internal nodes, :, :)
    !
    ! Cross products instead have additional layers (for boundaries and overlaps) along two
    ! directions. THus it is safer to write something like:
    !
    ! velvel(2, 1)%values(velvel(2, 1)%b(1, 1) : velvel(2, 1)%b(1, 2), &
    !                     velvel(2, 1)%b(2, 1) : velvel(2, 1)%b(2, 2), &
    !                     velvel(2, 1)%b(3, 1) : velvel(2, 1)%b(3, 2))
    !
    ! which translates into:
    !
    ! uv(internal nodes, internal nodes, internal nodes)

      !!!!!!!!!! x interpolation !!!!!!!!!!
      id = 1
      ! homogeneus products
      ic = id
      istag = 2
      !velvel(1, id)%values = NaN
      ! FIXME cercare di evitare tutte queste allocazioni e deallocazioni
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(1, id)%b(2, 1), velvel(1, id)%b(2, 2)
        DO k = velvel(1, id)%b(3, 1), velvel(1, id)%b(3, 2)

          q = uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      velvel(1, id)%values(velvel(1, id)%b(id, 1) : velvel(1, id)%b(id, 2), j, k))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, velvel(1, id)%values(velvel(1, id)%b(id, 1) : velvel(1, id)%b(id, 2), :, :))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(1, id)%values(velvel(1, id)%b_bo(1, 1), &
                             velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                             velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 1), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(1, id)%values(velvel(1, id)%b_bo(1, 2), &
                             velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                             velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 2), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))
      END IF

      ! v on u nodes (for uv calculation)
      ic = 2
      istag = 1
      !intvel(id)%values = NaN
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = intvel(id)%b(2, 1), intvel(id)%b(2, 2)
        DO k = intvel(id)%b(3, 1), intvel(id)%b(3, 2)

          q = uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      intvel(id)%values(intvel(id)%b(id, 1) : intvel(id)%b(id, 2), j, k))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, intvel(id)%values(intvel(id)%b(id, 1) : intvel(id)%b(id, 2), &
                                                           intvel(id)%b(2, 1)  : intvel(id)%b(2, 2), &
                                                           intvel(id)%b(3, 1)  : intvel(id)%b(3, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        intvel(id)%values(intvel(id)%b_bo(1, 1), &
                          intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                          intvel(id)%b(3, 1):intvel(id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 1), &
                        intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                        intvel(id)%b(3, 1):intvel(id)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        intvel(id)%values(intvel(id)%b_bo(1, 2), &
                          intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                          intvel(id)%b(3, 1):intvel(id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 2), &
                        intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                        intvel(id)%b(3, 1):intvel(id)%b(3, 2))
      END IF

      ! w on u nodes (for wu calculation)
      ic = 3
      istag = 1
      !velvel(2, ic)%values = NaN
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(2, ic)%b(2, 1), velvel(2, ic)%b(2, 2)
        DO k = velvel(2, ic)%b(3, 1), velvel(2, ic)%b(3, 2)

          q = uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      velvel(2, ic)%values(velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2), j, k))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, velvel(2, ic)%values(velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2), &
                                                              velvel(2, ic)%b(2, 1)  : velvel(2, ic)%b(2, 2), &
                                                              velvel(2, ic)%b(3, 1)  : velvel(2, ic)%b(3, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(2, ic)%values(velvel(2, ic)%b_bo(1, 1), &
                          velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                          velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 1), &
                        velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                        velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(2, ic)%values(velvel(2, ic)%b_bo(1, 2), &
                          velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                          velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 2), &
                        velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                        velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2))
      END IF



      !!!!!!!!!! y interpolation !!!!!!!!!!
      id = 2
      ! homogeneus products
      ic = id
      istag = 2
      !velvel(1, id)%values = NaN
      ! FIXME cercare di evitare tutte queste allocazioni e deallocazioni
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(1, id)%b(1, 1), velvel(1, id)%b(1, 2)
        DO k = velvel(1, id)%b(3, 1), velvel(1, id)%b(3, 2)

          q = uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      velvel(1, id)%values(j, velvel(1, id)%b(id, 1) :velvel(1, id)%b(id, 2), k))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, velvel(1, id)%values(:, velvel(1, id)%b(id, 1) : velvel(1, id)%b(id, 2), :))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(1, id)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                             velvel(1, id)%b_bo(2, 1), &
                             velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2)) = &
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 1), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(1, id)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                             velvel(1, id)%b_bo(2, 2), &
                             velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2)) = &
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 2), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))
      END IF

      ! w on v nodes (for vw calculation)
      ic = 3
      istag = 1
      !intvel(id)%values = NaN
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = intvel(id)%b(1, 1), intvel(id)%b(1, 2)
        DO k = intvel(id)%b(3, 1), intvel(id)%b(3, 2)

          q = uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      intvel(id)%values(j, intvel(id)%b(id, 1) : intvel(id)%b(id, 2), k))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, intvel(id)%values(intvel(id)%b(1, 1)  : intvel(id)%b(1, 2), &
                                                           intvel(id)%b(id, 1) : intvel(id)%b(id, 2), &
                                                           intvel(id)%b(3, 1)  : intvel(id)%b(3, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        intvel(id)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                             intvel(id)%b_bo(2, 1), &
                             intvel(id)%b(3, 1):intvel(id)%b(3, 2)) = &
        uvwp(ic)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 1), &
                        intvel(id)%b(3, 1):intvel(id)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        intvel(id)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                             intvel(id)%b_bo(2, 2), &
                             intvel(id)%b(3, 1):intvel(id)%b(3, 2)) = &
        uvwp(ic)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 2), &
                        intvel(id)%b(3, 1):intvel(id)%b(3, 2))
      END IF

      ! u on v nodes (for uv calculation)
      ic = 1
      istag = 1
      !velvel(2, ic)%values = NaN
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(2, ic)%b(1, 1), velvel(2, ic)%b(1, 2)
        DO k = velvel(2, ic)%b(3, 1), velvel(2, ic)%b(3, 2)

          q = uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      velvel(2, ic)%values(j, velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2), k))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, velvel(2, ic)%values(velvel(2, ic)%b(1, 1)  : velvel(2, ic)%b(1, 2), &
                                                              velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2), &
                                                              velvel(2, ic)%b(3, 1)  : velvel(2, ic)%b(3, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(2, ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                             velvel(2, ic)%b_bo(2, 1), &
                             velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2)) = &
        uvwp(ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 1), &
                        velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(2, ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                             velvel(2, ic)%b_bo(2, 2), &
                             velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2)) = &
        uvwp(ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 2), &
                        velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2))
      END IF



      !!!!!!!!!! z interpolation !!!!!!!!!!
      id = 3
      ! homogeneus products
      ic = id
      istag = 2
      !velvel(1, id)%values = NaN
      ! FIXME cercare di evitare tutte queste allocazioni e deallocazioni
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(1, id)%b(1, 1), velvel(1, id)%b(1, 2)
        DO k = velvel(1, id)%b(2, 1), velvel(1, id)%b(2, 2)

          q = uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      velvel(1, id)%values(j, k, velvel(1, id)%b(id, 1) : velvel(1, id)%b(id, 2)))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, velvel(1, id)%values(:, :, velvel(1, id)%b(id, 1) : velvel(1, id)%b(id, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(1, id)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                             velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                             velvel(1, id)%b_bo(3, 1)) = &
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 1))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(1, id)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                             velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                             velvel(1, id)%b_bo(3, 2)) = &
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 2))
      END IF

      ! u on w nodes (for wu calculation)
      ic = 1
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = intvel(id)%b(1, 1), intvel(id)%b(1, 2)
        DO k = intvel(id)%b(2, 1), intvel(id)%b(2, 2)

          q = uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      intvel(id)%values(j, k, intvel(id)%b(id, 1) : intvel(id)%b(id, 2)))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, intvel(id)%values(intvel(id)%b(1, 1)  : intvel(id)%b(1, 2), &
                                                           intvel(id)%b(2, 1)  : intvel(id)%b(2, 2), &
                                                           intvel(id)%b(id, 1) : intvel(id)%b(id, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        intvel(id)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                             intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                             intvel(id)%b_bo(3, 1)) = &
        uvwp(ic)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                        intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 1))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        intvel(id)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                             intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                             intvel(id)%b_bo(3, 2)) = &
        uvwp(ic)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                        intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 2))
      END IF

      ! v on w nodes (for vw calculation)
      ic = 2
      istag = 1
      !velvel(2, ic)%values = NaN
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(2, ic)%b(1, 1), velvel(2, ic)%b(1, 2)
        DO k = velvel(2, ic)%b(2, 1), velvel(2, ic)%b(2, 2)

          q = uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))
          CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                      cmp(id, ider, istag)%B*q     , &
                      velvel(2, ic)%values(j, k, velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2)))

        END DO
      END DO
      DEALLOCATE(q)
      CALL SPIKE_solve2(id, ider, istag, velvel(2, ic)%values(velvel(2, ic)%b(1, 1) : velvel(2, ic)%b(1, 2), &
                                                              velvel(2, ic)%b(2, 1) : velvel(2, ic)%b(2, 2), &
                                                              velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2)))
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(2, ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                             velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                             velvel(2, ic)%b_bo(3, 1)) = &
        uvwp(ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                        velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 1))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(2, ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                             velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                             velvel(2, ic)%b_bo(3, 2)) = &
        uvwp(ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                        velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 2))
      END IF


  END SUBROUTINE conv_interp



  SUBROUTINE conv_exchange

    USE MPI,        ONLY : MPI_COMM_WORLD, MPI_STATUS_SIZE, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_PROC_NULL
    USE essentials

    IMPLICIT NONE

    INTEGER :: ic, id, ierr
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

    !!!!!!!!!! Homogeneus products !!!!!!!!!!
    id = 1
    ic = id
    ! x forward communication
    CALL MPI_SENDRECV(velvel(1, id)%values( &
                      velvel(1, id)%b(id, 2) + velvel(1, id)%b_ol(id, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 182, &
                      velvel(1, id)%values( &
                      velvel(1, id)%b_bo(id, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 182, &
                      procs_grid, status, ierr)
    ! x backward communication
    CALL MPI_SENDRECV(velvel(1, id)%values( &
                      velvel(1, id)%b(id, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 183, &
                      velvel(1, id)%values( &
                      velvel(1, id)%b(id, 2)+1, &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 183, &
                      procs_grid, status, ierr)

    id = 2
    ic = id
    ! y forward communication
    CALL MPI_SENDRECV(velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(id, 2) + velvel(1, id)%b_ol(id, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 184, &
                      velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b_bo(id, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 184, &
                      procs_grid, status, ierr)
    ! y backward communication
    CALL MPI_SENDRECV(velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(id, 1), &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 185, &
                      velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(id, 2)+1, &
                      velvel(1, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 185, &
                      procs_grid, status, ierr)

    id = 3
    ic = id
    ! z forward communication
    CALL MPI_SENDRECV(velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(id, 2) + velvel(1, id)%b_ol(id, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 186, &
                      velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b_bo(id, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 186, &
                      procs_grid, status, ierr)
    ! z backward communication
    CALL MPI_SENDRECV(velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(id, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 187, &
                      velvel(1, id)%values( &
                      velvel(1, id)%b(1, 1), &
                      velvel(1, id)%b(2, 1), &
                      velvel(1, id)%b(id, 2)+1), &
                      1, MPI_slab_p(ic, id), idp(id), 187, &
                      procs_grid, status, ierr)


    !!!!!!!!!! Cross products !!!!!!!!!!
    ic = 1

    !!!!! uv along y !!!!!
    id = 2
    ! y forward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(id, 2) + velvel(2, ic)%b_ol(id, 1)-velvel(2, ic)%b(id, 1)+1, &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 188, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b_bo(id, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 188, &
                      procs_grid, status, ierr)
    ! y backward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(id, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 189, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(id, 2)+1, &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 189, &
                      procs_grid, status, ierr)

    !!!!! wu along z !!!!!
    id = 3
    ! z forward communication
    CALL MPI_SENDRECV(velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(id, 2) + velvel(2, id)%b_ol(id, 1)-velvel(2, id)%b(id, 1)+1), &
                      1, MPI_slab_m(ic, id), idp(id), 190, &
                      velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b_bo(id, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 190, &
                      procs_grid, status, ierr)
    ! z backward communication
    CALL MPI_SENDRECV(velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(id, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 191, &
                      velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(id, 2)+1), &
                      1, MPI_slab_p(ic, id), idp(id), 191, &
                      procs_grid, status, ierr)

    ic = 2

    !!!!! vw along z !!!!!
    id = 3
    ! z forward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(id, 2) + velvel(2, ic)%b_ol(id, 1)-velvel(2, ic)%b(id, 1)+1), &
                      1, MPI_slab_m(ic, id), idp(id), 192, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b_bo(id, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 192, &
                      procs_grid, status, ierr)
    ! z backward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(id, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 193, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(id, 2)+1), &
                      1, MPI_slab_p(ic, id), idp(id), 193, &
                      procs_grid, status, ierr)

    !!!!! uv along x !!!!!
    id = 1
    ! x forward communication
    CALL MPI_SENDRECV(velvel(2, id)%values( &
                      velvel(2, id)%b(id, 2) + velvel(2, id)%b_ol(id, 1)-velvel(2, id)%b(id, 1)+1, &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 194, &
                      velvel(2, id)%values( &
                      velvel(2, id)%b_bo(id, 1), &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 194, &
                      procs_grid, status, ierr)
    ! x backward communication
    CALL MPI_SENDRECV(velvel(2, id)%values( &
                      velvel(2, id)%b(id, 1), &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 195, &
                      velvel(2, id)%values( &
                      velvel(2, id)%b(id, 2)+1, &
                      velvel(2, id)%b(2, 1), &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 195, &
                      procs_grid, status, ierr)

    ic = 3

    !!!!! wu along x !!!!!
    id = 1
    ! x forward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(id, 2) + velvel(2, ic)%b_ol(id, 1)-velvel(2, ic)%b(id, 1)+1, &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 196, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b_bo(id, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 196, &
                      procs_grid, status, ierr)
    ! x backward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(id, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 197, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(id, 2)+1, &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 197, &
                      procs_grid, status, ierr)

    !!!!! vw along y !!!!!
    id = 2
    ! y forward communication
    CALL MPI_SENDRECV(velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(id, 2) + velvel(2, id)%b_ol(id, 1)-velvel(2, id)%b(id, 1)+1, &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 198, &
                      velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b_bo(id, 1), &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 198, &
                      procs_grid, status, ierr)
    ! y backward communication
    CALL MPI_SENDRECV(velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(id, 1), &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 199, &
                      velvel(2, id)%values( &
                      velvel(2, id)%b(1, 1), &
                      velvel(2, id)%b(id, 2)+1, &
                      velvel(2, id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 199, &
                      procs_grid, status, ierr)


  END SUBROUTINE conv_exchange



  SUBROUTINE conv_calc

    ! Subroutine that derives previously calculated squared velocities and stores
    ! the results in the convvel block
    USE SPIKE
    USE essentials
    USE bandedmatrix
    USE Thomas_suite
    USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

    IMPLICIT NONE

    INTEGER                               :: ic, id, ider = 1, istag
    INTEGER                               :: j, k, iii
    REAL, DIMENSION(:), ALLOCATABLE       :: q, psi  ! due to ANOMALIA a solution vector has to be defined
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: temp   ! used to temporarily store block-diagonal solution
    REAL                                  :: r = 0.0, NaN  ! r is a dummy real used to define a NaN of type real
    NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

    !!!!!!!!!! Velocity interpolations !!!!!!!!!!
    CALL conv_interp


    !!!!!!!!!! Velocity multiplication !!!!!!!!!!
    DO ic = 1, ndims
      velvel(1, ic)%values = velvel(1, ic)%values*velvel(1, ic)%values
      velvel(2, ic)%values = velvel(2, ic)%values*intvel(ic)%values
    END DO


    !!!!!!!!!! Velocity*velocity communications !!!!!!!!!!
    CALL conv_exchange


    !!!!!!!!!! Debugging inizialization !!!!!!!!!!
    DO ic = 1, ndims
      convvel(ic)%values = NaN
    END DO

    !!!!!!!!!! Velocity*velocity derivation !!!!!!!!!!

                    !!!!! x component !!!!!

    ! d(u^2)/dx
    id = 1
    istag = 1
    ALLOCATE(q(velvel(1, 1)%b_bo(id, 1):velvel(1, 1)%b_bo(id, 2)))
    ALLOCATE(temp(cmp(id, ider, istag)%A%ub(2)-cmp(id, ider, istag)%A%lb(2)+1, &
                  convvel(1)%b(2, 1):convvel(1)%b(2, 2), &
                  convvel(1)%b(3, 1):convvel(1)%b(3, 2)))
    DO j = convvel(1)%b(2, 1), convvel(1)%b(2, 2)
      DO k = convvel(1)%b(3, 1), convvel(1)%b(3, 2)

        q = velvel(1, 1)%values(:, j, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(:, j, k))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(1)%values = temp(convvel(1)%b(id, 1):convvel(1)%b(id, 2), :, :)
    DEALLOCATE(q)
    DEALLOCATE(temp)


    ! d(uv)/dy
    id = 2
    istag = 2
    ALLOCATE(q(velvel(2, 1)%b_bo(id, 1):velvel(2, 1)%b_bo(id, 2)))
    ALLOCATE(temp(convvel(1)%b(1, 1):convvel(1)%b(1, 2), &
                  convvel(1)%b(2, 1):convvel(1)%b(2, 2), &
                  convvel(1)%b(3, 1):convvel(1)%b(3, 2)))
    temp = NaN
    DO j = convvel(1)%b(1, 1), convvel(1)%b(1, 2)
      DO k = convvel(1)%b(3, 1), convvel(1)%b(3, 2)

        q = velvel(2, 1)%values(j, :, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, :, k))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(1)%values = convvel(1)%values + temp
    DEALLOCATE(q)


    ! d(wu)/dz
    id = 3
    istag = 2
    ALLOCATE(q(velvel(2, 3)%b_bo(id, 1):velvel(2, 3)%b_bo(id, 2)))
    DO j = convvel(1)%b(1, 1), convvel(1)%b(1, 2)
      DO k = convvel(1)%b(2, 1), convvel(1)%b(2, 2)

        q = velvel(2, 3)%values(j, k, :)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, k, :))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(1)%values = convvel(1)%values + temp
    DEALLOCATE(q)
    DEALLOCATE(temp)


                      !!!!! y component !!!!!

    ! d(v^2)/dy
    id = 2
    istag = 1
    ALLOCATE(q(velvel(1, 2)%b_bo(id, 1):velvel(1, 2)%b_bo(id, 2)))
    ALLOCATE(temp(convvel(2)%b(1, 1):convvel(2)%b(1, 2), &
                  cmp(id, ider, istag)%A%ub(2)-cmp(id, ider, istag)%A%lb(2)+1, &
                  convvel(2)%b(3, 1):convvel(2)%b(3, 2)))
    DO j = convvel(2)%b(1, 1), convvel(2)%b(1, 2)
      DO k = convvel(2)%b(3, 1), convvel(2)%b(3, 2)

        q = velvel(1, 2)%values(j, :, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, :, k))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(2)%values = temp(:, convvel(2)%b(id, 1):convvel(2)%b(id, 2), :)
    DEALLOCATE(q)
    DEALLOCATE(temp)


    ! d(uv)/dx
    id = 1
    istag = 2
    ALLOCATE(q(velvel(2, 1)%b_bo(id, 1):velvel(2, 1)%b_bo(id, 2)))
    ALLOCATE(temp(convvel(2)%b(1, 1):convvel(2)%b(1, 2), &
                  convvel(2)%b(2, 1):convvel(2)%b(2, 2), &
                  convvel(2)%b(3, 1):convvel(2)%b(3, 2)))
    DO j = convvel(2)%b(2, 1), convvel(2)%b(2, 2)
      DO k = convvel(2)%b(3, 1), convvel(2)%b(3, 2)

        q = velvel(2, 1)%values(:, j, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(:, j, k))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(2)%values = convvel(2)%values + temp
    DEALLOCATE(q)


    ! d(vw)/dz
    id = 3
    istag = 2
    ALLOCATE(q(velvel(2, 2)%b_bo(id, 1):velvel(2, 2)%b_bo(id, 2)))
    DO j = convvel(2)%b(1, 1), convvel(2)%b(1, 2)
      DO k = convvel(2)%b(2, 1), convvel(2)%b(2, 2)

        q = velvel(2, 2)%values(j, k, :)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, k, :))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(2)%values = convvel(2)%values + temp
    DEALLOCATE(q)
    DEALLOCATE(temp)


                    !!!!! z component !!!!!

    ! d(w^2)/dz
    id = 3
    istag = 1
    ALLOCATE(q(velvel(1, 3)%b_bo(id, 1):velvel(1, 3)%b_bo(id, 2)))
    ALLOCATE(temp(convvel(3)%b(1, 1):convvel(3)%b(1, 2), &
                  convvel(3)%b(2, 1):convvel(3)%b(2, 2), &
                  cmp(id, ider, istag)%A%ub(2)-cmp(id, ider, istag)%A%lb(2)+1))
    DO j = convvel(3)%b(1, 1), convvel(3)%b(1, 2)
      DO k = convvel(3)%b(2, 1), convvel(3)%b(2, 2)

        q = velvel(1, 3)%values(j, k, :)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, k, :))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(3)%values = temp(:, :, convvel(3)%b(id, 1):convvel(3)%b(id, 2))
    DEALLOCATE(q)
    DEALLOCATE(temp)


    ! d(wu)/dx
    id = 1
    istag = 2
    ALLOCATE(q(velvel(2, 3)%b_bo(id, 1):velvel(2, 3)%b_bo(id, 2)))
    ALLOCATE(temp(convvel(3)%b(1, 1):convvel(3)%b(1, 2), &
                  convvel(3)%b(2, 1):convvel(3)%b(2, 2), &
                  convvel(3)%b(3, 1):convvel(3)%b(3, 2)))
    DO j = convvel(3)%b(2, 1), convvel(3)%b(2, 2)
      DO k = convvel(3)%b(3, 1), convvel(3)%b(3, 2)

        q = velvel(2, 3)%values(:, j, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(:, j, k))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(3)%values = convvel(3)%values + temp
    DEALLOCATE(q)


    ! d(vw)/dy
    id = 2
    istag = 2
    ALLOCATE(q(velvel(2, 2)%b_bo(id, 1):velvel(2, 2)%b_bo(id, 2)))
    DO j = convvel(3)%b(1, 1), convvel(3)%b(1, 2)
      DO k = convvel(3)%b(3, 1), convvel(3)%b(3, 2)

        q = velvel(2, 2)%values(j, :, k)
        CALL Thomas(cmp(id, ider, istag)%A%matrix, &
                    cmp(id, ider, istag)%B*q     , &
                    temp(j, :, k))

      END DO
    END DO
    CALL SPIKE_solve2(id, ider, istag, temp)
    convvel(3)%values = convvel(3)%values + temp
    DEALLOCATE(q)
    DEALLOCATE(temp)


  END SUBROUTINE conv_calc



END MODULE convective_term
