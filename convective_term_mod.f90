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
          !IF (myid==0 .AND. id==1) THEN
          !  PRINT *, array_of_subsizes_m
          !  PRINT *, array_of_subsizes_p
          !END IF

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
                                       + (1-KronDelta(id, id2))*(velvel(2, ic)%b_bo(id2, 2)-velvel(2, ic)%b_bo(id2, 1)+1)
            array_of_subsizes_p(id2) = KronDelta(id, id2)*(velvel(2, ic)%b_ol(id2, 2)-velvel(2, ic)%b(id2, 2)) &
                                       + (1-KronDelta(id, id2))*(velvel(2, ic)%b_bo(id2, 2)-velvel(2, ic)%b_bo(id2, 1)+1)
          END DO
          array_of_starts = [0, 0, 0]
          !IF (myid==0 .AND. id==1) THEN
          !  PRINT *, array_of_subsizes_m
          !  PRINT *, array_of_subsizes_p
          !END IF

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
          array_of_sizes = SHAPE(intvel(id)%values)
          DO id2 = 1, ndims
            array_of_subsizes_m(id2) = KronDelta(id, id2)*(intvel(id)%b(id2, 1)-intvel(id)%b_ol(id2, 1)) &
                                       + (1-KronDelta(id, id2))*(intvel(id)%b_bo(id2, 2)-intvel(id)%b_bo(id2, 1)+1)
            array_of_subsizes_p(id2) = KronDelta(id, id2)*(intvel(id)%b_ol(id2, 2)-intvel(id)%b(id2, 2)) &
                                       + (1-KronDelta(id, id2))*(intvel(id)%b_bo(id2, 2)-intvel(id)%b_bo(id2, 1)+1)
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

    IMPLICIT NONE

    INTEGER :: ic, id, ider = 0, istag
    INTEGER :: j, k
    REAL, DIMENSION(:), ALLOCATABLE :: q

    !for_each_component: DO ic = 1, nn

      !!!!!!!!!! x interpolation !!!!!!!!!!
      id = 1
      ! homogeneus products
      ic = id
      istag = 2
      ! FIXME cercare di evitare tutte queste allocazioni e deallocazioni
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(1, id)%b(2, 1), velvel(1, id)%b(2, 2)
        DO k = velvel(1, id)%b(3, 1), velvel(1, id)%b(3, 2)
          q = uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)*&
              uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)    ! product is faster than power, even if less readable...
          CALL SPIKE_solve(id, ider, istag, &
                           velvel(1, id)%values(velvel(1, id)%b(id, 1) : velvel(1, id)%b(id, 2), j, k), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(1, id)%values(velvel(1, id)%b_bo(1, 1), &
                             velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                             velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 1), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))*&
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
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))*&
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 2), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))
      END IF

      ! v on u nodes (for uv calculation)
      ic = 2
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = intvel(id)%b_bo(2, 1), intvel(id)%b_bo(2, 2)
        DO k = intvel(id)%b(3, 1), intvel(id)%b(3, 2)
          q = uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)
          CALL SPIKE_solve(id, ider, istag, &
                           intvel(id)%values(intvel(id)%b(id, 1) : intvel(id)%b(id, 2), j, k), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        intvel(id)%values(intvel(id)%b_bo(1, 1), &
                          intvel(id)%b_bo(2, 1):intvel(id)%b_bo(2, 2), &
                          intvel(id)%b(3, 1):intvel(id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 1), &
                        intvel(id)%b_bo(2, 1):intvel(id)%b_bo(2, 2), &
                        intvel(id)%b(3, 1):intvel(id)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        intvel(id)%values(intvel(id)%b_bo(1, 2), &
                          intvel(id)%b_bo(2, 1):intvel(id)%b_bo(2, 2), &
                          intvel(id)%b(3, 1):intvel(id)%b(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 2), &
                        intvel(id)%b_bo(2, 1):intvel(id)%b_bo(2, 2), &
                        intvel(id)%b(3, 1):intvel(id)%b(3, 2))
      END IF

      ! w on u nodes (for wu calculation)
      ic = 3
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(2, ic)%b(2, 1), velvel(2, ic)%b(2, 2)
        DO k = velvel(2, ic)%b_bo(3, 1), velvel(2, ic)%b_bo(3, 2)
          q = uvwp(ic)%values(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), j, k)
          CALL SPIKE_solve(id, ider, istag, &
                           velvel(2, ic)%values(velvel(2, ic)%b(id, 1) : velvel(2, ic)%b(id, 2), j, k), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(2, ic)%values(velvel(2, ic)%b_bo(1, 1), &
                          velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                          velvel(2, ic)%b_bo(3, 1):velvel(2, ic)%b_bo(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 1), &
                        velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                        velvel(2, ic)%b_bo(3, 1):velvel(2, ic)%b_bo(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(2, ic)%values(velvel(2, ic)%b_bo(1, 2), &
                          velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                          velvel(2, ic)%b_bo(3, 1):velvel(2, ic)%b_bo(3, 2)) = &
        uvwp(ic)%values(uvwp(ic)%b_bo(1, 2), &
                        velvel(2, ic)%b(2, 1):velvel(2, ic)%b(2, 2), &
                        velvel(2, ic)%b_bo(3, 1):velvel(2, ic)%b_bo(3, 2))
      END IF



      !!!!!!!!!! y interpolation !!!!!!!!!!
      id = 2
      ! homogeneus products
      ic = id
      istag = 2
      ! FIXME cercare di evitare tutte queste allocazioni e deallocazioni
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(1, id)%b(1, 1), velvel(1, id)%b(1, 2)
        DO k = velvel(1, id)%b(3, 1), velvel(1, id)%b(3, 2)
          q = uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)*&
              uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)    ! product is faster than power, even if less readable...
          CALL SPIKE_solve(id, ider, istag, &
                           velvel(1, id)%values(j, velvel(1, id)%b(id, 1) :velvel(1, id)%b(id, 2), k), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(1, id)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                             velvel(1, id)%b_bo(2, 1), &
                             velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2)) = &
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 1), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))*&
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
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))*&
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 2), &
                        velvel(1, id)%b(3, 1):velvel(1, id)%b(3, 2))
      END IF

      ! w on v nodes (for vw calculation)
      ic = 3
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = intvel(id)%b(1, 1), intvel(id)%b(1, 2)
        DO k = intvel(id)%b_bo(3, 1), intvel(id)%b_bo(3, 2)
          q = uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)
          CALL SPIKE_solve(id, ider, istag, &
                           intvel(id)%values(j, intvel(id)%b(id, 1) : intvel(id)%b(id, 2), k), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        intvel(id)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                             intvel(id)%b_bo(2, 1), &
                             intvel(id)%b_bo(3, 1):intvel(id)%b_bo(3, 2)) = &
        uvwp(ic)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 1), &
                        intvel(id)%b_bo(3, 1):intvel(id)%b_bo(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        intvel(id)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                             intvel(id)%b_bo(2, 2), &
                             intvel(id)%b_bo(3, 1):intvel(id)%b_bo(3, 2)) = &
        uvwp(ic)%values(intvel(id)%b(1, 1):intvel(id)%b(1, 2), &
                        uvwp(ic)%b_bo(2, 2), &
                        intvel(id)%b_bo(3, 1):intvel(id)%b_bo(3, 2))
      END IF

      ! u on v nodes (for uv calculation)
      ic = 1
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(2, ic)%b_bo(1, 1), velvel(2, ic)%b_bo(1, 2)
        DO k = velvel(2, ic)%b(3, 1), velvel(2, ic)%b(3, 2)
          q = uvwp(ic)%values(j, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2), k)
          CALL SPIKE_solve(id, ider, istag, &
                           velvel(2, ic)%values(j, velvel(2, ic)%b(id, 1) :velvel(2, ic)%b(id, 2), k), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(2, ic)%values(velvel(2, ic)%b_bo(1, 1):velvel(2, ic)%b_bo(1, 2), &
                             velvel(2, ic)%b_bo(2, 1), &
                             velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2)) = &
        uvwp(ic)%values(velvel(2, ic)%b_bo(1, 1):velvel(2, ic)%b_bo(1, 2), &
                        uvwp(ic)%b_bo(2, 1), &
                        velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(2, ic)%values(velvel(2, ic)%b_bo(1, 1):velvel(2, ic)%b_bo(1, 2), &
                             velvel(2, ic)%b_bo(2, 2), &
                             velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2)) = &
        uvwp(ic)%values(velvel(2, ic)%b_bo(1, 1):velvel(2, ic)%b_bo(1, 2), &
                        uvwp(ic)%b_bo(2, 2), &
                        velvel(2, ic)%b(3, 1):velvel(2, ic)%b(3, 2))
      END IF



      !!!!!!!!!! z interpolation !!!!!!!!!!
      id = 3
      ! homogeneus products
      ic = id
      istag = 2
      ! FIXME cercare di evitare tutte queste allocazioni e deallocazioni
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(1, id)%b(1, 1), velvel(1, id)%b(1, 2)
        DO k = velvel(1, id)%b(2, 1), velvel(1, id)%b(2, 2)
          q = uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))*&
              uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))    ! product is faster than power, even if less readable...
          CALL SPIKE_solve(id, ider, istag, &
                           velvel(1, id)%values(j, k, velvel(1, id)%b(id, 1) :velvel(1, id)%b(id, 2)), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(1, id)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                             velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                             velvel(1, id)%b_bo(3, 1)) = &
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 1))*&
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
                        uvwp(ic)%b_bo(3, 2))*&
        uvwp(ic)%values(velvel(1, id)%b(1, 1):velvel(1, id)%b(1, 2), &
                        velvel(1, id)%b(2, 1):velvel(1, id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 2))
      END IF

      ! u on w nodes (for wu calculation)
      ic = 1
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = intvel(id)%b_bo(1, 1), intvel(id)%b_bo(1, 2)
        DO k = intvel(id)%b(2, 1), intvel(id)%b(2, 2)
          q = uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))
          CALL SPIKE_solve(id, ider, istag, &
                           intvel(id)%values(j, k, intvel(id)%b(id, 1) :intvel(id)%b(id, 2)), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        intvel(id)%values(intvel(id)%b_bo(1, 1):intvel(id)%b_bo(1, 2), &
                             intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                             intvel(id)%b_bo(3, 1)) = &
        uvwp(ic)%values(intvel(id)%b_bo(1, 1):intvel(id)%b_bo(1, 2), &
                        intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 1))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        intvel(id)%values(intvel(id)%b_bo(1, 1):intvel(id)%b_bo(1, 2), &
                             intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                             intvel(id)%b_bo(3, 2)) = &
        uvwp(ic)%values(intvel(id)%b_bo(1, 1):intvel(id)%b_bo(1, 2), &
                        intvel(id)%b(2, 1):intvel(id)%b(2, 2), &
                        uvwp(ic)%b_bo(3, 2))
      END IF

      ! v on w nodes (for vw calculation)
      ic = 2
      istag = 1
      ALLOCATE(q(cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2)))
      DO j = velvel(2, ic)%b(1, 1), velvel(2, ic)%b(1, 2)
        DO k = velvel(2, ic)%b_bo(2, 1), velvel(2, ic)%b_bo(2, 2)
          q = uvwp(ic)%values(j, k, cmp(id, ider, istag)%B%lb(2):cmp(id, ider, istag)%B%ub(2))
          CALL SPIKE_solve(id, ider, istag, &
                           velvel(2, ic)%values(j, k, velvel(2, ic)%b(id, 1) :velvel(2, ic)%b(id, 2)), &
                           q)
        END DO
      END DO
      DEALLOCATE(q)
      ! need to add boundary values
      IF (mycoords(id)==0) THEN
        velvel(2, ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                             velvel(2, ic)%b_bo(2, 1):velvel(2, ic)%b_bo(2, 2), &
                             velvel(2, ic)%b_bo(3, 1)) = &
        uvwp(ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                        velvel(2, ic)%b_bo(2, 1):velvel(2, ic)%b_bo(2, 2), &
                        uvwp(ic)%b_bo(3, 1))
      END IF
      IF (mycoords(id)==dims(id)-1) THEN
        velvel(2, ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                             velvel(2, ic)%b_bo(2, 1):velvel(2, ic)%b_bo(2, 2), &
                             velvel(2, ic)%b_bo(3, 2)) = &
        uvwp(ic)%values(velvel(2, ic)%b(1, 1):velvel(2, ic)%b(1, 2), &
                        velvel(2, ic)%b_bo(2, 1):velvel(2, ic)%b_bo(2, 2), &
                        uvwp(ic)%b_bo(3, 2))
      END IF



      !!!!!!!!!! Exchange overlap values !!!!!!!!!!
      !velvel(2, 1)%values = myid
      CALL conv_exchange


    !END DO for_each_component

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
    id = 2
    ! y forward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b_bo(1, 1), &
                      velvel(2, ic)%b(id, 2) + velvel(2, ic)%b_ol(id, 1)-velvel(2, ic)%b(id, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 188, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b_bo(1, 1), &
                      velvel(2, ic)%b_bo(id, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 188, &
                      procs_grid, status, ierr)
    ! y backward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b_bo(1, 1), &
                      velvel(2, ic)%b(id, 1), &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 189, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b_bo(1, 1), &
                      velvel(2, ic)%b(id, 2)+1, &
                      velvel(2, ic)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 189, &
                      procs_grid, status, ierr)
    id = 3
    ! z forward communication
    CALL MPI_SENDRECV(intvel(id)%values( &
                      intvel(id)%b_bo(1, 1), &
                      intvel(id)%b(2, 1), &
                      intvel(id)%b(id, 2) + intvel(id)%b_ol(id, 1)-intvel(id)%b(id, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 190, &
                      intvel(id)%values( &
                      intvel(id)%b_bo(1, 1), &
                      intvel(id)%b(2, 1), &
                      intvel(id)%b_bo(id, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 190, &
                      procs_grid, status, ierr)
    ! z backward communication
    CALL MPI_SENDRECV(intvel(id)%values( &
                      intvel(id)%b_bo(1, 1), &
                      intvel(id)%b(2, 1), &
                      intvel(id)%b(id, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 191, &
                      intvel(id)%values( &
                      intvel(id)%b_bo(1, 1), &
                      intvel(id)%b(2, 1), &
                      intvel(id)%b(id, 2)+1), &
                      1, MPI_slab_p(ic, id), idp(id), 191, &
                      procs_grid, status, ierr)

    ic = 2
    id = 3
    ! z forward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b_bo(2, 1), &
                      velvel(2, ic)%b(id, 2) + velvel(2, ic)%b_ol(id, 1)-velvel(2, ic)%b(id, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 192, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b_bo(2, 1), &
                      velvel(2, ic)%b_bo(id, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 192, &
                      procs_grid, status, ierr)
    ! z backward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b_bo(2, 1), &
                      velvel(2, ic)%b(id, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 193, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(1, 1), &
                      velvel(2, ic)%b_bo(2, 1), &
                      velvel(2, ic)%b(id, 2)+1), &
                      1, MPI_slab_p(ic, id), idp(id), 193, &
                      procs_grid, status, ierr)
    id = 1
    ! x forward communication
    CALL MPI_SENDRECV(intvel(id)%values( &
                      intvel(id)%b(id, 2) + intvel(id)%b_ol(id, 1)-intvel(id)%b(id, 1), &
                      intvel(id)%b_bo(2, 1), &
                      intvel(id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 194, &
                      intvel(id)%values( &
                      intvel(id)%b_bo(id, 1), &
                      intvel(id)%b_bo(2, 1), &
                      intvel(id)%b(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 194, &
                      procs_grid, status, ierr)
    ! x backward communication
    CALL MPI_SENDRECV(intvel(id)%values( &
                      intvel(id)%b(id, 1), &
                      intvel(id)%b_bo(2, 1), &
                      intvel(id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 195, &
                      intvel(id)%values( &
                      intvel(id)%b(id, 2)+1, &
                      intvel(id)%b_bo(2, 1), &
                      intvel(id)%b(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 195, &
                      procs_grid, status, ierr)

    ic = 3
    id = 1
    ! x forward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(id, 2) + velvel(2, ic)%b_ol(id, 1)-velvel(2, ic)%b(id, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b_bo(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 196, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b_bo(id, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b_bo(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 196, &
                      procs_grid, status, ierr)
    ! x backward communication
    CALL MPI_SENDRECV(velvel(2, ic)%values( &
                      velvel(2, ic)%b(id, 1), &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b_bo(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 197, &
                      velvel(2, ic)%values( &
                      velvel(2, ic)%b(id, 2)+1, &
                      velvel(2, ic)%b(2, 1), &
                      velvel(2, ic)%b_bo(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 197, &
                      procs_grid, status, ierr)
    id = 2
    ! y forward communication
    CALL MPI_SENDRECV(intvel(id)%values( &
                      intvel(id)%b(1, 1), &
                      intvel(id)%b(id, 2) + intvel(id)%b_ol(id, 1)-intvel(id)%b(id, 1), &
                      intvel(id)%b_bo(3, 1)), &
                      1, MPI_slab_m(ic, id), idp(id), 198, &
                      intvel(id)%values( &
                      intvel(id)%b(1, 1), &
                      intvel(id)%b_bo(id, 1), &
                      intvel(id)%b_bo(3, 1)), &
                      1, MPI_slab_m(ic, id), idm(id), 198, &
                      procs_grid, status, ierr)
    ! y backward communication
    CALL MPI_SENDRECV(intvel(id)%values( &
                      intvel(id)%b(1, 1), &
                      intvel(id)%b(id, 1), &
                      intvel(id)%b_bo(3, 1)), &
                      1, MPI_slab_p(ic, id), idm(id), 199, &
                      intvel(id)%values( &
                      intvel(id)%b(1, 1), &
                      intvel(id)%b(id, 2)+1, &
                      intvel(id)%b_bo(3, 1)), &
                      1, MPI_slab_p(ic, id), idp(id), 199, &
                      procs_grid, status, ierr)


  END SUBROUTINE conv_exchange



  SUBROUTINE conv_calc

    ! Subroutine that derives previously calculated squared velocities and stores
    ! the results in the convvel block
    USE SPIKE
    USE essentials
    USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines

    IMPLICIT NONE

    INTEGER :: ic, id, ider = 1, istag
    INTEGER :: j, k
    REAL, DIMENSION(:), ALLOCATABLE :: q, psi  ! due to ANOMALIA a solution vector has to be defined
    REAL               :: r = 0.0, NaN  ! r is a dummy real used to define a NaN of type real
    NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

    !!!!!!!!!! Velocity interpolations !!!!!!!!!!
    CALL conv_interp

    !!!!!!!!!! Velocity multiplication !!!!!!!!!!
    DO ic = 1, ndims
      velvel(2, ic)%values = velvel(2, ic)%values*intvel(Ic)%values
    END DO

    !!!!!!!!!! Debugging inizialization !!!!!!!!!!
    DO ic = 1, ndims
      convvel(ic)%values = NaN
    END DO

    !!!!!!!!!! Squared velocity derivation !!!!!!!!!!

    ! x derivatives
    id = 1
    ! d(u^2)/dx
    istag = 1
    ALLOCATE(q(velvel(1, 1)%b_bo(id, 1):velvel(1, 1)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(1)%b(2, 1), convvel(1)%b(2, 2)
      DO k = convvel(1)%b(3, 1), convvel(1)%b(3, 2)

        q = velvel(1, 1)%values(:, j, k)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(1)%values(:, j, k) = psi(convvel(1)%b(id, 1):convvel(1)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)
    ! d(uv)/dx
    istag = 2
    ALLOCATE(q(velvel(2, 1)%b_bo(id, 1):velvel(2, 1)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(2)%b(2, 1), convvel(2)%b(2, 2)
      DO k = convvel(2)%b(3, 1), convvel(2)%b(3, 2)

        q = velvel(2, 1)%values(:, j, k)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(2)%values(:, j, k) = psi(convvel(2)%b(id, 1):convvel(2)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)
    ! d(wu)/dx
    istag = 2
    ALLOCATE(q(velvel(2, 3)%b_bo(id, 1):velvel(2, 3)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(3)%b(2, 1), convvel(3)%b(2, 2)
      DO k = convvel(3)%b(3, 1), convvel(3)%b(3, 2)

        q = velvel(2, 3)%values(:, j, k)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(3)%values(:, j, k) = psi(convvel(3)%b(id, 1):convvel(3)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)

    ! y derivatives
    id = 2
    ! d(uv)/dy
    istag = 2
    ALLOCATE(q(velvel(2, 1)%b_bo(id, 1):velvel(2, 1)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(1)%b(1, 1), convvel(1)%b(1, 2)
      DO k = convvel(1)%b(3, 1), convvel(1)%b(3, 2)

        q = velvel(2, 1)%values(j, :, k)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(1)%values(j, :, k) = convvel(1)%values(j, :, k) + &
                                     psi(convvel(1)%b(id, 1):convvel(1)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)
    ! d(v^2)/dy
    istag = 1
    ALLOCATE(q(velvel(1, 2)%b_bo(id, 1):velvel(1, 2)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(2)%b(1, 1), convvel(2)%b(1, 2)
      DO k = convvel(2)%b(3, 1), convvel(2)%b(3, 2)

        q = velvel(1, 2)%values(j, :, k)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(2)%values(j, :, k) = convvel(2)%values(j, :, k) + &
                                     psi(convvel(2)%b(id, 1):convvel(2)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)
    ! d(vw)/dy
    istag = 2
    ALLOCATE(q(velvel(2, 2)%b_bo(id, 1):velvel(2, 2)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(3)%b(1, 1), convvel(3)%b(1, 2)
      DO k = convvel(3)%b(3, 1), convvel(3)%b(3, 2)

        q = velvel(2, 2)%values(j, :, k)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(3)%values(j, :, k) = convvel(3)%values(j, :, k) + &
                                     psi(convvel(3)%b(id, 1):convvel(3)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)

    ! z derivatives
    id = 3
    ! d(wu)/dz
    istag = 2
    ALLOCATE(q(velvel(2, 3)%b_bo(id, 1):velvel(2, 3)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(1)%b(1, 1), convvel(1)%b(1, 2)
      DO k = convvel(1)%b(2, 1), convvel(1)%b(2, 2)

        q = velvel(2, 3)%values(j, k, :)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(1)%values(j, k, :) = convvel(1)%values(j, k, :)  + &
                                     psi(convvel(1)%b(id, 1):convvel(1)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)
    ! d(vw)/dz
    istag = 2
    ALLOCATE(q(velvel(2, 2)%b_bo(id, 1):velvel(2, 2)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(2)%b(1, 1), convvel(2)%b(1, 2)
      DO k = convvel(2)%b(2, 1), convvel(2)%b(2, 2)

        q = velvel(2, 2)%values(j, k, :)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(2)%values(j, k, :) = convvel(2)%values(j, k, :)  + &
                                     psi(convvel(2)%b(id, 1):convvel(2)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)
    ! d(w^2)/dz
    istag = 1
    ALLOCATE(q(velvel(1, 3)%b_bo(id, 1):velvel(1, 3)%b_bo(id, 2)))
    ALLOCATE(psi(cmp(id, ider, istag)%A%lb(2):cmp(id, ider, istag)%A%ub(2)))
    DO j = convvel(3)%b(1, 1), convvel(3)%b(1, 2)
      DO k = convvel(3)%b(2, 1), convvel(3)%b(2, 2)

        q = velvel(1, 3)%values(j, k, :)
        CALL SPIKE_solve(id, ider, istag, psi, q)
        ! meglio se poi la fai diventare una somma, piuttosto che un'assegnazione
        convvel(3)%values(j, k, :) = convvel(3)%values(j, k, :)  + &
                                     psi(convvel(3)%b(id, 1):convvel(3)%b(id, 2))

      END DO
    END DO
    DEALLOCATE(q)
    DEALLOCATE(psi)



  END SUBROUTINE conv_calc



END MODULE convective_term
