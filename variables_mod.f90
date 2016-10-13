MODULE variables
! This module contains variables and functions related to the dependent
! physical variables, i.e. physical quantities, such as velocity, pressure, and
! maybe others.
!
! Procedures:
!    SUBROUTINE set_variables
!    SUBROUTINE set_nans

    IMPLICIT NONE
    
    ! questo deve essere spostato altrove (nell'input) TO DO
    ! in realtà dovrebbe essere collegato anche a ndims;
    ! qualcosa tipo ndims + 1 + altro_scalare + ndims*altro_vettoriale
    ! (ndims componenti di velocità, 1 pressione, altre grandezze)
    INTEGER, PARAMETER :: n_flow_variables = 4
    
    ! WORK IN PROGRESS
    ! TODO: il tipo pencil è proprio il tipo grid1D. Bisogna sostituire
    ! quest'ultimo.
    TYPE pencil
        REAL, DIMENSION(:), ALLOCATABLE :: values
    END TYPE pencil

    TYPE flat
        REAL, DIMENSION(:,:), POINTER :: values => NULL()
    END TYPE flat
    ! WORK IN PROGRESS

    ! Type used to contain the fluiddynamic variables.
    TYPE block
        REAL,    DIMENSION(:,:,:), ALLOCATABLE :: values
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: b    ! bounds not considering  borders  nor overlap
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: b_bc ! bounds considering      borders, not overlap
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: b_ol ! bounds not considering  borders, but overlap
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: b_bo ! bounds considering both borders  and overlap
    END TYPE block
    TYPE(block), DIMENSION(n_flow_variables), TARGET :: uvwp ! TO DO: meglio 0:n_flow_variables-1, associando 0 alla pressione e 1,... alle componenti di velocità

    ! MPI types to manage communications
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: MPI_flats_m, MPI_flats_p
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: rdispls_m, rdispls_p, sdispls_m, sdispls_p

CONTAINS

    SUBROUTINE set_variables
    ! This subroutine allocates space for variables used to contain fluid
    ! dynamic quantities pertaining to the current process.
    ! These multidimensional arrays are allocated with such dimensions that
    ! they contain values relative to overlapping points (both cells and faces)
    ! with neighboring processes.

        USE MPI,        ONLY : MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_ORDER_FORTRAN, MPI_PROC_NULL, MPI_ADDRESS_KIND
        USE MPI_module, ONLY : ndims, N, idp, idm, myid
        USE essentials, ONLY : log2int => logical2integer, KronDelta, printmatrix
        USE compact,    ONLY : cmp

        IMPLICIT NONE

        ! TODO: check for unused variables (a debugging compiler flag can be used)
        INTEGER :: errorcode, ierr, ic, id, ld1, ud1, ld2, ud2, i, id2
        INTEGER, DIMENSION(ndims) :: array_of_sizes, array_of_subsizes_m, array_of_subsizes_p, array_of_starts
        INTEGER :: dpsize!, sinb
        INTEGER :: address000, address1
        INTEGER, DIMENSION(ndims) :: s_start_coords_m, s_start_coords_p, r_start_coords_m, r_start_coords_p

        ! TODO: MPI_flats_m and MPI_flats_p can be declared and allocated in the MPI_module_mod.f90 and imported
        ! here from that module
        ALLOCATE(MPI_flats_m(ndims,ndims), MPI_flats_p(ndims,ndims))
        ALLOCATE(sdispls_m(ndims,ndims), sdispls_p(ndims,ndims), rdispls_m(ndims,ndims), rdispls_p(ndims,ndims))

        DO i = 1, n_flow_variables
            ALLOCATE(uvwp(i)%b(ndims, 2), uvwp(i)%b_bc(ndims, 2), uvwp(i)%b_ol(ndims, 2), uvwp(i)%b_bo(ndims, 2))
        END DO

        select_the_number_of_dimensions: SELECT CASE (ndims)

            CASE (1)

                PRINT *, 'ERROR in compact_mod.f90: una Navier-Stokes 1D (non quasi-1D) incompressile ha poco senso.&
                        & Esecuzione interrotta.'
                CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                STOP


            CASE (2)

                PRINT *, 'ERROR in compact_mod.f90: ancora non ho codificato il caso 2D. Esecuzione interrotta.'
                CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                STOP


            CASE (3)

                for_each_component: DO ic = 1, ndims

                    for_each_direction_compute_size: DO id = 1, ndims

                        ! bounds of the velocity array not considering boundaries nor overlap
                        uvwp(ic)%b(id,:) = [1 + KronDelta(id,ic)*log2int(idm(id) == MPI_PROC_NULL), N(id)]

                        ! retrieve the number of diagonals of the B matrix for
                        ! interpolation and second derivative along the current
                        ! direction for the current component (take the max between the two)
                        ld1 = MAX(cmp(id,0,1)%B%lid, cmp(id,2,1)%B%lid) ! number of lower diagonals from centers (to faces or to centers)
                        ud1 = MAX(cmp(id,0,1)%B%uid, cmp(id,2,1)%B%uid) ! number of upper diagonals from centers (to faces or to centers)
                        ld2 = MAX(cmp(id,0,2)%B%lid, cmp(id,2,2)%B%lid) ! number of lower diagonals from faces   (to faces or to centers)
                        ud2 = MAX(cmp(id,0,2)%B%uid, cmp(id,2,2)%B%uid) ! number of upper diagonals from faces   (to faces or to centers)

                        ! TODO: the following (till END DO) can be simplified a lot (the modifications
                        ! to obtain b_bc and b_ol are mutually exclusive!)
                        ! modify the lower bound of velocity array considering the boundaries
                        uvwp(ic)%b_bc(id,1) = uvwp(ic)%b(id,1) &
                                  & - log2int(idm(id) == MPI_PROC_NULL) ! minus 1 if touches the minus border along the current direction
                        ! modify the lower bound of velocity array considering the overlap (regardless of the possible non periodicity)
                        uvwp(ic)%b_ol(id,1) = uvwp(ic)%b(id,1) &
                                  & - ld1*(1 - KronDelta(id,ic)) &! minus the number of lower diagonals from centers for the non-ic-th direction 
                                  & - ld2*KronDelta(id,ic)        ! minus the number of lower diagonals from faces   for the     ic-th direction
                        ! modify the lower bound of velocity array considering the boundaries and the overlap
                        uvwp(ic)%b_bo(id,1) = uvwp(ic)%b(id,1) &
                                  & - log2int(idm(id) == MPI_PROC_NULL) & ! minus 1 if touches the minus border along the current direction
                                  & - log2int(idm(id) /= MPI_PROC_NULL)*ld1*(1 - KronDelta(id,ic)) &! minus the number of lower diagonals from centers for the non-ic-th direction 
                                  & - log2int(idm(id) /= MPI_PROC_NULL)*ld2*KronDelta(id,ic)        ! minus the number of lower diagonals from faces   for the     ic-th direction

                        ! modify the lower bound of velocity array considering the boundaries
                        uvwp(ic)%b_bc(id,2) = uvwp(ic)%b(id,2) &
                                  & + log2int(idp(id) == MPI_PROC_NULL) ! plus 1 if touches the plus border along the current direction
                        ! modify the lower bound of velocity array considering the overlap (regardless of the possible non periodicity)
                        uvwp(ic)%b_ol(id,2) = uvwp(ic)%b(id,2) &
                                  & + ud1*(1 - KronDelta(id,ic)) &! plus the number of upper diagonals from centers for the non-ic-th direction
                                  & + ud2*KronDelta(id,ic)        ! plus the number of upper diagonals from faces   for the     ic-th direction
                        ! modify the lower bound of velocity array considering the boundaries and the overlap
                        uvwp(ic)%b_bo(id,2) = uvwp(ic)%b(id,2) &
                                  & + log2int(idp(id) == MPI_PROC_NULL) & ! plus 1 if touches the plus border along the current direction
                                  & + log2int(idp(id) /= MPI_PROC_NULL)*ud1*(1 - KronDelta(id,ic)) &! plus the number of upper diagonals from centers for the non-ic-th direction
                                  & + log2int(idp(id) /= MPI_PROC_NULL)*ud2*KronDelta(id,ic)        ! plus the number of upper diagonals from faces   for the     ic-th direction

                    END DO for_each_direction_compute_size

                    ! allocate velocity component array
                    ALLOCATE(uvwp(ic)%values(uvwp(ic)%b_bo(1,1):uvwp(ic)%b_bo(1,2), &
                                           & uvwp(ic)%b_bo(2,1):uvwp(ic)%b_bo(2,2), &
                                           & uvwp(ic)%b_bo(3,1):uvwp(ic)%b_bo(3,2)))

                    ! shape of the velocity component array with ghost cells and boundaries
                    array_of_sizes = SHAPE(uvwp(ic)%values)

                    for_each_direction_compute_subsize: DO id = 1, ndims

                        ! shape of the velocity component sub-array
                        ! TODO: in a far future I could combine array_of_subsizes_m and array_of_subsizes_p
                        ! into one 2D array so that the following DO would contain just one assignation
                        ! in another DO = 1,2
                        DO id2 = 1, ndims

                            ! The size along the id2-th direction of the ndims-ranked subarray relative to the
                            ! exchange of the ic-th component along the id-th direction is equal to the
                            ! number of overlap layers if id2 == id, whereas, if id2 /= id, it is N(id2) diminished by 1
                            ! if the the process touches the minus boundary along the ic-th direction.
                            array_of_subsizes_m(id2) = KronDelta(id,id2)*(uvwp(ic)%b(id2,1) - uvwp(ic)%b_ol(id2,1)) + &
                                                & (1 - KronDelta(id,id2))*(N(id2) - &
                                                & log2int(idm(id2) == MPI_PROC_NULL)*KronDelta(id2,ic))

                            array_of_subsizes_p(id2) = KronDelta(id,id2)*(uvwp(ic)%b_ol(id2,2) - uvwp(ic)%b(id2,2)) + &
                                                & (1 - KronDelta(id,id2))*(N(id2) - &
                                                & log2int(idm(id2) == MPI_PROC_NULL)*KronDelta(id2,ic))

                            ! TODO: a causa della non periodicità, la dei displacement sembra essere un fatto
                            ! davvero tosto da gestire; a proposito del seguente controllo, sono abbstanza
                            ! sicuro del blocco ELSE, un po' meno del blocco IF ... THEN, soprattutto relativamente
                            ! alla definizione dei displacement dei send buffer (i.e. s_start_...)
                            IF (id2 == id) THEN
                                r_start_coords_m(id2) = uvwp(ic)%b_bo(id,1)
                                r_start_coords_p(id2) = uvwp(ic)%b(id,2) + 1
                                ! s_start_coords_m(id2) = r_start_coords_m(id2) + N(id2)
                                ! s_start_coords_p(id2) = r_start_coords_p(id2) - N(id2)
                                s_start_coords_m(id2) = N(id2) - (uvwp(ic)%b(id,1) - uvwp(ic)%b_ol(id,1)) + 1
                                s_start_coords_p(id2) = 1! + log2int(idm(ic) == MPI_PROC_NULL)
                            ELSE
                                ! the point from witch the 4 flats (two to be received and two to be sent) start
                                ! has the same coordinates for the 4 flats along the non-id-th direction
                                s_start_coords_m(id2) = uvwp(ic)%b(id,1)
                                s_start_coords_p(id2) = uvwp(ic)%b(id,1)
                                r_start_coords_m(id2) = uvwp(ic)%b(id,1)
                                r_start_coords_p(id2) = uvwp(ic)%b(id,1)
                            END IF

                        END DO

                        !if (myid == 5 .and. ic == 1) then
                        !    print *, 'ic', ic, 'id', id, 'array_of_sizes', array_of_sizes, &
                        !                               & 'array_of_subsizes_m', array_of_subsizes_m, &
                        !                               & 'array_of_subsizes_p', array_of_subsizes_p
                        !endif

                        array_of_starts = [0, 0, 0] ! We can't use negative starts in MPI calls (this is C legacy),
                                                    ! so we set them to zero and rely upon displacements

                        ! Create MPI_flats_m(ic,id) type, which represents the arrangement in memory of the values
                        ! of the ic-th component of velocity to be sent to (received from) the process on the plus
                        ! (minus) side along the id-th direction.
                        CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                                    & array_of_sizes, array_of_subsizes_m, array_of_starts, &
                                                    & MPI_ORDER_FORTRAN, &
                                                    & MPI_DOUBLE_PRECISION, MPI_flats_m(ic,id), &
                                                    & ierr)

                        ! Create MPI_flats_p(ic,id) type, which represents the arrangement in memory of the values
                        ! of the ic-th component of velocity to be received from (sent to) the process on the plus
                        ! (minus) side along the id-th direction.
                        CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                                                    & array_of_sizes, array_of_subsizes_p, array_of_starts, &
                                                    & MPI_ORDER_FORTRAN, &
                                                    & MPI_DOUBLE_PRECISION, MPI_flats_p(ic,id), &
                                                    & ierr)

                        ! TODO: remember to deallocate the stuff used only to build this types and not elsewhere
                        CALL MPI_TYPE_COMMIT(MPI_flats_m(ic,id),ierr)
                        CALL MPI_TYPE_COMMIT(MPI_flats_p(ic,id),ierr)

                        ! CALL MPI_TYPE_SIZE(MPI_flats_m(ic,id),sinb,ierr)
                        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dpsize, ierr)
                        !IF (myid == 5) THEN
                        !    print *, myid, ic, id, array_of_subsizes_m, array_of_subsizes_p!, sinb/dpsize
                        !ENDIF
                        !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
                        !IF (myid == 2) THEN
                        !    print *, myid, ic, id, array_of_subsizes_m, array_of_subsizes_p!, sinb/dpsize

                        !ENDIF
                        !print *, myid, ic, id, array_of_subsizes_m, array_of_subsizes_p

                        ! address in memory of the reference point (smallest coordinates along the 3 dimensions, uvwp(.)%b_bo(.,1))
                        CALL MPI_ADDRESS(uvwp(ic)%values(uvwp(ic)%b_bo(1,1), uvwp(ic)%b_bo(2,1), uvwp(ic)%b_bo(3,1)), &
                                                       & address000,ierr)
                        ! address in memory of of starting points of the two receiving and two sending flats
                        CALL MPI_ADDRESS(uvwp(ic)%values(r_start_coords_m(1), r_start_coords_m(2), r_start_coords_m(3)), &
                                                       & rdispls_m(ic,id), ierr)
                        CALL MPI_ADDRESS(uvwp(ic)%values(r_start_coords_p(1), r_start_coords_p(2), r_start_coords_p(3)), &
                                                       & rdispls_p(ic,id), ierr)
                        CALL MPI_ADDRESS(uvwp(ic)%values(s_start_coords_m(1), s_start_coords_m(2), s_start_coords_m(3)), &
                                                       & sdispls_m(ic,id), ierr)
                        CALL MPI_ADDRESS(uvwp(ic)%values(s_start_coords_p(1), s_start_coords_p(2), s_start_coords_p(3)), &
                                                       & sdispls_p(ic,id), ierr)
                        ! displacements with respect to the reference address000
                        rdispls_m(ic,id) = rdispls_m(ic,id) - address000
                        rdispls_p(ic,id) = rdispls_p(ic,id) - address000
                        sdispls_m(ic,id) = sdispls_m(ic,id) - address000
                        sdispls_p(ic,id) = sdispls_p(ic,id) - address000

                        IF (myid == 5 .and. ic == 2) THEN
                            ! print *, uvwp(2)%b_bo(1,1),uvwp(2)%b_bo(2,1),uvwp(2)%b_bo(3,1)
                            ! print *, uvwp(2)%b(1,1),uvwp(2)%b(2,2)+1,uvwp(2)%b(3,1)
                            ! print *, id, rdispls_m/dpsize, rdispls_p/dpsize
                            ! print *, id, sdispls_m/dpsize, sdispls_p/dpsize
                        END IF

                    END DO for_each_direction_compute_subsize

                END DO for_each_component

                ! (Here the handles  MPI_flats_m and MPI_flats_p are set)

                !CALL MPI_ADDRESS(uvwp(2)%values(uvwp(2)%b_bo(1,1),uvwp(2)%b_bo(2,1),uvwp(2)%b_bo(3,1)),address000,ierr)
                !CALL MPI_ADDRESS(uvwp(2)%values(uvwp(2)%b(1,1),uvwp(2)%b(2,2)+1,   uvwp(2)%b(3,1)),   address1,  ierr)
                !CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dpsize, ierr)
                !IF (myid == 5) THEN
                !    print *, uvwp(2)%b_bo(1,1),uvwp(2)%b_bo(2,1),uvwp(2)%b_bo(3,1)
                !    print *, uvwp(2)%b(1,1),uvwp(2)%b(2,2)+1,uvwp(2)%b(3,1)
                !    print *, (address1 - address000)/dpsize
                !END IF

        END SELECT select_the_number_of_dimensions

    END SUBROUTINE set_variables


    SUBROUTINE set_nans
    ! This subroutine insert NaNs in the grids in positions that should never be
    ! used. In particular it inserts quiet NaNs which propagate silently.

        USE MPI,        ONLY : MPI_PROC_NULL
        USE MPI_module, ONLY : ndims, idm, idp, N
        USE essentials, ONLY : KronDelta
        USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines
        ! such as IEEE_QUIET_NAN, IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF

        IMPLICIT NONE

        INTEGER :: ic

        ! IEEE
        REAL :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
        NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

        for_each_component: DO ic = 1, ndims

            if_touches_LD_edge: IF (idm(1) == MPI_PROC_NULL .and. idm(2) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(KronDelta(ic,1),KronDelta(ic,2),:) = NaN
            END IF if_touches_LD_edge

            if_touches_LB_edge: IF (idm(1) == MPI_PROC_NULL .and. idm(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(KronDelta(ic,1),:,KronDelta(ic,3)) = NaN
            END IF if_touches_LB_edge

            if_touches_DB_edge: IF (idm(2) == MPI_PROC_NULL .and. idm(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(:,KronDelta(ic,2),KronDelta(ic,3)) = NaN
            END IF if_touches_DB_edge

            if_touches_LU_edge: IF (idm(1) == MPI_PROC_NULL .and. idp(2) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(KronDelta(ic,1),N(2)+1,:) = NaN
            END IF if_touches_LU_edge

            if_touches_LF_edge: IF (idm(1) == MPI_PROC_NULL .and. idp(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(KronDelta(ic,1),:,N(3)+1) = NaN
            END IF if_touches_LF_edge

            if_touches_RD_edge: IF (idp(1) == MPI_PROC_NULL .and. idm(2) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(N(1)+1,KronDelta(ic,2),:) = NaN
            END IF if_touches_RD_edge

            if_touches_RB_edge: IF (idp(1) == MPI_PROC_NULL .and. idm(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(N(1)+1,:,KronDelta(ic,3)) = NaN
            END IF if_touches_RB_edge

            if_touches_DF_edge: IF (idm(2) == MPI_PROC_NULL .and. idp(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(:,KronDelta(ic,2),N(3)+1) = NaN
            END IF if_touches_DF_edge

            if_touches_UB_edge: IF (idp(2) == MPI_PROC_NULL .and. idm(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(:,N(2)+1,KronDelta(ic,3)) = NaN
            END IF if_touches_UB_edge

            if_touches_RU_edge: IF (idp(1) == MPI_PROC_NULL .and. idp(2) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(N(1)+1,N(2)+1,:) = NaN
            END IF if_touches_RU_edge

            if_touches_RF_edge: IF (idp(1) == MPI_PROC_NULL .and. idp(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(N(1)+1,:,N(3)+1) = NaN
            END IF if_touches_RF_edge

            if_touches_UF_edge: IF (idp(2) == MPI_PROC_NULL .and. idp(3) == MPI_PROC_NULL) THEN
                uvwp(ic)%values(:,N(2)+1,N(3)+1) = NaN
            END IF if_touches_UF_edge

        END DO for_each_component

    END SUBROUTINE set_nans

END MODULE variables
