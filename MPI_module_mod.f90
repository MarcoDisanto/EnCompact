MODULE MPI_module
! This module contains variables identifying the MPI processes in terms of rank,
! number of processes, cartesian coordinates inside a cartesian communicator,
! neighboring processes and so on. It furthermore contains the cartesian
! communicator itself, as well as informations about its dimensions, periodicity
! and so on.
! This module contains the following procedures
!    SUBROUTINE set_topology
!    SUBROUTINE distribute_cells

    USE MPI

    IMPLICIT NONE

    INTEGER                                :: nprocs     ! number of processes
    INTEGER                                :: myid       ! index  of the current process
    INTEGER, DIMENSION(:),     ALLOCATABLE :: idm, idp   ! rank of minus (Left, Upper, Forth) and plus (Right, Lower, Back) processes
    INTEGER                                :: procs_grid ! handle to the cartesian communicator
    INTEGER                                :: ndims      ! number of dimensions of the cartesian grid (it's the size of dims and periodic)
    INTEGER, DIMENSION(:),     ALLOCATABLE :: mycoords   ! cartesian coordinates of the current process
    LOGICAL                                :: reorder    ! no reordering (see the documentation)
    INTEGER, DIMENSION(:),     ALLOCATABLE :: dims       ! array of the ndims numbers of process per direction
    LOGICAL, DIMENSION(:),     ALLOCATABLE :: periodic   ! array of periodicity along each dimension
    INTEGER, DIMENSION(:),     ALLOCATABLE :: Ntot       ! total number of cells along each dimension
    INTEGER, DIMENSION(:),     ALLOCATABLE :: N          ! number of cells to the current process along each dimension
    INTEGER, DIMENSION(:),     ALLOCATABLE :: rem        ! remainder of the division Ntot(i)/dims(i)
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Neq        ! number of equations to the current process along each dimension, for each order of derivation, and for each staggering

    INTEGER, DIMENSION(:), ALLOCATABLE :: pencil_comm ! pencil communicators
    INTEGER, DIMENSION(:), ALLOCATABLE :: slab_comm   ! slab communicators (BEWARE, the index defines the normal to the slab)


CONTAINS

    SUBROUTINE set_topology
    ! This procedure determines the number of processes, the rank of the current one,
    ! then creates a new communicator with a cartesian topology attached on it.
    ! Based on this cartesian communicator, the cartesian coordinates of the
    ! current process are determined. Furthermore, the ranks of neighboring
    ! processes are determined (neighborhood = left right up down front back).

        USE essentials, ONLY : log2int => logical2integer

        IMPLICIT NONE

        INTEGER :: ierr, errorcode, i!, id, ip
        ! INTEGER, DIMENSION(ndims) :: temp

        ! determine the total number of processes
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

        ! determine the rank (i.e. the id) of the current process
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid,  ierr)

        ! determine the distribution of processes in a Cartesian grid
        CALL MPI_DIMS_CREATE(nprocs, ndims, dims, ierr)

        ! create the cartesian communicator (communicator with cartesian topology)
        CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periodic, reorder, procs_grid, ierr)

        ! determine the cartesian coordinates of the current process
        ALLOCATE(mycoords(ndims))
        CALL MPI_CART_COORDS(procs_grid, myid, ndims, mycoords, ierr)

        ! determine the ranks of neighboring processes along the dimensions
        ALLOCATE(idm(ndims), idp(ndims))
        DO i = 1, ndims
            ! NB: the second argument is the direction along which the shift
            ! is going to be performed, and it is indexed by 1 (like in C,
            ! not Fortran).
            CALL MPI_CART_SHIFT(procs_grid, i-1, 1, idm(i), idp(i), ierr)
        END DO

        !ALLOCATE(touchm(ndims), touchp(ndims), touch(ndims))
        !touchm = -log2int(idm == MPI_PROC_NULL)
        !touchp = +log2int(idp == MPI_PROC_NULL)
        !touch  = touchp - touchm
        ! touchm(i) == -1 (resp. touchp(i) == +1) if the current process touches
        ! the minus (resp. plus) boundary along the i-th direction, which can be
        ! L, D, or B (resp. R, U, or F), == 0 if it doesn't.
        ! touch(i) == 1 if the current process touches either minus or plus
        ! boundary along the i-th direction, == 2 if it touches both boundaries,
        ! == 0 if it doesn't touch any of them.

        ! create ndims "pencil" communicators for each process, pencil_comm(i)
        ! containing the process itself and the other processes whose coordinates
        ! only differs in the i-th direction
        ALLOCATE(pencil_comm(ndims))
        ALLOCATE(slab_comm(ndims))
        create_ndims_pencil_communicators: SELECT CASE (ndims)
        CASE (1)
            PRINT *, 'ERROR in MPI_module_mod.f90: per ndims == 1 non credo di aver fatto un check.'
            PRINT *, 'Inoltre una Navier-Stokes 1D (non quasi-1D) incompressile ha poco senso...'
            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            STOP

        CASE (2)
            ! pencils
            CALL MPI_CART_SUB(procs_grid, [.true., .false.], pencil_comm(1), ierr)
            CALL MPI_CART_SUB(procs_grid, [.false., .true.], pencil_comm(2), ierr)
            ! slabs
            CALL MPI_CART_SUB(procs_grid, [.false., .true.], slab_comm(1), ierr)
            CALL MPI_CART_SUB(procs_grid, [.true., .false.], slab_comm(2), ierr)

        CASE (3)
            ! pencils
            CALL MPI_CART_SUB(procs_grid, [.true., .false., .false.], pencil_comm(1), ierr)
            CALL MPI_CART_SUB(procs_grid, [.false., .true., .false.], pencil_comm(2), ierr)
            CALL MPI_CART_SUB(procs_grid, [.false., .false., .true.], pencil_comm(3), ierr)
            ! slabs
            CALL MPI_CART_SUB(procs_grid, [.false., .true., .true.], slab_comm(1), ierr)
            CALL MPI_CART_SUB(procs_grid, [.true., .false., .true.], slab_comm(2), ierr)
            CALL MPI_CART_SUB(procs_grid, [.true., .true., .false.], slab_comm(3), ierr)

        CASE (4:)
            PRINT *, "ERROR in MPI_module_mod.f90: 'Cube 2: Hypercube' is a great movie!"
            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            STOP

        CASE DEFAULT
            PRINT *, "ERROR in MPI_module_mod.f90: non-positive ndims."
            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            STOP

        END SELECT create_ndims_pencil_communicators

    END SUBROUTINE set_topology


    SUBROUTINE distribute_cells
    ! This procedure determines how many cells are assigned to the current
    ! process along each dimension (also in the case that nprocs doesn't divide
    ! Ntot evenly).

        USE essentials, ONLY : log2int => logical2integer

        IMPLICIT NONE

        INTEGER :: i

        ! Allocation
        ! ALLOCATE(N(ndims), rem(ndims), Nft(ndims), Nfi(ndims), Neq(ndims, 0:2, 2))
        ALLOCATE(N(ndims), rem(ndims), Neq(ndims, 0:2, 2))

        DO i = 1, ndims                    ! For each direction...

            ! determine the number of cells handled by the current process:
            N(i)  = Ntot(i)/dims(i)        ! distribute the total number of cells, ...
            rem(i) = mod(Ntot(i), dims(i)) ! ... consider the remainder, ...
            IF (mycoords(i) < rem(i)) THEN ! ... and give 1 cell more...
                N(i) = N(i) + 1            ! ... to the leading processes
            END IF

            ! set the number of equation for each direction, order of drivation, and staggering
            Neq(i, 0, 1) = N(i) - log2int(idm(i) == MPI_PROC_NULL)
            Neq(i, 0, 2) = N(i)
            Neq(i, 1, 1) = N(i) + log2int(idp(i) == MPI_PROC_NULL) ! incognite anche le facce di bordo, poi scartate (questo è ANOMALO [cerca questa parola nei file per individuare le dipendenze])
            Neq(i, 1, 2) = N(i)
            Neq(i, 2, 1) = N(i)
            Neq(i, 2, 2) = N(i) - log2int(idm(i) == MPI_PROC_NULL)

        END DO

        ! Now rem can be used for other purposes, and so it will be... (it will
        ! contain the remainder of the division of the number of unknowns of the
        ! various schemes by the number of processes along a direction..)

    END SUBROUTINE distribute_cells

END MODULE MPI_module

            ! originariamente in distribure_cells
            ! determine the number of internal faces handled by the current process:
            ! the same, in number, as the cells ...
            !Nfi(i) = N(i)
            !! ... except for the first process in the case of non-periodic direction ...
            !! print *, (((mycoords(i) == 0) .AND. (periodic(i) .EQV. .FALSE.)) .eqv. (idm(i) == MPI_PROC_NULL))
            !IF (idm(i) == MPI_PROC_NULL) THEN
            !! ... that handles one less face (the first phisical one)
            !    Nfi(i) = Nfi(i) - 1
            !END IF

            ! determine the number of total faces handled by the current process:
            ! the same, in number, as the cells ...
            !Nft(i) = N(i)
            !! ... except for the last process in the case of non-periodic direction ...
            !! print *, (((mycoords(i) == dims(i) - 1) .AND. (periodic(i) .EQV. .FALSE.)) .eqv. (idp(i) == MPI_PROC_NULL))
            !IF (idp(i) == MPI_PROC_NULL) THEN
            !! ... that handles one more face (the last phisical one)
            !    Nft(i) = Nft(i) + 1
            !END IF
