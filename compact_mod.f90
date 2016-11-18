MODULE compact
! This module contains the following procedures
!    SUBROUTINE set_compacts
!    SUBROUTINE compact_coeffs(LHS_s,RHS_s,LHS_o,RHS_o,c,coeffs)
!    SUBROUTINE calc_compact_matrices(xL, xR, Adummy, Bdummy, schemes)
! where the first calls the second, and the second calls the third

    USE bandedmatrix, ONLY : CDS

    IMPLICIT NONE

    TYPE compact_type
        TYPE(CDS)                                        :: A
        TYPE(CDS)                                        :: B
        CHARACTER(len = 20), DIMENSION(:),   ALLOCATABLE :: sch ! character string identifying the scheme (central and boundary)
        INTEGER                                          :: N ! number of equations (i.e. number of rows)
    END TYPE compact_type

    TYPE(compact_type), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: cmp
    ! The three indices are supposed to refer to:
    !  -- direction (ndims, that is 1, 2, or 3)
    !  -- order of derivation (0, 1, or 2)
    !  -- relative staggering of LHS and RHS (1: c2c or c2f, or 2: f2c or f2f)

    TYPE compact_schemes
        CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE :: sch
    END TYPE compact_schemes
    TYPE(compact_schemes), DIMENSION(:,:,:), ALLOCATABLE :: all_schemes ! structure containing all the schemes used

    ! global temporary matrices and schemes
    TYPE(CDS) :: Ag
    TYPE(CDS) :: Bg
    CHARACTER(len = 20), DIMENSION(:), ALLOCATABLE :: sch

    ! TYPE(CDS), POINTER :: A ! per usare questi anziché Ag e Bg,
    ! TYPE(CDS), POINTER :: B ! bisogna mettere l'attributo target a cmp

CONTAINS

    SUBROUTINE set_compacts
    ! This procedure sets the matrices A and B of the various compact schemes
    ! (first and second derivatives schemes, as well as interpolations schemes,
    ! and so on...).

        USE MPI_module
        USE grids
        USE library
        USE essentials

        IMPLICIT NONE

        ! MPI-concerning variables
        INTEGER :: ierr, errorcode

        INTEGER :: i, j, k ! generic indices

        ! *********************************************************************
        ! The following internal variables are replaced time by time depending on the scheme considered
        ! NB: this can be done since the couples of compact matrices A and B are
        ! computed one by one (sequentially) and scattered, so that a lot of
        ! master memory is saved.
        INTEGER :: il, iu   ! index of first and last equations of the generic process (il can be different from 1 and similarly iu ...)
        INTEGER                            :: pindex ! generic process' index
        INTEGER, DIMENSION(:), ALLOCATABLE :: pcoord ! generic process' coordinates in cartesian communicator (not those of the current process, which are stored in the array mycoords of the MPI_module module)
        INTEGER, DIMENSION(:), ALLOCATABLE :: nequ_v ! array containing the number of equations of the various processes
        INTEGER, DIMENSION(2) :: Bband_tot ! lower and upper band (of the banded matrix B) for all the schemes
        INTEGER, DIMENSION(2) :: Bband_int ! lower and upper band (of the banded matrix B) for the central scheme
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: nequ_m ! array containing the number of equations of the various processes
        ! The preceding internal variables are replaced time by time depending on the scheme considered
        ! *********************************************************************


        ! *********************************************************************
        ! The matrix B of each compact scheme, because of the boundaries, has a
        ! band that is wider near the boundaries than in the inner part of the
        ! mesh. These matrices are built as a whole by the master process and
        ! are stored in their CDS (Compressed Diagonal Storage) form. This means
        ! that the most external diagonals have just a few non-null elements
        ! (few upper elements of upper diagonals, few lower elements of lower
        ! diagonals), but these diagonals must be stored in the CDS form.
        ! When the matrix B has to be scattered from master to slaves, the
        ! problem of storing (parts of) diagonals made up of a lot of zeros only
        ! concernes the two processes that are near the two boundaries; for the
        ! internal processes, such (parts of) diagonals are identically zero,
        ! and can be removed.
        ! *********************************************************************

        IF (myid == 0) THEN
            ! allocates the number of equation of each process (a gather will be used)
            ALLOCATE(nequ_v(0:nprocs-1))

            SELECT CASE (ndims)
            CASE (1)
                ALLOCATE(nequ_m(0:dims(1)-1, 0:0, 0:0))
            CASE (2)
                ALLOCATE(nequ_m(0:dims(1)-1, 0:dims(2)-1, 0:0))
            CASE (3)
                ALLOCATE(nequ_m(0:dims(1)-1, 0:dims(2)-1, 0:dims(3)-1))
            END SELECT

            ! allocates the coords of the generic process
            ALLOCATE(pcoord(ndims))

        ELSE
            ! puppet allocation
            ALLOCATE(nequ_v(1))

        END IF

        ! allocates directions (ndims), order of derivation (0, 1, 2), and staggering of LHS and RHS (1: c2c or c2f, 2: f2c or f2f)
        ALLOCATE(cmp(ndims, 0:2, 2))

        for_each_direction: DO i = 1, ndims

            for_each_order_of_derivation: DO j = 0, 2

                for_each_staggering: DO k = 1, 2

                    ! The matrices A and B are built by the master process only
                    ! TO DO: questa parte può essere cambiata considerando che
                    ! esistono i pencil-master, il che consentirebbe uno sgravio
                    ! al master globale.
                    master_proc_builds_compacts: IF (myid == 0) THEN

                        ! ALLOCATE(sch(lbound(all_schemes(i,j,k)%sch:ubound(all_schemes(i,j,k)%sch)) ! the allocation is automatic
                        sch = all_schemes(i,j,k)%sch

                        ! lower and upper bounds of full counterparts of matrices A and B along the two dimensions:
                        ! - both matrices A and B have the rows indexed as the unknowns, that is
                        !   the lower and upper indices are those of the LHS grid, that is ( , , ,1)
                        ! - matrix A is square, so the columns are indexed as the rows
                        ! - matrix B has columns indexed as the the knowns, that is the lower and
                        !   upper indices are those of the RHS grid, that is ( , , ,2)
                        Ag%lb = [lbound(all_grids(i,j,k,1)%g), lbound(all_grids(i,j,k,1)%g)]
                        Ag%ub = [ubound(all_grids(i,j,k,1)%g), ubound(all_grids(i,j,k,1)%g)]
                        Bg%lb = [lbound(all_grids(i,j,k,1)%g), lbound(all_grids(i,j,k,2)%g)]
                        Bg%ub = [ubound(all_grids(i,j,k,1)%g), ubound(all_grids(i,j,k,2)%g)]

                        ! lower and upper band of A and B matrix...
                        Ag%ld = 1
                        Ag%ud = 1

                        Bband_tot = detBband(sch)      ! ... considering all  the      rows
                        Bband_int = detBband(sch(0:0)) ! ... considering only internal rows
                        Bg%ld = -Bband_tot(1) ! number of lower diagonals (positive if lower diagonals are present)
                        Bg%ud = +Bband_tot(2) ! number of upper diagonals (positive if upper diagonals are present)
                        ! TO DO: dopo le due precedenti assegnazioni, le
                        ! variabili Bband_tot e Bband_int non dovrebbero essere
                        ! più usate, il che significa che devo cercarle e
                        ! sostituirle (anche nel modulo set_pressure_mod.f90)

                        ! allocate sparse counterparts of matrices A and B based on number of rows and diagonals (CDS)
                        ALLOCATE(Ag%matrix(Ag%lb(1):Ag%ub(1),     -1:+1))
                        ALLOCATE(Bg%matrix(Bg%lb(1):Bg%ub(1), -Bg%ld:+Bg%ud))

                        ! fill the matrices A and B
                        CALL calc_compact_matrices(all_grids(i,j,k,1)%g, all_grids(i,j,k,2)%g, &
                          & Ag%matrix, Bg%matrix, sch)

                        ! ******************* da eliminare
                        ! per usare printmatrix devi usare il modulo essentials
                        !IF (i == 1 .and. j == 1 .and. k == 2) THEN
                        !    print *, 'Lhs'
                        !    call printmatrix(all_grids(i,j,k,1)%g)
                        !    print *, 'Rhs'
                        !    call printmatrix(all_grids(i,j,k,2)%g)
                        !    print *, 'A (compact_mod.f90)'
                        !    CALL printmatrix(Ag%matrix)
                        !    print *, 'B (compact_mod.f90)'
                        !    CALL printmatrix(Bg%matrix)
                        !END IF
                        ! ******************* da eliminare

                    END IF master_proc_builds_compacts

                    ! the master process broadcasts Bband_int to every process
                    CALL MPI_BCAST(Bband_int, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

                    ! the master process broadcasts Bband_tot to every process
                    CALL MPI_BCAST(Bband_tot, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

                    ! the master process gathers the number of equations of all the processes
                    CALL MPI_GATHER(Neq(i,j,k), 1, MPI_INTEGER, nequ_v, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

                    master_proc_scatters_compacts: IF (myid == 0) THEN ! master process sends chunks of A to other processes

                        ! the master process arranges the number of equations of the
                        ! processes in an array sized as procs_grid
                        SELECT CASE (ndims)
                        CASE (1)
                            nequ_m = RESHAPE(nequ_v, [dims(1:1), 1, 1])
                            PRINT *, 'WARNING in compact_mod.f90: per ndims == 1 non credo di aver fatto un check.'
                            PRINT *, 'Inoltre una Navier-Stokes 1D (non quasi-1D) incompressile ha poco senso...'
                            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                            STOP

                        CASE (2)
                            nequ_m = RESHAPE(nequ_v, [dims(1:2), 1], order = [3, 2, 1]) ! NB: questo 'order =' ci deve stare, perché il compilatore è un figlio di buona donna
                            PRINT *, 'WARNING in compact_mod.f90: per ndims == 2 non credo di aver fatto un check.'
                            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                            STOP

                        CASE (3)
                            nequ_m = RESHAPE(nequ_v, dims(1:3), order = [3, 2, 1]) ! NB: questo 'order =' ci deve stare, perché il compilatore è un figlio di buona donna
                        END SELECT

                        ! cycle on other processes (included the master itself)
                        send_to_each_process: DO pindex = 0, nprocs-1

                            ! determination of cartesian coordinates of the pindex process
                            CALL MPI_CART_COORDS(procs_grid, pindex, ndims, pcoord, ierr)
                            ! lower and upper indices of the chunk to be sent to the pindex process
                            SELECT CASE (i)
                            CASE (1)
                                il = SUM(nequ_m(0:pcoord(1)-1, pcoord(2), pcoord(3))) + Ag%lb(1)
                            CASE (2)
                                il = SUM(nequ_m(pcoord(1), 0:pcoord(2)-1, pcoord(3))) + Ag%lb(1)
                            CASE (3)
                                il = SUM(nequ_m(pcoord(1), pcoord(2), 0:pcoord(3)-1)) + Ag%lb(1)
                            END SELECT
                            iu = il + nequ_v(pindex) - 1

                            ! send the chunks of A
                            CALL MPI_SEND(Ag%matrix(il:iu,:), 3*nequ_v(pindex), &
                                    & MPI_DOUBLE_PRECISION, pindex, 123, MPI_COMM_WORLD, ierr)

                            ! send the chunks of B
                            IF (pcoord(i) == 0) THEN ! process at the 'minus' border

                                CALL MPI_SEND(Bg%matrix(il:iu, Bband_int(1):Bband_tot(2)), &
                                        & (Bband_tot(2) - Bband_int(1) + 1)*nequ_v(pindex), &
                                        & MPI_DOUBLE_PRECISION, pindex, 12, MPI_COMM_WORLD, ierr)

                            ELSE IF (pcoord(i) == dims(i) - 1) THEN ! process at the 'plus' border

                                CALL MPI_SEND(Bg%matrix(il:iu, Bband_tot(1):Bband_int(2)), &
                                        & (Bband_int(2) - Bband_tot(1) + 1)*nequ_v(pindex), &
                                        & MPI_DOUBLE_PRECISION, pindex, 22, MPI_COMM_WORLD, ierr)

                            ELSE ! process in the center

                                CALL MPI_SEND(Bg%matrix(il:iu, Bband_int(1):Bband_int(2)), &
                                        & (Bband_int(2) - Bband_int(1) + 1)*nequ_v(pindex), &
                                        & MPI_DOUBLE_PRECISION, pindex, 23, MPI_COMM_WORLD, ierr)

                            END IF

                        END DO send_to_each_process

                        ! deallocate matrices A and B and sch
                        DEALLOCATE(Ag%matrix, Bg%matrix, sch)

                    ENDIF master_proc_scatters_compacts

                    ! each process allocates A and receives the proper chunk of it
                    IF (mycoords(i) == 0) THEN

                      cmp(i,j,k)%A%lb  = [1, (2-j)/2*(3-k) + &
                                             (j/2)*k       + &
                                             (2-j)*j]                ! not intuitive, it exploits the limitations on integer division
                      cmp(i,j,k)%A%lb(1) = cmp(i,j,k)%A%lb(2)        ! this corrects column indices
                      cmp(i,j,k)%A%ub  = [Neq(i,j,k), cmp(i,j,k)%A%lb(2) + Neq(i,j,k) - 1]
                      cmp(i,j,k)%A%ub(1) = cmp(i,j,k)%A%ub(2)        ! this corrects column indices
                      cmp(i,j,k)%A%ld  = 1
                      cmp(i,j,k)%A%ud  = 1
                      cmp(i,j,k)%A%lid = 1
                      cmp(i,j,k)%A%uid = 1
                      ALLOCATE(cmp(i,j,k)%A%matrix(cmp(i,j,k)%A%lb(1):cmp(i,j,k)%A%ub(1), -cmp(i,j,k)%A%ld:+cmp(i,j,k)%A%ud))
                      CALL MPI_RECV(cmp(i,j,k)%A%matrix, &
                                      & 3*Neq(i,j,k),MPI_DOUBLE_PRECISION,0,123,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

                    ELSE

                      cmp(i,j,k)%A%lb  = [1, 1]
                      cmp(i,j,k)%A%ub  = [Neq(i,j,k), Neq(i,j,k)]
                      cmp(i,j,k)%A%ld  = 1
                      cmp(i,j,k)%A%ud  = 1
                      cmp(i,j,k)%A%lid = 1
                      cmp(i,j,k)%A%uid = 1
                      ALLOCATE(cmp(i,j,k)%A%matrix(cmp(i,j,k)%A%lb(1):cmp(i,j,k)%A%ub(1), -cmp(i,j,k)%A%ld:+cmp(i,j,k)%A%ud))
                      CALL MPI_RECV(cmp(i,j,k)%A%matrix, &
                                      & 3*Neq(i,j,k),MPI_DOUBLE_PRECISION,0,123,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

                    END IF

                    ! each process allocates B and receives the proper chunk of it
                    IF (mycoords(i) == 0) THEN ! process at the 'minus' border

                        cmp(i,j,k)%B%ld  = -Bband_int(1)
                        cmp(i,j,k)%B%ud  = +Bband_tot(2)
                        cmp(i,j,k)%B%lid = -Bband_int(1)
                        cmp(i,j,k)%B%uid = +Bband_int(2)
                        cmp(i,j,k)%B%lb  = [1, k-1]
                        cmp(i,j,k)%B%lb(1) = cmp(i,j,k)%A%lb(2)        ! this corrects column indices
                        cmp(i,j,k)%B%ub  = [Neq(i,j,k), Neq(i,j,k) + cmp(i,j,k)%B%uid + cmp(i,j,k)%B%lid - 1 &
                                                        - logical2integer(j==1)*logical2integer(k==1)        &
                                                        - logical2integer(j==2)*logical2integer(k==1)] ! FIXME: una modifica è richiesta dalla derivata 2a per via del duplice significato
                                                                                                       ! dell'indice di staggering 2 (f2f e f2c); l'altra, invece, è dovuta all'inclusione
                                                                                                       ! dei valori di bordo tra le incognite della derivazione prima (ANOMALO).
                        cmp(i,j,k)%B%ub(1) = cmp(i,j,k)%A%ub(2)        ! this corrects column indices
                        ALLOCATE(cmp(i,j,k)%B%matrix(cmp(i,j,k)%B%lb(1):cmp(i,j,k)%B%ub(1), -cmp(i,j,k)%B%ld:+cmp(i,j,k)%B%ud))
                        CALL MPI_RECV(cmp(i,j,k)%B%matrix, (Bband_tot(2) - Bband_int(1) + 1)*Neq(i,j,k), &
                                        & MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    ELSE IF (mycoords(i) == dims(i) - 1) THEN ! process at the 'plus' border

                        cmp(i,j,k)%B%ld  = -Bband_tot(1)
                        cmp(i,j,k)%B%ud  = +Bband_int(2)
                        cmp(i,j,k)%B%lid = -Bband_int(1)
                        cmp(i,j,k)%B%uid = +Bband_int(2)
                        cmp(i,j,k)%B%lb  = [1,       1-cmp(i,j,k)%B%lid]
                        cmp(i,j,k)%B%ub  = [Neq(i,j,k), Neq(i,j,k)+1-logical2integer(j==1)*logical2integer(k==1)]
                        ALLOCATE(cmp(i,j,k)%B%matrix(cmp(i,j,k)%B%lb(1):cmp(i,j,k)%B%ub(1), -cmp(i,j,k)%B%ld:+cmp(i,j,k)%B%ud))
                        CALL MPI_RECV(cmp(i,j,k)%B%matrix, (Bband_int(2) - Bband_tot(1) + 1)*Neq(i,j,k), &
                                        & MPI_DOUBLE_PRECISION, 0, 22, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    ELSE ! process in the center

                        cmp(i,j,k)%B%ld  = -Bband_int(1)
                        cmp(i,j,k)%B%ud  = +Bband_int(2)
                        cmp(i,j,k)%B%lid = -Bband_int(1)
                        cmp(i,j,k)%B%uid = +Bband_int(2)
                        cmp(i,j,k)%B%lb  = [1,       1-cmp(i,j,k)%B%ld]
                        cmp(i,j,k)%B%ub  = [Neq(i,j,k), Neq(i,j,k)+cmp(i,j,k)%B%ud]
                        ALLOCATE(cmp(i,j,k)%B%matrix(cmp(i,j,k)%B%lb(1):cmp(i,j,k)%B%ub(1), -cmp(i,j,k)%B%ld:+cmp(i,j,k)%B%ud))
                        CALL MPI_RECV(cmp(i,j,k)%B%matrix, (Bband_int(2) - Bband_int(1) + 1)*Neq(i,j,k), &
                                        & MPI_DOUBLE_PRECISION, 0, 23, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    END IF

                    cmp(i,j,k)%N = Neq(i,j,k)

                END DO for_each_staggering

            END DO for_each_order_of_derivation

        END DO for_each_direction

        ! deallocation of what's no more used
        IF (myid == 0) THEN
            DEALLOCATE(nequ_v, pcoord)
        END IF
        DEALLOCATE(all_schemes)

    END SUBROUTINE set_compacts


    SUBROUTINE compact_coeffs(LHS_s,RHS_s,LHS_o,RHS_o,c,coeffs)
    ! This subroutine determines the coefficients of a scheme collocated in c and
    ! with the left hand side stencil, where derivatives of order LHS_o are unknown,
    ! defined by LHS_s, and with the right hand side stencil, where derivatives
    ! of order RHS_o are known, defined by RHS_s. The output is in the array
    ! coeffs.
    !
    ! For instance, the following call
    !
    !  CALL compact_coeffs([0.0,1.0],[-.5,0.0,1.0,2.0,3.0],[2,2],[1,0,0,0,0],0.0,coeffs)
    !
    ! Gives the coefficients of the compact scheme whose stencil is the
    ! following
    !
    !  coeffs(2)  coeffs(3)     coeffs(4)     coeffs(5)     coeffs(6)
    !        \      o             o             o             o
    !        |______|_____________|_____________|_____________|
    !               |             |
    !              \\            \\                   IV ordine
    !               1           coeffs(1)
    !

        IMPLICIT NONE

        ! input/output variables
        REAL,    DIMENSION(:),              INTENT(IN)  :: LHS_s,RHS_s ! LHS and RHS stencils
        INTEGER, DIMENSION(:),              INTENT(IN)  :: LHS_o,RHS_o ! LHS and RHS order of derivation
        REAL,                               INTENT(IN)  :: c           ! collocation point of the stencil/scheme
        REAL,    DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: coeffs      ! coefficients of the scheme

        ! internal variables
        INTEGER :: ncoeff, i, ic, sz, k, ii, kk, info
        REAL    :: xmin,xc
        INTEGER, DIMENSION(:),   ALLOCATABLE :: d
        REAL,    DIMENSION(:),   ALLOCATABLE :: x,temp,b,xL,xR,ipiv
        REAL,    DIMENSION(:,:), ALLOCATABLE :: Ab
        REAL,    DIMENSION(:,:), ALLOCATABLE :: A

        ! --------------------------------------------------
        sz = SIZE(LHS_s) + SIZE(RHS_s)

        ALLOCATE(x(sz))
        ALLOCATE(d(sz))
        ALLOCATE(temp(sz))
        ALLOCATE(xL(sz),xR(sz))
        ! ALLOCATE(Ab(sz-1,sz),A(sz-1,sz-1),b(sz-1),coeffs(sz-1))
        ALLOCATE(Ab(sz-1,sz),A(sz-1,sz-1),b(sz-1))
        ALLOCATE(ipiv(sz-1))

        xmin = MIN(MINVAL(LHS_s),MINVAL(RHS_s))
        xL = LHS_s - xmin + 1
        xR = RHS_s - xmin + 1
        xc = c - xmin + 1

        ncoeff = sz - 1
        DO i = 1,SIZE(xL)
            IF (xL(i) == xc) ic = i
        END DO

        x(:SIZE(xL))   = xL
        x(SIZE(xL)+1:) = xR

        d(:SIZE(xL))   = LHS_o
        d(SIZE(xL)+1:) = RHS_o

        temp = 0.0

        ! costruzione del sistema
        DO k = 0,ncoeff-1
            DO ii = 1,sz
                temp(ii) = PRODUCT([ (kk, kk=k,k-d(ii)+1,-1) ])
            END DO
            temp(1:SIZE(xL)) = -temp(1:SIZE(xL))
            Ab(k+1,:) = (x**(k - d))*temp
        END DO

        b = - Ab(:,ic)
        A(:,1:ic-1) = Ab(:,1:ic-1); A(:,ic:sz-1) = Ab(:,ic+1:sz)

        CALL DGESV(sz-1,1,A,sz-1,ipiv,b,sz-1,info)

        coeffs = b

        RETURN

    END SUBROUTINE compact_coeffs


    SUBROUTINE calc_compact_matrices(xL, xR, Adummy, Bdummy, schemes)
    ! The present procedure fills the matrices Adummy and Bdummy of the compact
    ! schemes, based on the LHS mesh xL and RHS mesh xR, using boundary and
    ! central schemes described by the variable schemes (cf. comments below).

        USE MPI, ONLY : MPI_COMM_WORLD

        IMPLICIT NONE

        ! in/out/inout variables
        REAL,                DIMENSION(:),   ALLOCATABLE, INTENT(IN)    :: xL, xR         ! LHS and RHS meshes
        REAL,                DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Adummy, Bdummy ! LHS and RHS matrix
        CHARACTER(len = 20), DIMENSION(:),   ALLOCATABLE, INTENT(IN)    :: schemes        ! array of strings identifying the schemes (see below)

        ! The variable 'schemes' is an array of type CHARACTER(len = 20), that is,
        ! its elements are strings 20 characters long. The string is supposed to
        ! contain the tag of the scheme (e.g. it can contain ____________ to point
        ! to the eigth order compact scheme, ...). This variable is an array, so
        ! it can contain more strings. Such strings (i.e. the elements of the
        ! array 'schemes') are supposed to contain the tag of the the different
        ! schemes used on different locations of the mesh (boundary schemes,
        ! near-boundary schemes, internal schemes, and so on).
        ! In particular 'schemes' is supposed to have a non-positive lower bound and
        ! a non-negative upper bound; schemes(0) should point to the interior scheme,
        ! schemes(-1) to the most internal left boundary scheme, schemes(+1) to
        ! the most internal right boundary scheme, schemes(-2) to the second-to-most
        ! internal left boundary scheme, schemes(+2) to the second-to-most
        ! internal right boundary scheme, and so on...
        ! See the following scheme relative to the case of 3 boundary schemes on
        ! the left and two on the right.
        !  index of schemes( ) used   -3 -2 -1  0  0  0  0 . . . .  0   0   0   +1  +2
        !   collocation point          1  2  3  4  5  6  7 . . . . N-4 N-3 N-2  N-1  N

        ! internal variables
        INTEGER :: ierr, errorcode      ! MPI error
        INTEGER :: i        ! generic indices
        INTEGER :: nLs, nRs ! number of Left and Right boundary schemes
        ! INTEGER :: m, n      ! number of knowns, number of unknowns; Adummy matrix is
                             ! sized (m,m), Bdummy matrix is sized (m,n).
        INTEGER :: mfirst, mlast
        REAL, DIMENSION(:), ALLOCATABLE :: coeffs ! coefficients of the scheme

        nLs = -lbound(schemes, 1)
        nRs = +ubound(schemes, 1)

        mfirst = lbound(Adummy, 1)
        mlast  = ubound(Adummy, 1)

        ! initialize to 0 both matrices
        Adummy = 0; Bdummy = 0

        ! determination of coefficients for internal points
        internal_scheme: SELECT CASE (schemes(0))

            CASE ('D0_c2f_6C_C')
                ! 6-th order compact interpolation scheme from centers to faces
                ! The array of the centers is supposed to be completed with the
                ! first and last faces.
                ! The array of faces doesn't contain the first and last faces.
                !   o   o       o       o                 o       o   o
                !   |___|_______|_______|__  .....  ______|_______|___|
                !           |       |                 |       |
                !           o       o                 o       o

                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-2:i+1),[0,0,0],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1), 1.0, coeffs(2)]
                    Bdummy(i,-2:+1) = coeffs(3:)
                END DO
                DEALLOCATE(coeffs)


            CASE ('D0_f2c_6C_C')
                ! 6-th order compact interpolation scheme from faces to centers
                !   o       o       o                o       o       o
                !   |_______|_______|_  .....  ______|_______|_______|
                !       |       |                |       |       |
                !       o       o                o       o       o

                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-1:i+2),[0,0,0],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1),1.0,coeffs(2)]
                    Bdummy(i,-1:+2) = coeffs(3:)
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_c2f_2E_C')
                !              o       o
                !              |_______|
                !                  |
                !                  \
                ALLOCATE(coeffs(2))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i:i),xR(i-1:i),[1],[0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1: 0) = coeffs
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_c2f_4E_C')
                !              o       o       o       o
                !              |_______|_______|_______|
                !                          |
                !                          \
                ALLOCATE(coeffs(4))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+1),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+1) = coeffs
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_c2f_6E_C')
                ! 6-th order explicit first derivative scheme from centers to faces
                !      o       o       o       o       o       o
                !      |_______|_______|_______|_______|_______|
                !                          |
                !                          \
                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i:i),xR(i-3:i+2),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-3:+2) = coeffs
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_f2c_2E_C')
                !              o       o
                !              |_______|
                !                  |
                !                  \
                ALLOCATE(coeffs(2))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i:i),xR(i:i+1),[1],[0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+1) = coeffs
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_f2c_4E_C')
                !              o       o       o       o
                !              |_______|_______|_______|
                !                          |
                !                          \
                ALLOCATE(coeffs(4))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+2),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+2) = coeffs
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_f2c_6E_C')
                !      o       o       o       o       o       o
                !      |_______|_______|_______|_______|_______|
                !                          |
                !                          \
                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+3),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+3) = coeffs
                END DO
                DEALLOCATE(coeffs)

            CASE ('D1_f2c_6C_C')
                ! 6-th order compact first derivative scheme from faces to centers
                !   o       o       o                o       o       o
                !   |_______|_______|_  .....  ______|_______|_______|
                !       |       |                |       |       |
                !       \       \                \       \       \

                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-1:i+2),[1,1,1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1),1.0,coeffs(2)]
                    Bdummy(i,-1:+2) = coeffs(3:)
                END DO
                DEALLOCATE(coeffs)


            CASE ('D1_c2f_6C_C')
                ! 6-th order compact first derivative scheme from centers to faces
                ! The array of the centers is supposed to be completed with the
                ! first and last faces.
                ! The array of the faces includes the first and last faces (the
                ! unknowns on these points are not used).
                !   o   o       o                o       o   o
                !   |___|_______|_  .....  ______|_______|___|
                !   |       |                |       |       |
                !   \       \                \       \       \

                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-2:i+1),[1,1,1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1),1.0,coeffs(2)]
                    Bdummy(i,-2:+1) = coeffs(3:)
                END DO
                DEALLOCATE(coeffs)


            CASE ('D2_c2c_6C_C','D2_f2f_6C_C')
                ! 6-th order compact second derivative collocated scheme (from
                ! faces to faces)
                !   o       o       o                   o       o       o
                !   |_______|_______|____  .....  ______|_______|_______|
                !           |       |                   |       |
                !           \\      \\                  \\      \\

                ALLOCATE(coeffs(6))
                coeffs = 0.0
                DO i = mfirst+nLs, mlast-nRs
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-2:i+2),[2,2,2],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1),1.0,coeffs(2)]
                    Bdummy(i,-2:+2) = coeffs(3:)
                END DO
                DEALLOCATE(coeffs)


            CASE DEFAULT

                PRINT *, 'Error in compact_mod.f90: no such internal scheme defined.'
                PRINT *, schemes(0)
                CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                STOP

        END SELECT internal_scheme


        ! determination of coefficients for points near the left boundary
        left_points: DO i = mfirst,mfirst+nLs-1

            left_schemes: SELECT CASE (schemes(i-mfirst+1-nLs-1))

                CASE ('D1_f2c_2E_C')
                    !              o       o
                    !              |_______|
                    !                  |
                    !                  \
                    ALLOCATE(coeffs(2))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+1),[1],[0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+1) = coeffs
                    DEALLOCATE(coeffs)

                CASE ('D1_f2c_4E_C')
                    !              o       o       o       o
                    !              |_______|_______|_______|
                    !                          |
                    !                          \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+2),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+2) = coeffs
                    DEALLOCATE(coeffs)

                CASE ('D1_c2f_4E_C')
                    !              o       o       o       o
                    !              |_______|_______|_______|
                    !                          |
                    !                          \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+1),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+1) = coeffs
                    DEALLOCATE(coeffs)

                CASE ('D1_f2c_3E_L')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+3),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+3) = coeffs(:)
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_4E_L')
                    ! 4-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o
                    !      |_______|_______|_______|_______|
                    !          |
                    !          \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+4),[1],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)    = 1
                    Bdummy(i,0:+4) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_5E_L')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                  |
                    !                  \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+4),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+4) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_5E_LL')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !          |
                    !          \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+5),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)    = 1
                    Bdummy(i,0:+5) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_6E_L')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                  |
                    !                  \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+5),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+5) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_6E_LL')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !          |
                    !          \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+6),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)   = 1
                    Bdummy(i,0:6) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_1E_L')
                    ! 1-th order explicit first derivative scheme from faces to centers
                    !           o       o
                    !        ___|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(2))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+1),[1],[0,0],xL(i),coeffs)
                    Adummy(i,0)    = 1
                    Bdummy(i,0:+1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_2E_L')
                    ! 2-th order explicit first derivative scheme from faces to centers
                    !           o       o       o
                    !        ___|_______|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(3))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+2),[1],[0,0,0],xL(i),coeffs)
                    Adummy(i,0)    = 1
                    Bdummy(i,0:+2) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_3E_L')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+2),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i, 0)    = 1.0
                    Bdummy(i,-1:+2) = coeffs(:)
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_3E_LL')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o
                    !        ___|_______|_______|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+3),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+3) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4E_L')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !           |_______|_______|_______|_______|
                    !               |
                    !               \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+3),[1],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+3) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4E_LL')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !        ___|_______|_______|_______|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+4),[1],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+4) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_5E_L')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                       |
                    !                       \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+3),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+3) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_5E_LL')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !               |
                    !               \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+4),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+4) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_5E_LLL')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !        ___|_______|_______|_______|_______|_______|
                    !       |
                    !       \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+5),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+5) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_6E_L')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                  |
                    !                  \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+4),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+4) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_6E_LL')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !          |
                    !          \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+5),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+5) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_6E_LLL')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !      o       o       o       o       o       o       o
                    !   ___|_______|_______|__IF (mycoords(i) == 0) THEN_____|_______|_______|_______|
                    !  |
                    !  \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+6),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)   = 1
                    Bdummy(i,0:6) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D0_c2f_4C_L')
                    !   o   o       o
                    !   |___|_______|____
                    !           |       |
                    !           o       o
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i+1),xR(i-2:i),[0,0],[0,0,0],xL(i),coeffs)
                    Adummy(i, 0:+1) = [1.0,coeffs(1)]
                    Bdummy(i,-2: 0) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('D0_f2c_5C_L')
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !       |       |
                    !       o       o
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i+1),xR(i:i+3),[0,0],[0,0,0,0],xL(i),coeffs)
                    Adummy(i, 0:+1) = [1.0,coeffs(1)]
                    Bdummy(i, 0:+3) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4C_L')
                    !       o       o
                    !   ____|_______|____
                    !   |       |       |
                    !   \       \       \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-1:i),[1,1,1],[0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1),1.0,coeffs(2)]
                    Bdummy(i,-1: 0) = coeffs(3:)
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4C_LL')
                    !   o   o       o       o
                    !   |___|_______|_______|
                    !   |       |
                    !   \       \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i+1),xR(i-1:i+2),[1,1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i, 0:+1) = [1.0,coeffs(1)]
                    Bdummy(i,-1:+2) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('D2_c2c_5C_L','D2_f2f_5C_L')
                    !   o       o       o       o       o       o
                    !   |_______|_______|_______|_______|_______|
                    !           |       |
                    !           \\      \\
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i+1),xR(i-1:i+4),[2,2],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i, 0:+1) = [1.0,coeffs(1)]
                    Bdummy(i,-1:+4) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('zero')
                    ! No boundary scheme
                    Bdummy(i,:) = 0


                CASE DEFAULT

                    PRINT *, 'Error in compact_mod.f90: no such left scheme defined.'
                    print *, schemes(i-mfirst+1-nLs-1)
                    CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                    STOP

            END SELECT left_schemes

        ENDDO left_points

        ! determination of coefficients for points near the right boundary
        right_points: DO i = mlast-nRs+1,mlast

            right_schemes: SELECT CASE (schemes(i-mlast+nRs))

                CASE ('D1_f2c_2E_C')
                    !              o       o
                    !              |_______|
                    !                  |
                    !                  \
                    ALLOCATE(coeffs(2))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i:i+1),[1],[0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i, 0:+1) = coeffs
                    DEALLOCATE(coeffs)

                CASE ('D1_f2c_4E_C')
                    !              o       o       o       o
                    !              |_______|_______|_______|
                    !                          |
                    !                          \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-1:i+2),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-1:+2) = coeffs
                    DEALLOCATE(coeffs)

                CASE ('D1_c2f_4E_C')
                    !              o       o       o       o
                    !              |_______|_______|_______|
                    !                          |
                    !                          \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+1),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+1) = coeffs
                    DEALLOCATE(coeffs)

                CASE ('D1_f2c_3E_R')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !                       |
                    !                       \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-2:i+1),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:+1) = coeffs(:)
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_4E_R')
                    ! 4-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o
                    !      |_______|_______|_______|_______|
                    !                                  |
                    !                                  \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-3:i+1),[1],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-3:+1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_5E_R')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                                  |
                    !                                  \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-3:i+2),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-3:+2) = coeffs
                    DEALLOCATE(coeffs)



                CASE ('D1_f2c_5E_RR')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                                          |
                    !                                          \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-4:i+1),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-4:+1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_6E_R')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                                          |
                    !                                          \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-4:i+2),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-4:+2) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_f2c_6E_RR')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                                                  |
                    !                                                  /
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-5:i+1),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)    = 1
                    Bdummy(i,-5:1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_1E_R')
                    ! 1-th order explicit first derivative scheme from faces to centers
                    !                   o       o
                    !                   |_______|___
                    !                               |
                    !                               \
                    ALLOCATE(coeffs(2))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-2:i-1),[1],[0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-2:-1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_2E_R')
                    ! 2-th order explicit first derivative scheme from faces to centers
                    !           o       o       o
                    !           |_______|_______|___
                    !                               |
                    !                               \
                    ALLOCATE(coeffs(3))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-3:i-1),[1],[0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-3:-1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_3E_R')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !                   o       o       o       o
                    !                   |_______|_______|_______|
                    !                                       |
                    !                                       \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-3:i),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)    = 1
                    Bdummy(i,-3:0) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_3E_RR')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !                   o       o       o       o
                    !                   |_______|_______|_______|___
                    !                                               |
                    !                                               \
                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-4:i-1),[1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-4:-1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4E_R')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !           |_______|_______|_______|_______|
                    !                                       |
                    !                                       \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-4:i),[1],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-4: 0) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4E_RR')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !           |_______|_______|_______|_______|___
                    !                                               |
                    !                                               \
                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-5:i-1),[1],[0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-5:-1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_5E_R')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                                       |
                    !                                       \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-4:i+1),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-4:+1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_5E_RR')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                                               |
                    !                                               \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-5:i),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-5: 0) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_5E_RRR')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|___
                    !                                                       |
                    !                                                       \
                    ALLOCATE(coeffs(6))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-6:i-1),[1],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-6:-1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_6E_R')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !    o       o       o       o       o       o       o
                    !    |_______|_______|_______|_______|_______|_______|
                    !                                        |
                    !                                        \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-5:i+1),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-5:+1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_6E_RR')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !    o       o       o       o       o       o       o
                    !    |_______|_______|_______|_______|_______|_______|
                    !                                                |
                    !                                                \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-6:i),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-6:0) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_6E_RRR')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !    o       o       o       o       o       o       o
                    !    |_______|_______|_______|_______|_______|_______|___
                    !                                                        |
                    !                                                        \
                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i:i),xR(i-7:i-1),[1],[0,0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,0)     = 1
                    Bdummy(i,-7:-1) = coeffs
                    DEALLOCATE(coeffs)


                CASE ('D0_c2f_4C_R')
                    !       o       o   o
                    !    ___|_______|___|
                    !   |       |
                    !   o       o

                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i-1:i),xR(i-1:i+1),[0,0],[0,0,0],xL(i),coeffs)
                    Adummy(i,-1: 0) = [coeffs(1),1.0]
                    Bdummy(i,-1:+1) = coeffs(2:)
                    DEALLOCATE(coeffs)

                CASE ('D0_f2c_5C_R')
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !               |       |
                    !               o       o

                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i-1:i),xR(i-2:i+1),[0,0],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1: 0) = [coeffs(1),1.0]
                    Bdummy(i,-2:+1) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4C_R')
                    !       o       o
                    !   ____|_______|____
                    !   |       |       |
                    !   \       \       \

                    ALLOCATE(coeffs(4))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i-1:i+1),xR(i-1:i),[1,1,1],[0,0],xL(i),coeffs)
                    Adummy(i,-1:+1) = [coeffs(1),1.0,coeffs(2)]
                    Bdummy(i,-1:0) = coeffs(3:)
                    DEALLOCATE(coeffs)


                CASE ('D1_c2f_4C_RR')
                    !   o       o       o   o
                    !   |_______|_______|___|
                    !               |       |
                    !               \       \

                    ALLOCATE(coeffs(5))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i-1:i),xR(i-3:i),[1,1],[0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1: 0) = [coeffs(1),1.0]
                    Bdummy(i,-3:0) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('D2_c2c_5C_R','D2_f2f_5C_R')
                    !   o       o       o       o       o       o
                    !   |_______|_______|_______|_______|_______|
                    !                           |       |
                    !                           \\      \\

                    ALLOCATE(coeffs(7))
                    coeffs = 0.0
                    CALL compact_coeffs(xL(i-1:i),xR(i-4:i+1),[2,2],[0,0,0,0,0,0],xL(i),coeffs)
                    Adummy(i,-1: 0) = [coeffs(1),1.0]
                    Bdummy(i,-4:+1) = coeffs(2:)
                    DEALLOCATE(coeffs)


                CASE ('zero')
                    ! No boundary scheme
                    Bdummy(i,:) = 0


                CASE DEFAULT

                    PRINT *, 'Error in compact_mod.f90: no such right scheme defined.'
                    PRINT *, schemes(i-mlast+nRs)
                    CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                    STOP

            END SELECT right_schemes

        ENDDO right_points

    END SUBROUTINE calc_compact_matrices


END MODULE compact
