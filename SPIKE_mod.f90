MODULE SPIKE
! This module contains the following procedures
!    SUBROUTINE SPIKE_init
!    SUBROUTINE SPIKE_exchange

    USE bandedmatrix, ONLY : CDS

    IMPLICIT NONE

    TYPE SPIKE_type
        REAL, DIMENSION(:,:), ALLOCATABLE :: a1end ! first and last column of the inverse of A
        TYPE(CDS)                         :: C
    END TYPE SPIKE_type

    TYPE(SPIKE_type), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: SPK

    REAL, DIMENSION(:,:), ALLOCATABLE :: a1end         ! first and last columns of the inverse of A
    ! REAL, DIMENSION(:),   ALLOCATABLE :: DphiC, Dphi0  ! coupling unknowns and knowns (master)
    ! REAL, DIMENSION(2)                :: DphiCp        ! coupling unknowns of the current process

CONTAINS

    SUBROUTINE SPIKE_init
    ! This procedure "initializes" the SPIKE, in that it calculates the first and
    ! last columns of the inverse of A. Furthermore the last and first elements
    ! of these columns (4 scalars) are sent to the "pencil"-master process which
    ! arranges them in the coupling matrix C.

        USE MPI
        USE Thomas_suite
        USE MPI_module, ONLY : ndims, dims, mycoords, periodic, pencil_comm
        USE essentials, ONLY : log2int => logical2integer, printmatrix
        USE compact,    ONLY : cmp
        USE, INTRINSIC :: IEEE_ARITHMETIC

        IMPLICIT NONE

        !internal variables
        INTEGER :: i, j, k, ierr
        REAL, DIMENSION(:,:), ALLOCATABLE :: delta ! first and last column of the identity matrix
        TYPE(CDS),            POINTER :: A     => NULL()
        INTEGER,              POINTER :: N     => NULL()
        TYPE(CDS),            POINTER :: C     => NULL()
        REAL, DIMENSION(:,:), POINTER :: a1end => NULL()

        INTEGER                           :: MPI_zig, MPI_zigzag , MPI_zigzag_res ! handles to MPI derived data types
        INTEGER                           :: dpsize ! size in byte of the MPI_DOUBLE_PRECISION
        INTEGER(KIND = MPI_ADDRESS_KIND)  :: start, extent ! lower bound and extent used in MPI_TYPE_CREATE_RESIZED

        ! Define
        REAL :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
        NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

        ALLOCATE(SPK(ndims, 0:2, 2))

        for_each_direction: DO i = 1, ndims

            for_each_order_of_derivation: DO j = 0, 2

                for_each_staggering: DO k = 1, 2

                    ! Set the pointers of already allocated variables
                    A => cmp(i,j,k)%A
                    N => cmp(i,j,k)%N

                    ! Calculate the spikes (i.e. first and last column of the inverse of A)
                    ALLOCATE(SPK(i,j,k)%a1end(N, 2))! Allocate spikes
                    a1end => SPK(i,j,k)%a1end   ! pointer to the spikes
                    ALLOCATE(delta(N, 2))       ! first and last column of the identity matrix
                    delta = 0                   ! ... already multiplied...
                    delta(1,1) = A%matrix(A%lb(1),-1) ! ... by alpha...
                    delta(N,2) = A%matrix(A%ub(1),+1) ! ... and beta
                    CALL Thomas(A%matrix, delta(:,1), a1end(:,1)) ! solution for left  spike
                    CALL Thomas(A%matrix, delta(:,2), a1end(:,2)) ! solution for right spike
                    DEALLOCATE(delta)           ! deallocation

                    ! Allocate the coupling matrix
                    pencil_master_allocates_C: IF (mycoords(i) == 0) THEN

                        ! ALLOCATE(SPK(i,j,k)%C%matrix(2*(dims(i) - log2int(.NOT.periodic(i))),-2:+2))
                        ALLOCATE(SPK(i,j,k)%C%matrix(2*dims(i),-2:+2))

                    ELSE
                        ! Puppet allocation for non-master.
                        ! Note that even if this array is never used
                        ! by non-master processes, it has to make sense
                        ! in MPI calls. As a consequence I cannot leave
                        ! it unallocated, nor allocate it with
                        !       ALLOCATE(SPK(i,j,k)%C%matrix(1,1))
                        ! since in the MPI_GATHER it is referenced by
                        !       C%matrix(2,-2)
                        ! and, in turn, at least the element (2,-2) has to be
                        ! present. In light of this, the proper allocation is
                        ALLOCATE(SPK(i,j,k)%C%matrix(2:2,-2:-2))

                    END IF pencil_master_allocates_C

                    ! Set the pointer to the coupling matrix
                    C => SPK(i,j,k)%C

                    ! To build the coupling matrix, three MPI derived data types are defined,
                    ! each built upon the preceding, as follows:
                    !    - the first type (MPI_zig, which is temporary) consists of two scalars
                    !      positioned in the C matrix in the positions represented by x
                    !      in the following pattern (the 'i' represents the "insertion point",
                    !      i.e. where the subsequent element would be inserted):
                    !
                    !                                 x
                    !                                xi
                    !
                    CALL MPI_TYPE_VECTOR(2, 1, 2*dims(i)-1, MPI_DOUBLE_PRECISION, MPI_zig, ierr)
                    !
                    !    - the second type (MPI_zigzag, which is also temporary) is a sequence
                    !      of two MPI_zig elements properly stridden, so as to
                    !      be positioned in the C matrix in these positions:
                    !
                    !                                 x  x
                    !                                x  xi
                    !
                    CALL MPI_TYPE_VECTOR(2, 1, 3, MPI_zig, MPI_zigzag, ierr)
                    !
                    !    - the third type (MPI_zigzag_res), differs from MPI_zigzag only in that
                    !      its extent is changed in order to allow the insertion of the following
                    !      element in the same column where the preceding has started, so as
                    !      to change the pattern in the following:
                    !
                    !                                 x  x
                    !                                x  x
                    !
                    !                                i
                    !
                    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dpsize, ierr)
                    extent = 2*dpsize; start = 0
                    CALL MPI_TYPE_CREATE_RESIZED(MPI_zigzag, start, extent, MPI_zigzag_res, ierr)
                    ! Only the third type is committed
                    CALL MPI_TYPE_COMMIT(MPI_zigzag_res, ierr)

                    ! The coupling matrix Cdummy is built
                    pencil_master_build_C_matrix: IF (mycoords(i) == 0) THEN

                        ! initialization of the coupling matrix
                        C%matrix      = 0.0
                        C%matrix(:,0) = 1.0
                        C%ld = 2; C%ud = 2
                        ! C%lb = [1, 1]; C%ub = [1, 1]*2*(dims(i) - log2int(.NOT.periodic(i)))
                        C%lb = [1, 1]; C%ub = [1, 1]*2*dims(i)

                    ELSE
                        ! puppet filling of C in other processes
                        C%matrix = NaN
                        C%ld = NaN; C%ud = NaN
                        C%lb = NaN; C%ub = NaN

                    END IF pencil_master_build_C_matrix

                    ! pencil-master process gathers the 4 values (base and tips
                    ! of spikes) from each process on the same pencil and puts
                    ! them in two rows corresponding to each process (included
                    ! itself)
                    CALL MPI_GATHER([a1end(N,1), a1end(1,1), a1end(N,2), a1end(1,2)], 4, MPI_DOUBLE_PRECISION, &
                                 & C%matrix(2,-2), 1, MPI_zigzag_res, 0, pencil_comm(i), ierr)

                    ! free the used type
                    CALL MPI_TYPE_FREE(MPI_zigzag_res, ierr)

                    ! TO DO:
                    ! in the case of non-periodic direction, first and last rows
                    ! of C should not be used; for now I simply insert NaNs for
                    ! pencil-master
                    non_periodic_direction: IF (.NOT.periodic(i) .AND. mycoords(i) == 0) THEN
                        C%matrix(1,:)         = NaN
                        C%matrix(2*dims(i),:) = NaN
                    END IF non_periodic_direction

                END DO for_each_staggering

            END DO for_each_order_of_derivation

        END DO for_each_direction

    END SUBROUTINE SPIKE_init


    SUBROUTINE SPIKE_exchange_uvw

        USE MPI,        ONLY : MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_PROC_NULL
        USE MPI_module, ONLY : ndims, myid, procs_grid, idm, idp
        USE variables,  ONLY : uvwp, MPI_flats_m, MPI_flats_p
        USE essentials, ONLY : log2int => logical2integer, KronDelta

        IMPLICIT NONE

        INTEGER :: ierr, ii, ic, id

        ! boundary values of are exchanged between adjacent processes

        for_each_component: DO ic = 1, ndims

            id = 1
            CALL MPI_SENDRECV(uvwp(ic)%values( &
                            & uvwp(ic)%b(id,2) + uvwp(ic)%b_ol(id,1) - KronDelta(ic,id)*log2int(idm(id) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_m(ic,id), idp(id), 111, &
                            & uvwp(ic)%values( &
                            & uvwp(ic)%b_bo(id,1), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_m(ic,id), idm(id), 111, &
                            & procs_grid, MPI_STATUS_IGNORE, ierr)
            CALL MPI_SENDRECV(uvwp(ic)%values( &
                            & uvwp(ic)%b(id,1), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_p(ic,id), idm(id), 1111, &
                            & uvwp(ic)%values( &
                            & uvwp(ic)%b(id,2) + 1, &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_p(ic,id), idp(id), 1111, &
                            & procs_grid, MPI_STATUS_IGNORE, ierr)

            id = 2
            CALL MPI_SENDRECV(uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & uvwp(ic)%b(id,2) + uvwp(ic)%b_ol(id,1) - KronDelta(ic,id)*log2int(idm(id) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_m(ic,id), idp(id), 111, &
                            & uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & uvwp(ic)%b_bo(id,1), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_m(ic,id), idm(id), 111, &
                            & procs_grid, MPI_STATUS_IGNORE, ierr)
            CALL MPI_SENDRECV(uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & uvwp(ic)%b(id,1), &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_p(ic,id), idm(id), 1111, &
                            & uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & uvwp(ic)%b(id,2) + 1, &
                            & 1 + KronDelta(ic,3)*log2int(idm(3) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_p(ic,id), idp(id), 1111, &
                            & procs_grid, MPI_STATUS_IGNORE, ierr)

            id = 3
            CALL MPI_SENDRECV(uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & uvwp(ic)%b(id,2) + uvwp(ic)%b_ol(id,1) - KronDelta(ic,id)*log2int(idm(id) == MPI_PROC_NULL)), &
                            & 1, MPI_flats_m(ic,id), idp(id), 111, &
                            & uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & uvwp(ic)%b_bo(id,1)), &
                            & 1, MPI_flats_m(ic,id), idm(id), 111, &
                            & procs_grid, MPI_STATUS_IGNORE, ierr)
            CALL MPI_SENDRECV(uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & uvwp(ic)%b(id,1)), &
                            & 1, MPI_flats_p(ic,id), idm(id), 1111, &
                            & uvwp(ic)%values( &
                            & 1 + KronDelta(ic,1)*log2int(idm(1) == MPI_PROC_NULL), &
                            & 1 + KronDelta(ic,2)*log2int(idm(2) == MPI_PROC_NULL), &
                            & uvwp(ic)%b(id,2) + 1), &
                            & 1, MPI_flats_p(ic,id), idp(id), 1111, &
                            & procs_grid, MPI_STATUS_IGNORE, ierr)

            END DO for_each_component

    END SUBROUTINE SPIKE_exchange_uvw






SUBROUTINE SPIKE_solve(id, ider, istag, psi, q)
! This subroutine solves the tridiagonal system necessary to retrieve the derivatives
! It finds the necessary matrices by itself, accessing compact_mod by means of id, ider
! and istag.

USE compact, ONLY: compact_type, cmp
USE essentials ! dopo cancella
USE bandedmatrix
USE MPI_module
USE Thomas_suite
USE, INTRINSIC :: IEEE_ARITHMETIC

IMPLICIT NONE

INTEGER, INTENT(IN)                           :: id, ider, istag    ! indices for direction, order of derivative and staggering (see compact_mod)
REAL, DIMENSION(:), ALLOCATABLE, INTENT(IN)   :: q                  ! velocity vector in input
REAL, DIMENSION(:), INTENT(OUT)               :: psi                ! derivative, solution of the system in output
TYPE(CDS), POINTER                            :: Ap, Bp             ! just to simplify notations
REAL, DIMENSION(SIZE(psi))                    :: bq                 ! RHS vector (B*q), same size as psi
REAL, DIMENSION(:), ALLOCATABLE               :: psiC               ! vectors of coupled unknowns (only master process has to allocate it)
REAL, DIMENSION(:), ALLOCATABLE               :: psi0C              ! RHS of reduced system (only master process has to allocate it)
INTEGER                                       :: FirstLast_type     ! MPI type that contains first and last elements of a vector, ONLY: dims, ierr, ONLY: dims, ierr
REAL, DIMENSION(:, :), POINTER                :: Cp                 ! SPIKE matrix pointer (just to simplify notation)
INTEGER                                       :: ierr
REAL                                          :: dum                ! temporary storage used to support elements reordering
INTEGER                                       :: i

! Define NaN
REAL :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

Ap => cmp(id, ider, istag)%A
Bp => cmp(id, ider, istag)%B

! Dimensions check
IF      (Ap%ub(2)-Ap%lb(2)+1/=SIZE(psi)) THEN
  PRINT *, 'ERROR in SPIKE_mod. Output vector dimensions not coherent'
  STOP
ELSE IF (Bp%ub(2)-Bp%lb(2)+1/=SIZE(q))   THEN
  PRINT *, 'ERROR in SPIKE_mod. Input vector dimensions not coherent'
  STOP
END IF

!!!!!!!!! Partial solution (often denoted as psi0) !!!!!!!!!!
bq = Bp*q
CALL Thomas(Ap%matrix, bq, psi)

!!!!!!!!!! Master process gathers required data !!!!!!!!!!
! Type creation
CALL MPI_TYPE_VECTOR(2, 1, SIZE(psi)-1, MPI_DOUBLE_PRECISION, FirstLast_type, ierr)
CALL MPI_TYPE_COMMIT(FirstLast_type, ierr)
! Reduced system vectors allocation
ALLOCATE(psiC(2*dims(id)))
ALLOCATE(psi0C(2*dims(id)))
! master gathers RHS terms consisting of first and last elements of partial solution vector
CALL MPI_gather(psi, 1, FirstLast_type, psi0C, 2, MPI_DOUBLE_PRECISION, 0, pencil_comm(id), ierr)
psiC = NaN ! debug-purpose allocation

!!!!!!!!!! Master process solves the penta-diagonal system !!!!!!!!!!
coupled_system: IF (mycoords(id)==0) THEN
  Cp => SPK(id, ider, istag)%C%matrix(2:SIZE(SPK(id, ider, istag)%C%matrix, 1)-1, :)  ! NOTE : no need for NaN rows
  CALL pentdag(Cp(:, 1), Cp(:, 2), Cp(:, 3), Cp(:, 4), Cp(:, 5), &
               psi0C(2:SIZE(psi0C)-1), psiC(2:SIZE(psiC)-1), SIZE(psiC)-2) ! must exclude first and last equations and unknowns to avoid generating more NaNs
END IF coupled_system

!!!!!!!!!! Reordering the reduced system's solution vector !!!!!!!!!!
! Simple swap algorithm
DO i = 2, SIZE(psiC)-2, 2
  dum = psiC(i)
  psiC(i) = psiC(i+1)
  psiC(i+1) = dum
END DO

!!!!!!!!!! Master process distributes required data !!!!!!!!!!
CALL MPI_SCATTER(psiC, 2, MPI_DOUBLE_PRECISION, psiC, 2, MPI_DOUBLE_PRECISION, 0, pencil_comm(id), ierr)

!!!!!!!!!! Update unknowns !!!!!!!!!!
IF (mycoords(id)==0) THEN
  ! add only right value
  psi = psi - SPK(id, ider, istag)%a1end(:, 2)*psiC(2)
ELSE IF (mycoords(id)==dims(id)-1) THEN
  ! add only left value
  psi = psi - SPK(id, ider, istag)%a1end(:, 1)*psiC(1)
ELSE
  ! add both left and right values
  psi = psi - SPK(id, ider, istag)%a1end(:, 1)*psiC(1) &
            - SPK(id, ider, istag)%a1end(:, 2)*psiC(2)
END IF

END SUBROUTINE SPIKE_solve







SUBROUTINE SPIKE_solve2(id, ider, istag, PSI)
! This is an evolution of previous "SPIKE_solve" subroutine, that works on blocks instead
! of lines of nodes. Its role is to correct the previously calculated (input argument
! PSI) solution of the block-diagonal system. Compared to the previous version, this
! one is much faster since it communicates a bunch of values with a single call,
! thus using the BUS as little as possible.

USE compact, ONLY: compact_type, cmp
USE essentials
USE bandedmatrix
USE MPI_module
USE Thomas_suite
USE, INTRINSIC :: IEEE_ARITHMETIC

IMPLICIT NONE

INTEGER, INTENT(IN)                           :: id, ider, istag    ! indices for direction, order of derivative and staggering (see compact_mod)
REAL, DIMENSION(:, :, :), INTENT(INOUT)       :: PSI                ! solution block to be updated after calculating coupling terms
REAL, DIMENSION(:, :, :), ALLOCATABLE         :: PSIC               ! vectors of coupled unknowns (only master process has to allocate it)
REAL, DIMENSION(:, :, :), ALLOCATABLE         :: PSI0C              ! RHS of reduced system (only master process has to allocate it)
REAL, DIMENSION(:, :), POINTER                :: Cp                 ! SPIKE matrix pointer (just to simplify notation)
REAL, DIMENSION(:, :), ALLOCATABLE            :: dum                ! temporary storage used to support elements reordering
INTEGER                                       :: i, j, k            ! just indices used to move along 3 directions
INTEGER, DIMENSION(ndims)                     :: dmn                ! dimensions used for various allocations
! type construction stuff
INTEGER                                       :: MPI_1slab_type, MPI_1slab_res, MPI_2slab_type
INTEGER                                       :: MPI_gather_type, MPI_gather_type_res
INTEGER, DIMENSION(ndims)                     :: array_of_sizes
INTEGER, DIMENSION(ndims)                     :: array_of_subsizes
INTEGER, DIMENSION(ndims)                     :: array_of_starts
INTEGER(MPI_ADDRESS_KIND)                     :: lb, ext            ! lower bound and extent of MPI types
INTEGER                                       :: dblsz              ! size in bytes of a DOUBLE PRECISION type
INTEGER                                       :: ierr

! Define NaN
REAL :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

!!!!!!!!!! Type creation and allocations !!!!!!!!!!
! Sender side (slave processes)
array_of_sizes = SHAPE(PSI)
array_of_starts = [0, 0, 0]
DO i = 1, ndims
  array_of_subsizes(i) = KronDelta(i, id) + (1-KronDelta(i, id))*SIZE(PSI, i)
END DO
! 1 element thick slab definition
CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                            & array_of_sizes, array_of_subsizes, array_of_starts, &
                            & MPI_ORDER_FORTRAN, &
                            & MPI_DOUBLE_PRECISION, MPI_1slab_type, &
                            & ierr)
! the following resizing deceives the code into believing the type to be shorter
! than it actually is. Thank MPI guys for this s**t
CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblsz, ierr)
lb = 0
ext = dblsz*( KronDelta(id, 1) + KronDelta(id, 2)*SIZE(PSI, 1) + KronDelta(id, 3)*SIZE(PSI, 1)*SIZE(PSI, 2) )
CALL MPI_TYPE_CREATE_RESIZED(MPI_1slab_type, lb, ext, MPI_1slab_res, ierr)
! put together first and last plates of the block
CALL MPI_TYPE_VECTOR(2, 1, SIZE(PSI, id)-1, MPI_1slab_res, MPI_2slab_type, ierr)
CALL MPI_TYPE_COMMIT(MPI_2slab_type, ierr)

! Receiver side (master process)
! only the master process needs all the values of psiC (reduced system solution);
! the slave ones require only 2 values (or plates in this case)
DO i = 1, ndims
  dmn(i) = KronDelta(i, id)*dims(id)*2*logical2integer(mycoords(id)==0) + &
           KronDelta(i, id)*2*logical2integer(mycoords(id)/=0) + &
          (1-KronDelta(i, id))*SIZE(PSI, i)
END DO
ALLOCATE(PSIC(dmn(1), dmn(2), dmn(3)))
ALLOCATE(PSI0C(dmn(1), dmn(2), dmn(3)))

! cannot recycle previous plate type, for PSIC is smaller than PSI
array_of_sizes = SHAPE(PSI0C)
DO i = 1, ndims
  array_of_subsizes(i) = KronDelta(i, id) + (1-KronDelta(i, id))*SIZE(PSI, i)
END DO
! 1 element thick slab definition
CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
                            & array_of_sizes, array_of_subsizes, array_of_starts, &
                            & MPI_ORDER_FORTRAN, &
                            & MPI_DOUBLE_PRECISION, MPI_gather_type, &
                            & ierr)
! another deception like the previous one
CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblsz, ierr)
lb = 0
ext = dblsz*( KronDelta(id, 1) + KronDelta(id, 2)*SIZE(PSI, 1) + KronDelta(id, 3)*SIZE(PSI, 1)*SIZE(PSI, 2) )
CALL MPI_TYPE_CREATE_RESIZED(MPI_gather_type, lb, ext, MPI_gather_type_res, ierr)
CALL MPI_TYPE_COMMIT(MPI_gather_type_res, ierr)

! Debug-purpose initialization
PSI0C = NaN
PSIC = NaN


!!!!!!!!!! DATA gathering !!!!!!!!!!
CALL MPI_gather(PSI, 1, MPI_2slab_type, PSI0C, 2, MPI_gather_type_res, 0, pencil_comm(id), ierr)


! just to simplify notation
Cp => SPK(id, ider, istag)%C%matrix(2:SIZE(SPK(id, ider, istag)%C%matrix, 1)-1, :)  ! NOTE : no need for NaN rows

! depending on direction the DO cycle must be different
SELECT CASE(id)

  CASE(1) !!!!!!!!!!!!!!!!!!!! x direction !!!!!!!!!!!!!!!!!!!!

    master_only: IF (mycoords(id)==0) THEN
      !!!!!!!!!! 5-diag system resolution !!!!!!!!!!
      DO j = 1, SIZE(PSI0C, 2)
        DO k = 1, SIZE(PSI0C, 3)
          ! must exclude first and last equations and unknowns to avoid generating more NaNs
          CALL pentdag(Cp(:, 1), Cp(:, 2), Cp(:, 3), Cp(:, 4), Cp(:, 5), &
                       PSI0C(2:SIZE(PSI0C, id)-1, j, k), PSIC(2:SIZE(PSIC, id)-1, j, k), SIZE(PSIC, id)-2)
        END DO
      END DO

      !!!!!!!!!! Reordering the reduced system's solution vector !!!!!!!!!!
      ! Simple swap algorithm
      ALLOCATE(dum(SIZE(PSIC, 2), SIZE(PSIC, 3)))
      DO i = 2, SIZE(PSIC, id)-2, 2
        dum = PSIC(i, :, :)
        PSIC(i, :, :) = PSIC(i+1, :, :)
        PSIC(i+1, :, :) = dum
      END DO

    END IF master_only

    !!!!!!!!!!! Master process distributes required data !!!!!!!!!!
    CALL MPI_SCATTER(PSIC, 2, MPI_gather_type_res, &    ! sender side
    PSIC, 2, MPI_gather_type_res, &    ! receiver side
    0, pencil_comm(id), ierr)

    !!!!!!!!!! Update solution !!!!!!!!!!
    IF (mycoords(id)==0) THEN
      ! add only right value
      DO j = 1, SIZE(PSI, 2)
        DO k = 1, SIZE(PSI, 3)
          PSI(:, j, k) = PSI(:, j, k) - SPK(id, ider, istag)%a1end(:, 2)*PSIC(2, j, k)
        END DO
      END DO

    ELSE IF (mycoords(id)==dims(id)-1) THEN
      ! add only left value
      DO j = 1, SIZE(PSI, 2)
        DO k = 1, SIZE(PSI, 3)
          PSI(:, j, k) = PSI(:, j, k) - SPK(id, ider, istag)%a1end(:, 1)*PSIC(1, j, k)
        END DO
      END DO

    ELSE
      ! add both left and right values
      DO j = 1, SIZE(PSI, 2)
        DO k = 1, SIZE(PSI, 3)
          PSI(:, j, k) = PSI(:, j, k) - SPK(id, ider, istag)%a1end(:, 1)*PSIC(1, j, k) &
                                      - SPK(id, ider, istag)%a1end(:, 2)*PSIC(2, j, k)
        END DO
      END DO
    END IF


  CASE(2) !!!!!!!!!!!!!!!!!!!! y direction !!!!!!!!!!!!!!!!!!!!

    master_only2: IF (mycoords(id)==0) THEN

      !!!!!!!!!! 5-diag system resolution !!!!!!!!!!
      DO j = 1, SIZE(PSI0C, 1)
        DO k = 1, SIZE(PSI0C, 3)
          ! must exclude first and last equations and unknowns to avoid generating more NaNs
          CALL pentdag(Cp(:, 1), Cp(:, 2), Cp(:, 3), Cp(:, 4), Cp(:, 5), &
                       PSI0C(j, 2:SIZE(PSI0C, id)-1, k), PSIC(j, 2:SIZE(PSIC, id)-1, k), SIZE(PSIC, id)-2)
        END DO
      END DO

      !!!!!!!!!! Reordering the reduced system's solution vector !!!!!!!!!!
      ! Simple swap algorithm
      ALLOCATE(dum(SIZE(PSIC, 1), SIZE(PSIC, 3)))
      DO i = 2, SIZE(PSIC, id)-2, 2
        dum = PSIC(:, i, :)
        PSIC(:, i, :) = PSIC(:, i+1, :)
        PSIC(:, i+1, :) = dum
      END DO

    END IF master_only2

    !!!!!!!!!!! Master process distributes required data !!!!!!!!!!
    CALL MPI_SCATTER(PSIC, 2, MPI_gather_type_res, &    ! sender side
                     PSIC, 2, MPI_gather_type_res, &    ! receiver side
                     0, pencil_comm(id), ierr)

    !!!!!!!!!! Update solution !!!!!!!!!!
    IF (mycoords(id)==0) THEN
      ! add only right value
      DO j = 1, SIZE(PSI, 1)
        DO k = 1, SIZE(PSI, 3)
          PSI(j, :, k) = PSI(j, :, k) - SPK(id, ider, istag)%a1end(:, 2)*PSIC(j, 2, k)
        END DO
      END DO

    ELSE IF (mycoords(id)==dims(id)-1) THEN
      ! add only left value
      DO j = 1, SIZE(PSI, 1)
        DO k = 1, SIZE(PSI, 3)
          PSI(j, :, k) = PSI(j, :, k) - SPK(id, ider, istag)%a1end(:, 1)*PSIC(j, 1, k)
        END DO
      END DO

    ELSE
      ! add both left and right values
      DO j = 1, SIZE(PSI, 1)
        DO k = 1, SIZE(PSI, 3)
          PSI(j, :, k) = PSI(j, :, k) - SPK(id, ider, istag)%a1end(:, 1)*PSIC(j, 1, k) &
                                      - SPK(id, ider, istag)%a1end(:, 2)*PSIC(j, 2, k)
        END DO
      END DO
    END IF


  CASE(3) !!!!!!!!!!!!!!!!!!!! z direction !!!!!!!!!!!!!!!!!!!!

    master_only3: IF (mycoords(id)==0) THEN

      !!!!!!!!!! 5-diag system resolution !!!!!!!!!!
      DO j = 1, SIZE(PSI0C, 1)
        DO k = 1, SIZE(PSI0C, 2)
          ! must exclude first and last equations and unknowns to avoid generating more NaNs
          CALL pentdag(Cp(:, 1), Cp(:, 2), Cp(:, 3), Cp(:, 4), Cp(:, 5), &
                       PSI0C(j, k, 2:SIZE(PSI0C, id)-1), PSIC(j, k, 2:SIZE(PSIC, id)-1), SIZE(PSIC, id)-2)
        END DO
      END DO

      !!!!!!!!!! Reordering the reduced system's solution vector !!!!!!!!!!
      ! Simple swap algorithm
      ALLOCATE(dum(SIZE(PSIC, 1), SIZE(PSIC, 2)))
      DO i = 2, SIZE(PSIC, id)-2, 2
        dum = PSIC(:, :, i)
        PSIC(:, :, i) = PSIC(:, :, i+1)
        PSIC(:, :, i+1) = dum
      END DO

    END IF master_only3

    !!!!!!!!!!! Master process distributes required data !!!!!!!!!!
    CALL MPI_SCATTER(PSIC, 2, MPI_gather_type_res, &    ! sender side
                     PSIC, 2, MPI_gather_type_res, &    ! receiver side
                     0, pencil_comm(id), ierr)

    !!!!!!!!!! Update solution !!!!!!!!!!
    IF (mycoords(id)==0) THEN
      ! add only right value
      DO j = 1, SIZE(PSI, 1)
        DO k = 1, SIZE(PSI, 2)
          PSI(j, k, :) = PSI(j, k, :) - SPK(id, ider, istag)%a1end(:, 2)*PSIC(j, k, 2)
        END DO
      END DO

    ELSE IF (mycoords(id)==dims(id)-1) THEN
      ! add only left value
      DO j = 1, SIZE(PSI, 1)
        DO k = 1, SIZE(PSI, 2)
          PSI(j, k, :) = PSI(j, k, :) - SPK(id, ider, istag)%a1end(:, 1)*PSIC(j, k, 1)
        END DO
      END DO

    ELSE
      ! add both left and right values
      DO j = 1, SIZE(PSI, 1)
        DO k = 1, SIZE(PSI, 2)
          PSI(j, k, :) = PSI(j, k, :) - SPK(id, ider, istag)%a1end(:, 1)*PSIC(j, k, 1) &
                                      - SPK(id, ider, istag)%a1end(:, 2)*PSIC(j, k, 2)
        END DO
      END DO
    END IF

END SELECT


CALL MPI_TYPE_FREE(MPI_1slab_type, ierr)
CALL MPI_TYPE_FREE(MPI_1slab_res, ierr)
CALL MPI_TYPE_FREE(MPI_2slab_type, ierr)
CALL MPI_TYPE_FREE(MPI_gather_type, ierr)
CALL MPI_TYPE_FREE(MPI_gather_type_res, ierr)


END SUBROUTINE SPIKE_solve2













!    SUBROUTINE SPIKE_step(Adummy,Bdummy,phidummy,Dphidummy)
!
!        USE MPI_module
!        USE bandedmatrix
!        USE Thomas_suite
!
!        IMPLICIT NONE
!
!        ! in/out/inout variables
!        REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)    :: Adummy
!        TYPE(CDS),                         INTENT(IN)    :: Bdummy
!        REAL, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: phidummy
!        REAL, DIMENSION(N),                INTENT(OUT)   :: Dphidummy
!
!        ! internal variables
!        INTEGER :: i
!        REAL, DIMENSION(N) :: q
!        REAL, DIMENSION(N) :: Dphi0
!
!        ! boundary values of phi are exchanged between adjacent processes in
!        ! order to compute q
!        CALL MPI_SENDRECV(phidummy(N-Bdummy%lb+1:N),  Bdummy%lb,MPI_DOUBLE_PRECISION,idR,123,&
!                        & phidummy(:0),  Bdummy%lb,MPI_DOUBLE_PRECISION,idL,123,&
!                        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!        CALL MPI_SENDRECV(phidummy(1:Bdummy%ub),      Bdummy%ub,MPI_DOUBLE_PRECISION,idL,123,&
!                        & phidummy(N+1:),Bdummy%ub,MPI_DOUBLE_PRECISION,idR,123,&
!                        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!
!        ! computation of q based on the exchanged values
!        ! print *, lbound(phidummy), ubound(phidummy)
!        q = Bdummy*phidummy
!        ! print *, 'ciaoo'
!        ! stop
!
!        ! solution of the block-diagonal system (provisional solution)
!        CALL Thomas(Adummy, q, Dphidummy)
!
!        ! construction of the rhs of the coupling system
!        CALL MPI_GATHER([Dphidummy(1), Dphidummy(N)], &
!                          & 2, MPI_DOUBLE_PRECISION, &
!                          & Dphi0, &
!                          & 2, MPI_DOUBLE_PRECISION, &
!                          & 0, MPI_COMM_WORLD, ierr)
!        ! This is the input to the pentadiagonal solver and contains the provisional
!        ! block-to-block values of Dphi, Dphi0 in the order first of proc 0, last of
!        ! proc 0, first of proc 1, last of proc 1, first of proc 2, ... and so on.
!
!        ! Solution of the reduced coupling system
!        IF (myid == 0) THEN
!                CALL cypent(C(:,-2), C(:,-1), C(:,-0), C(:,+1), C(:,+2), Dphi0, &
!                          & C(1,-2), C(1,-1), C(2,-2), C(2*dims(i)-1,+2), C(2*dims(i),+1), C(2*dims(i),+2), &
!                          & DphiC, 2*dims(i))
!                ! The output contains the exact values of the unknowns in the order:
!                ! first of proc 0, last of proc 0, first of proc 1, last of proc 1,
!                ! first of proc 2, ... and so on (the same as the input array).
!                ! These values should be sent by master process to other processes
!                ! (included the master itself). The task is accomplished using
!                ! MPI_SCATTER. This routine, however, splits the message (i.e. DphiC)
!                ! into continuous segments (of size 2) and sends them in rank order.
!                ! So it's necessary to edit the vector DphiC such that each couple of
!                ! odd-(odd+1) elments is substituted by the preceding and following
!                ! elements.
!                !
!                ! odd-index values are stored in the temporary array
!                temp1D = DphiC(1:2*dims(i)-1:2)
!                ! the odd-index values are substituted with the (odd-1)-index values
!                ! (cyclically)
!                DphiC(1:2*dims(i)-1:2) = DphiC([2*dims(i),(i, i = 2,2*dims(i)-2,2)])
!                ! the even-index values are substituted with the (even+1)-index values
!                ! that were stored in the temporary array
!                DphiC(2:2*dims(i):2)   = temp1D([(i, i = 2,dims(i)),1])
!         END IF
!
!        ! master process 'scatters' the array DphiC across the processes
!        CALL MPI_SCATTER(DphiC, 2, MPI_DOUBLE_PRECISION, DphiCp, 2, &
!                       & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!
!        ! provisional values of Dphi are corrected adding the SPIKES
!        Dphidummy   = Dphidummy - a1end(:,1)*DphiCp(1) - a1end(:,2)*DphiCp(2)
!
!    END SUBROUTINE SPIKE_step




!    ! LA SEGUENTE PROCEDURA NON È PIÙ UTILIZZATA NEL CASO PARALLELO,
!    ! in quanto inglobata nella procedura SPIKE_step, che applica tutto
!    ! l'algoritmo SPIKE (opportunamente inizializzato con la procedura SPIKE_init).
!    !
!    ! È invece usata nel caso seriale, in cui lo scambio di icognite deve
!    ! comunque essere effettuato per tenere conto della periodicità
!    SUBROUTINE SPIKE_q(phidummy,Bdummy,qdummy)
!    ! This procedure performs the exchange of the values of phi at block-to-block
!    ! boundaries and then correctly computes the product B*y.
!
!        USE MPI_module
!        ! USE compact
!        USE bandedmatrix
!
!        IMPLICIT NONE
!
!        REAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: phidummy
!        TYPE(CDS),                       INTENT(IN)    :: Bdummy
!        REAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: qdummy
!
!        CALL MPI_SENDRECV(phidummy(N),  1,MPI_DOUBLE_PRECISION,idR,123,&
!                        & phidummy(0),  1,MPI_DOUBLE_PRECISION,idL,123,&
!                        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!        CALL MPI_SENDRECV(phidummy(1),  1,MPI_DOUBLE_PRECISION,idL,123,&
!                        & phidummy(N+1),1,MPI_DOUBLE_PRECISION,idR,123,&
!                        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!
!        ! computation of q based on the exchanged values
!        qdummy = Bdummy*phidummy
!
!    END SUBROUTINE SPIKE_q

END MODULE SPIKE
