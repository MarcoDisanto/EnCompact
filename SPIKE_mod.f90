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
    REAL, DIMENSION(:),   ALLOCATABLE :: DphiC, Dphi0  ! coupling unknowns and knowns (master)
    REAL, DIMENSION(2)                :: DphiCp        ! coupling unknowns of the current process

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
                    delta(1,1) = A%matrix(1,-1) ! ... by alpha...
                    delta(N,2) = A%matrix(N,+1) ! ... and beta
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


    SUBROUTINE SPIKE_exchange

        USE MPI,          ONLY : MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_INTEGER, MPI_DOUBLE_PRECISION
        USE variables, ONLY : uvwp
        USE MPI_module, ONLY : ndims, myid, procs_grid
        USE variables, ONLY : MPI_flats_m, MPI_flats_p, rdispls_m, rdispls_p, sdispls_m, sdispls_p

        IMPLICIT NONE

        INTEGER :: ierr, ii
        INTEGER :: dpsize ! size in byte of the MPI_DOUBLE_PRECISION
        !INTEGER, DIMENSION(2*ndims) :: rdispls, sdispls

        !IF (myid == 13) THEN
        !    CALL MPI_SEND(uvwp(2)%values(uvwp(2)%b(1,2)+uvwp(2)%b_ol(1,1),1,1),1,MPI_flats_m(2,1),22,777,MPI_COMM_WORLD,ierr)
        !    CALL MPI_RECV(uvwp(2)%values(uvwp(2)%b(1,2)+1,1,1),1,MPI_flats_p(2,1),22,888,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        !    CALL MPI_SEND(uvwp(2)%values(1,1,uvwp(2)%b(3,2)+uvwp(2)%b_ol(1,1)),1,MPI_flats_m(2,3),14,666,MPI_COMM_WORLD,ierr)
        !    CALL MPI_RECV(uvwp(2)%values(1,1,uvwp(2)%b(3,2)+1),1,MPI_flats_p(2,3),14,555,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        !ELSE IF (myid == 2) THEN
        !    CALL MPI_RECV(uvwp(3)%values(1,uvwp(3)%b(2,2)+1,1),1,MPI_flats_p(3,2),5,111,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        !    CALL MPI_SEND(uvwp(3)%values(1,uvwp(3)%b(2,2)+uvwp(3)%b_ol(2,1),1),1,MPI_flats_m(3,2),5,222,MPI_COMM_WORLD,ierr)

        !ELSE IF (myid == 5) THEN
        !    CALL MPI_SEND(uvwp(3)%values(1,uvwp(3)%b(2,1),1),1,MPI_flats_p(3,2),2,111,MPI_COMM_WORLD,ierr)
        !    CALL MPI_RECV(uvwp(3)%values(1,uvwp(3)%b_ol(2,1),1),1,MPI_flats_m(3,2),2,222,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        !ELSE IF (myid == 22) THEN
        !    CALL MPI_RECV(uvwp(2)%values(uvwp(2)%b_ol(1,1),1,1),1,MPI_flats_m(2,1),13,777,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        !    CALL MPI_SEND(uvwp(2)%values(uvwp(2)%b(1,1),1,1),1,MPI_flats_p(2,1),13,888,MPI_COMM_WORLD,ierr)

        !ELSE IF (myid == 14) THEN
        !    CALL MPI_RECV(uvwp(2)%values(1,1,uvwp(2)%b_ol(3,1)),1,MPI_flats_m(2,3),13,666,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        !    CALL MPI_SEND(uvwp(2)%values(1,1,uvwp(2)%b(3,1)),1,MPI_flats_p(2,3),13,555,MPI_COMM_WORLD,ierr)
        !END IF


        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dpsize, ierr)

        !rdispls = [size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)), &
        !     
        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,2) - uvwp(2)%b_bo(1,1) + 1), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,2) - uvwp(2)%b_bo(2,1) + 1) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,2) - uvwp(2)%b_bo(3,1) + 1) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1))]*dpsize

        !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !IF (myid == 5) THEN
        !    print *, "uvwp(2)%values' shape", shape(uvwp(2)%values)
        !    print *, "uvwp(2)%values' b    ", uvwp(2)%b
        !    print *, "uvwp(2)%values' b_ol ", uvwp(2)%b_ol
        !    print *, "uvwp(2)%values' b_bc ", uvwp(2)%b_bc
        !    print *, "uvwp(2)%values' b_bo ", uvwp(2)%b_bo
        !    print *, rdispls/dpsize
        !END IF
        !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        !sdispls = [size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &
        !     
        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & - (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,2) - uvwp(2)%b_bo(2,1) + 1) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,1) - uvwp(2)%b_bo(3,1)) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1)), &

        !         & size(uvwp(2)%values,1)*size(uvwp(2)%values,2)*(uvwp(2)%b(3,2) - uvwp(2)%b_bo(3,1) + 1) + &
        !         & size(uvwp(2)%values,1)*(uvwp(2)%b(2,1) - uvwp(2)%b_bo(2,1)) + &
        !         & (uvwp(2)%b(1,1) - uvwp(2)%b_bo(1,1))]*dpsize

        !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !IF (myid == 5) THEN
        !    print *, "uvwp(2)%values' shape", shape(uvwp(2)%values)
        !    print *, "uvwp(2)%values' b    ", uvwp(2)%b
        !    print *, "uvwp(2)%values' b_ol ", uvwp(2)%b_ol
        !    print *, "uvwp(2)%values' b_bc ", uvwp(2)%b_bc
        !    print *, "uvwp(2)%values' b_bo ", uvwp(2)%b_bo
        !    print *, sdispls/dpsize
        !END IF
        !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ! IF (myid == 6) THEN
            ! print *, [sdispls_p(2,1), sdispls_m(2,1), sdispls_p(2,2), sdispls_m(2,2), sdispls_p(2,3), sdispls_m(2,3)]/dpsize
            ! print *, [rdispls_m(2,1), rdispls_p(2,1), rdispls_m(2,2), rdispls_p(2,2), rdispls_m(2,3), rdispls_p(2,3)]/dpsize
        ! END IF

        DO ii = 0, 26
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
            IF (myid == ii) THEN
                print *,myid,[sdispls_p(2,1), sdispls_m(2,1), sdispls_p(2,2), sdispls_m(2,2), sdispls_p(2,3), sdispls_m(2,3)]/dpsize
                ! print *,myid,[rdispls_m(2,1), rdispls_p(2,1), rdispls_m(2,2), rdispls_p(2,2), rdispls_m(2,3), rdispls_p(2,3)]/dpsize
                ! print *, 'prima', myid
                ! CALL MPI_NEIGHBOR_ALLTOALLW(uvwp(2)%values(uvwp(2)%b_bo(1,1),uvwp(2)%b_bo(2,1),uvwp(2)%b_bo(3,1)), &
                    ! & [1,1,1,1,1,1], &
                    ! & [sdispls_p(2,1), sdispls_m(2,1), sdispls_p(2,2), sdispls_m(2,2), sdispls_p(2,3), sdispls_m(2,3)],&
                    ! & [MPI_flats_p(2,1), MPI_flats_m(2,1), MPI_flats_p(2,2), MPI_flats_m(2,2), MPI_flats_p(2,3), MPI_flats_m(2,3)],&
                                          ! & uvwp(2)%values(uvwp(2)%b_bo(1,1),uvwp(2)%b_bo(2,1),uvwp(2)%b_bo(3,1)), &
                    ! & [1,1,1,1,1,1], &
                    ! & [rdispls_m(2,1), rdispls_p(2,1), rdispls_m(2,2), rdispls_p(2,2), rdispls_m(2,3), rdispls_p(2,3)],&
                    ! & [MPI_flats_m(2,1), MPI_flats_p(2,1), MPI_flats_m(2,2), MPI_flats_p(2,2), MPI_flats_m(2,3), MPI_flats_p(2,3)],&
                                           ! &procs_grid, ierr)                       
                ! print *, 'dopo', myid
            END IF
        END DO
                                                                            
        ! MPI_NEIGHBOR_ALLTOALLW(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPES,
        !                        RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPES, COMM, IERROR)
        !<type>    SENDBUF(*), RECVBUF(*)
        !INTEGER    SENDCOUNTS(*), SDISPLS(*), SENDTYPES(*)
        !INTEGER    RECVCOUNTS(*), RDISPLS(*), RECVTYPES(*)
        !INTEGER    COMM, IERROR
        !sendbuf
        !       Starting address of send buffer.
        !sendcounts
        !       Integer array, where entry i specifies the number of elements to send to neighbor i.
        !sdispls
        !       Integer array, where entry i specifies the displacement (in bytes, offset from sendbuf) from which to send data to neighbor i.
        !sendtypes
        !       Datatype array, where entry i specifies the datatype to use when sending data to neighbor i.
        !recvcounts
        !       Integer array, where entry j specifies the number of elements to receive from neighbor j.
        !rdispls
        !       Integer array, where entry j specifies the displacement (in bytes, offset from recvbuf) to which data from neighbor j should be written.
        !recvtypes
        !       Datatype array, where entry j specifies the datatype to use when receiving data from neighbor j.
        !comm
        !       Communicator over which data is to be exchanged.

    END SUBROUTINE SPIKE_exchange


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
