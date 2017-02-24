MODULE set_pressure
! This module contains the following procedures
!    SUBROUTINE Set_p
!    SUBROUTINE set_pres_weig
!    SUBROUTINE build_GraDiv_glob(i,k)
!    SUBROUTINE spread_CDS(Matr_loc,i,k)
!    SUBROUTINE set_var
!    SUBROUTINE Weigc_gr

    USE MPI_module
    USE essentials
    USE compact, ONLY : calc_compact_matrices
    USE library
    USE grids
    USE variables
    USE bandedmatrix, ONLY : CDS, OPERATOR(*)

    IMPLICIT NONE

    TYPE pressure_schemes
        CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE :: sch
    END TYPE pressure_schemes
    TYPE(pressure_schemes), DIMENSION(:,:),   ALLOCATABLE :: pres_schemes ! structure containing all the schemes used
    !First index== direction, Second Index=saggering (1 = Grad, 2 = Div)

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: OD, b, weigC_grid, res

    INTEGER, DIMENSION(3,2) :: indv
    INTEGER, DIMENSION(3,2) :: ovl_rec, ovl_sen ! ANOMALO: beware the -5 and +5 alias PORCATE

    REAL, DIMENSION(:,:,:), POINTER :: p

    TYPE P_comm
        REAL, DIMENSION (:,:,:), POINTER :: p_m, p_p
    END TYPE

    TYPE(P_comm), ALLOCATABLE, DIMENSION (:) :: p_rec, p_sen

    CHARACTER(4) :: It_type
    REAL :: tol

    !!!!!!AGGIUNTE!!!!
    TYPE(CDS) :: Bg
    Type(CDS), ALLOCATABLE, DIMENSION(:,:) :: GraDiv, GraDiv_glob
    !First index == direction, Second Index = saggering (1 = Grad, 2 = Div)

    Type(CDS), ALLOCATABLE, DIMENSION(:) :: Lapl, Lapl_glob
    !Array dimension = ndims

    INTEGER, DIMENSION(2) :: Bband_tot ! lower and upper band (of the banded matrix B) for all the schemes
    INTEGER, DIMENSION(2) :: Bband_int ! lower and upper band (of the banded matrix B) for the central scheme

CONTAINS

    SUBROUTINE Set_p

        IMPLICIT NONE

        CALL set_pres_weig

        CALL set_var

    END SUBROUTINE Set_p


    SUBROUTINE set_pres_weig
    ! Alfonsina l'ha scritta modificando set_compacts del modulo compact_mod.f90

        IMPLICIT NONE

        INTEGER :: i,k

        IF(myid == 0)THEN
            ALLOCATE(GraDiv_glob(ndims,2))
        END IF

        ALLOCATE(GraDiv(ndims, 2))
        ALLOCATE(Lapl(ndims))

        DO k = 1, 2
            DO i = 1, ndims
                CALL build_GraDiv_glob(i,k)
                CALL spread_CDS(GraDiv(i,k),i, k)
            END DO
        END DO

        DO i = 1, ndims
            IF (myid == 0) THEN
                !!!!!!!!!!!!!!!!!!!!!CALCOLO PESI LAPLACIANO GLOBALE!!!!!!!!!!!!
                ! GraDiv_glob(i,1)%matrix(        1,:) = 0
                ! GraDiv_glob(i,1)%matrix(Ntot(i)+1,:) = 0
                Bg = GraDiv_glob(i,2)*GraDiv_glob(i,1)
                ! Bg%matrix(1,0) = Bg%matrix(1,0) - 1
                ! Bg%matrix(Ntot(i),0) = Bg%matrix(Ntot(i),0) - 1
                !IF (i == 1) THEN
                !    print *, ''
                !    print *, 'divergenza'
                !    PRINT *, GraDiv_glob(i,2)%lb, GraDiv_glob(i,2)%ub
                !!    CALL printmatrix(GraDiv_glob(i,2)%matrix)
                !    print *, ''
                !    print *, 'gradiente'
                !    PRINT *, GraDiv_glob(i,1)%lb, GraDiv_glob(i,1)%ub
                !!    CALL printmatrix(GraDiv_glob(i,1)%matrix)
                !    print *, ''
                !    print *, 'laplaciano'
                !    PRINT *, Bg%lb, Bg%ub
                !!    CALL printmatrix(Bg%matrix)
                !END IF
                ! TO DO: le quattro seguenti assegnazioni vanno eliminate assieme alla
                ! modifica da fare al modulo compact_mod.f90 e alla routine spread_CDS e
                ! build_GraDiv_glob del presente modulo.
                ! Inoltre il prodotto tra matrici CDS (funzione CDS_mat_x_CDS_mat del
                ! modulo bandedmatrix_mod.f90) assegna automaticamente la banda
                ! totale della matrice prodotto, ma non quella interna. Dovrei
                ! codificare l'assegnazione di questa informazione ad un nuovo campo
                ! del tipo CDS, in maniera tale da poter eliminare le PORCATE.
                Bband_tot(1) = -Bg%ld
                Bband_tot(2) = +Bg%ud
                Bband_int(1) = -1!5 ! PORCATE
                Bband_int(2) = +1!5 ! PORCATE

            END IF

        !!!!!!!Comunico Pesi Laplaciano
        CALL spread_CDS(Lapl(i), i, 2)
        END DO

        IF(myid == 0)THEN
            DEALLOCATE(GraDiv_glob)
            DEALLOCATE(pres_schemes)
        END IF


    END SUBROUTINE set_pres_weig

    !L'idea è quella di mettere il k fuori!!!! così per ogni k fornito, calcola, e scattera!!!

    SUBROUTINE build_GraDiv_glob(i,k)

        USE essentials, ONLY: logical2integer

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: i, k

        REAL, ALLOCATABLE, DIMENSION(:) :: grid_appoL, grid_appoR
        CHARACTER(len = 20), DIMENSION(:), ALLOCATABLE :: sch
        INTEGER :: j = 1
        REAL                            :: h ! spacing needed to create fictional faces and cells in the periodic case

        TYPE(CDS) :: Afake ! this is bound to be the identity matrix, since schemes are explicit

        IF (myid == 0) THEN

            j = 1

            ! ALLOCATE(sch(___:___)) ! sembra che l'allocazione sia automatica....
            sch = pres_schemes(i,k)%sch  !!!!!!!!!!!!

            ! Qui Alfonsina utilizza degli array di appoggio, invece di useare
            ! il mio array all_grids. La motivazione è che all_grids( , ,1,2) ha
            ! due celle fittizie in più: la cella 0, che coincide con la prima
            ! faccia, e la cella N(i)+1, che coincide con l'ultima faccia.
            ! Negli altri 3 casi, cioè all_grids( , ,1,1), all_grids( , ,2,1) e
            ! all_grids( , ,2,2), le griglie coincidono.
            IF (k == 1) THEN
                grid_appoL = all_grids(i,j,k,1)%g(1:Ntot(i)+1-logical2integer(periodic(i)))!Facce
                grid_appoR = all_grids(i,j,k,2)%g(1:Ntot(i))!Celle
                ! periodicity correction (mesh widening)
                IF (periodic(i)) THEN
                  h = all_grids(i,j,k,2)%g(2)-all_grids(i,j,k,2)%g(1)
                  grid_appoL = [grid_appoL, grid_appoL(Ntot(i)+1-logical2integer(periodic(i)))+h]
                  grid_appoR = [grid_appoR(1)-h, grid_appoR, grid_appoR(Ntot(i))+h]
                END IF
                !print *, 'k =', k, 'grid_appoL =', grid_appoL
                !print *, 'k =', k, 'grid_appoR =', grid_appoR
                !print *, LBOUND(grid_appoL), UBOUND(grid_appoL)
                !print *, LBOUND(grid_appoR), UBOUND(grid_appoR)
            ELSE
                grid_appoL = all_grids(i,j,k,1)%g!Celle
                grid_appoR = all_grids(i,j,k,2)%g!Facce
                ! periodicity correction (mesh widening)
                IF (periodic(i)) THEN
                  h = all_grids(i,j,k,2)%g(2)-all_grids(i,j,k,2)%g(1)
                  grid_appoR = [grid_appoR, grid_appoR(Ntot(i))+h]
                  !print *, 'k =', k, 'grid_appoR =', grid_appoR
                  !print *, LBOUND(grid_appoL), UBOUND(grid_appoL)
                  !print *, LBOUND(grid_appoR), UBOUND(grid_appoR)
                END IF
            ENDIF

            ! lower and upper bounds of full counterparts of matrices A and B along the two dimensions:
            ! - both matrices A and B have the rows indexed as the unknowns, that is
            !   the lower and upper indices are those of the LHS grid, that is ( , , ,1)
            ! - matrix A is square, so the columns are indexed as the rows
            ! - matrix B has columns indexed as the the knowns, that is the lower and
            !   upper indices are those of the RHS grid, that is ( , , ,2)
            Afake%lb = [lbound(grid_appoL), lbound(grid_appoL)]
            Afake%ub = [ubound(grid_appoL), ubound(grid_appoL)]
            Bg%lb    = [lbound(grid_appoL), lbound(grid_appoR)]
            Bg%ub    = [ubound(grid_appoL), ubound(grid_appoR)]

            IF (k==1 .AND. periodic(i)) THEN
              Bg%lb(2) = Bg%lb(2) - 1
              Bg%ub(2) = Bg%ub(2) - 1

              !PRINT *, Bg%lb
              !PRINT *, Bg%ub
            END IF

            ! lower and upper band of B matrix...
            Bband_tot = detBband(sch)      ! ... considering all  the      rows
            Bband_int = detBband(sch(0:0)) ! ... considering only internal rows
            Bg%ld = -Bband_tot(1) ! number of lower diagonals (positive if lower diagonals are present)
            Bg%ud = +Bband_tot(2) ! number of upper diagonals (positive if upper diagonals are present)
            ! allocate matrices A and B
            ALLOCATE(Afake%matrix(Afake%lb(1):Afake%ub(1),      0:0))
            ALLOCATE(Bg%matrix(Bg%lb(1):Bg%ub(1),          -Bg%ld:+Bg%ud))

            ! fill the matrices A and B
            CALL calc_compact_matrices(grid_appoL, grid_appoR, Afake%matrix, Bg%matrix, sch, periodic(i))

            DEALLOCATE(Afake%matrix,sch, grid_appoL, grid_appoR)
            GraDiv_glob(i,k) = Bg

        END IF

    END SUBROUTINE build_GraDiv_glob


    SUBROUTINE spread_CDS(Matr_loc,i,k)

        IMPLICIT NONE

        INTEGER :: j ! generic indices
        INTEGER :: ierr
        INTEGER :: il, iu   ! index of first and last equations of the generic process (il can be different from 1 and similarly iu ...)
        INTEGER                             :: pindex ! generic process' index
        INTEGER, DIMENSION(:), ALLOCATABLE  :: pcoord ! generic process' coordinates in cartesian communicator (not those of the current process, which are stored in the array mycoords of the MPI_module module)
        INTEGER, DIMENSION(:), ALLOCATABLE  :: nequ_v ! array containing the number of equations of the various processes, which changes depending on the varibales considered
        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpistatus
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: nequ_m ! array containing the number of equations of the various processes
        ! The preceding internal variables are replaced time by time depending on the scheme considered
        !!!!!
        !!!!!Aggiunte!!!!!
        TYPE(CDS), INTENT(OUT) ::Matr_loc
        INTEGER, INTENT(IN) ::k,i

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

        j = 1
        !!!!INIZIO COMUNICAZIONI!!!!!!!!!!!!!!!!!!!
        ! the master process broadcasts Bband_int to every process
        CALL MPI_BCAST(Bband_int, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! the master process broadcasts Bband_tot to every process
        CALL MPI_BCAST(Bband_tot, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! the master process gathers the number of equations of all the processes
        CALL MPI_GATHER(Neq(i,j,k), 1, MPI_INTEGER, nequ_v, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        master_process_scatters_matrix: IF (myid == 0) THEN

            ! the master process arranges the number of equations of the
            ! processes in an array sized as procs_grid
            SELECT CASE (ndims)
                CASE (1)
                    nequ_m = RESHAPE(nequ_v, [dims(1:1), 1, 1])
                    PRINT *, 'WARNING in compact_mod.f90: per ndims == 1 non credo di aver fatto un check.'
                    PRINT *, 'Inoltre una Navier-Stokes 1D (non quasi-1D) incompressile ha poco senso...'
                CASE (2)
                    nequ_m = RESHAPE(nequ_v, [dims(1:2), 1], order = [3, 2, 1]) ! NB: questo 'order = ' ci deve stare, perché il compilatore è un figlio di buona donna
                    PRINT *, 'WARNING in compact_mod.f90: per ndims == 2 non credo di aver fatto un check.'
                CASE (3)
                    nequ_m = RESHAPE(nequ_v, dims(1:3), order = [3, 2, 1]) ! NB: questo 'order = ' ci deve stare, perché il compilatore è un figlio di buona donna
            END SELECT

            ! remaining processes' loop  (included the master itself)
            send_to_each_process: DO pindex = 0, nprocs-1

                ! determination of cartesian coordinates of the pindex process
                CALL MPI_CART_COORDS(procs_grid, pindex, ndims, pcoord, ierr)
                ! lower and upper indices of the chunk to be sent to the pindex process
                SELECT CASE (i)
                    CASE (1)
                        il = SUM(nequ_m(0:pcoord(1)-1, pcoord(2), pcoord(3))) + Bg%lb(1)
                    CASE (2)
                        il = SUM(nequ_m(pcoord(1), 0:pcoord(2)-1, pcoord(3))) + Bg%lb(1)
                    CASE (3)
                        il = SUM(nequ_m(pcoord(1), pcoord(2), 0:pcoord(3)-1)) + Bg%lb(1)
                END SELECT

                iu = il + nequ_v(pindex) - 1

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

        END IF master_process_scatters_matrix

        ! each process allocates B and receives the proper chunk of it
        IF (mycoords(i) == 0) THEN ! process at the 'minus' border

            Matr_loc%ld = -Bband_int(1)
            Matr_loc%ud = +Bband_tot(2)
            Matr_loc%lb = [1,       1-Matr_loc%ld]
            Matr_loc%ub = [Neq(i,j,k), Neq(i,j,k)+Matr_loc%ud]

            ALLOCATE(Matr_loc%matrix(Matr_loc%lb(1):Matr_loc%ub(1),-Matr_loc%ld:+Matr_loc%ud))

            CALL MPI_RECV(Matr_loc%matrix, (Bband_tot(2) - Bband_int(1) + 1)*Neq(i,j,k), &
            & MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, mpistatus, ierr)

        ELSE IF (mycoords(i) == dims(i) - 1) THEN ! process at the 'plus' border

            Matr_loc%ld = -Bband_tot(1)
            Matr_loc%ud = +Bband_int(2)
            Matr_loc%lb = [1,       1-Matr_loc%ld]
            Matr_loc%ub = [Neq(i,j,k), Neq(i,j,k)+Matr_loc%ud]

            ALLOCATE(Matr_loc%matrix(Matr_loc%lb(1):Matr_loc%ub(1), -Matr_loc%ld:+Matr_loc%ud))

            CALL MPI_RECV(Matr_loc%matrix, (Bband_int(2) - Bband_tot(1) + 1)*Neq(i,j,k), &
            & MPI_DOUBLE_PRECISION, 0, 22, MPI_COMM_WORLD, mpistatus, ierr)

        ELSE ! process in the center

            Matr_loc%ld = -Bband_int(1)
            Matr_loc%ud = +Bband_int(2)
            Matr_loc%lb = [1,       1-Matr_loc%ld]
            Matr_loc%ub = [Neq(i,j,k), Neq(i,j,k)+Matr_loc%ud]
            ALLOCATE(Matr_loc%matrix(Matr_loc%lb(1):Matr_loc%ub(1), -Matr_loc%ld:+Matr_loc%ud))

            CALL MPI_RECV(Matr_loc%matrix, (Bband_int(2) - Bband_int(1) + 1)*Neq(i,j,k), &
            & MPI_DOUBLE_PRECISION, 0, 23, MPI_COMM_WORLD, mpistatus, ierr)

        END IF

        IF (myid == 0) THEN
            DEALLOCATE(Bg%matrix)
            DEALLOCATE(nequ_v, pcoord)
        END IF

    END SUBROUTINE spread_CDS


    SUBROUTINE set_var

        IMPLICIT NONE

        INTEGER :: i

        indv(:,1) = 1
        indv(:,2) = N(:)

        DO i = 1, ndims
            ovl_rec(i,1) = indv(i,1)-1!5 PORCATE
            ovl_rec(i,2) = indv(i,2)+1!5 PORCATE
            !Dovrebbe essere 'ampiezza' di quello a fianco!
            ovl_sen(i,1) = indv(i,1)-1+1!5 PORCATE
            ovl_sen(i,2) = indv(i,2)+1-1!5 PORCATE
        END DO

        ALLOCATE(uvwp(4)%values(ovl_rec(1,1):ovl_rec(1,2),ovl_rec(2,1):ovl_rec(2,2),ovl_rec(3,1):ovl_rec(3,2)))

        p => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2))

        ALLOCATE(p_rec(ndims), p_sen(ndims))

        If(ndims>= 1) then
            !NOTA :: il PUNTAMENTO a una porzione di matrice "probabilemnte" richiede
            !sempre vettori di indici!!!
            p_rec(1)%p_m => uvwp(4)%values(ovl_rec(1,1):indv(1,1)-1,indv(2,1):indv(2,2), indv(3,1):indv(3,2))
            p_rec(1)%p_p => uvwp(4)%values(indv(1,2)+1:ovl_rec(1,2),indv(2,1):indv(2,2), indv(3,1):indv(3,2))

            !Bordo Ovest
            p_sen(1)%p_m => uvwp(4)%values(indv(1,1):ovl_sen(1,1),indv(2,1):indv(2,2), indv(3,1):indv(3,2))
            !Bordo Est
            p_sen(1)%p_p => uvwp(4)%values(ovl_sen(1,2):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2))
        end if

        If(ndims>= 2) then
            p_rec(2)%p_m => uvwp(4)%values(indv(1,1):indv(1,2),ovl_rec(2,1):indv(2,1)-1, indv(3,1):indv(3,2))
            p_rec(2)%p_p => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,2)+1:ovl_rec(2,2), indv(3,1):indv(3,2))

            !Bordo Down
            p_sen(2)%p_m => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,1):ovl_sen(2,1), indv(3,1):indv(3,2))
            !Bordo Up
            p_sen(2)%p_p => uvwp(4)%values(indv(1,1):indv(1,2),ovl_sen(2,2):indv(2,2), indv(3,1):indv(3,2))
        end if

        if(ndims == 3) then
            p_rec(3)%p_m => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,1):indv(2,2),ovl_rec(3,1):indv(3,1)-1)
            p_rec(3)%p_p => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,2)+1:ovl_rec(3,2))

            !Bordo Back
            p_sen(3)%p_m => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,1):indv(2,2), indv(3,1):ovl_sen(3,1))
            !Bordo Forth
            p_sen(3)%p_p => uvwp(4)%values(indv(1,1):indv(1,2),indv(2,1):indv(2,2), ovl_sen(3,2):indv(3,2))
        end if

        ALLOCATE(OD(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2)))
        ALLOCATE(res(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2)))
        ALLOCATE(weigC_grid(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2)))
        CALL Weigc_gr


    END SUBROUTINE set_var


    SUBROUTINE Weigc_gr

        IMPLICIT NONE

        INTEGER :: i,j,k
        !!!!Calcolo la giglia di weigC!!!Forse posso anche evitarla

        !ATTENZIONE
        !QUESTA PARTE  ANCORA VINCOLATA AD ESSERE 3D!!!!!!!!!!
        !Fai in modo che weigC sia come b e p
        DO k = indv(3,1),indv(3,2)
            DO j = indv(2,1),indv(2,2)
                DO i = indv(1,1),indv(1,2)
                    weigC_grid(i,j,k) = Lapl(1)%matrix(i,0)+ Lapl(2)%matrix(j,0)+ Lapl(3)%matrix(k,0)
                END Do
            END DO
        END DO

    END SUBROUTINE Weigc_gr

END MODULE set_pressure
