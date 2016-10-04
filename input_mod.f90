MODULE input
    
    USE MPI,          ONLY : MPI_COMM_WORLD
    USE MPI_module,   ONLY : ndims, Ntot, periodic, reorder, dims, myid
    USE grids,        ONLY : L, spacings
    USE IC_and_BC,    ONLY : uvw_BC
    USE compact,      ONLY : all_schemes
    USE set_pressure, ONLY : pres_schemes

    IMPLICIT NONE

CONTAINS

    SUBROUTINE readin
    
    IMPLICIT NONE
    
    ! internal variables
    INTEGER :: ierr, errorcode ! error of MPI routines
    INTEGER :: ii, jj, i, j, k ! generic indices
    INTEGER :: stat ! status of read statement
    INTEGER :: f_in = 10 ! input file unit number
    INTEGER, PARAMETER :: Nv = 9 ! Number  of variables
    INTEGER            :: Nmv    ! Number  of mandatory variables (dynamically changed based on ndims)
    INTEGER            :: Cmv    ! Counter of mandatory variables
    CHARACTER (LEN=20)                :: Vname  ! Variable's name
    CHARACTER (LEN=20), DIMENSION(Nv) :: Vnames ! Array of variables' names
    REAL                                          :: rscalar      ! real     value   to be read
    REAL,               DIMENSION(:), ALLOCATABLE :: rarray       ! array of values  to be read
    INTEGER                                       :: iscalar      ! integer  value   to be read
    INTEGER,            DIMENSION(:), ALLOCATABLE :: iarray       ! array of values  to be read
    CHARACTER (LEN=20)                            :: str1, str2   ! dummy    strings
    CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE :: strings      ! array of strings to be read
    CHARACTER (LEN=20)                            :: dir_der_stag ! scalar   str1 containing direction, order of derivation, and staggering
    INTEGER :: idf
    INTEGER, DIMENSION(2) :: tesa ! 2-elements support array

    ! open the stream
    OPEN(f_in, FILE = 'IN', ACCESS = 'sequential', FORM = 'formatted', ACTION = 'read')

    ! set default values of optional variables
    reorder = .TRUE.
    Nmv = 1 ! it's changed based on ndims

    ! here set the names of variables to be read in IN file
    Vnames( 1) = 'ndims'   ! number of dimensions
    Vnames( 2) = 'N'       ! number of total cells per direction
    Vnames( 3) = 'lengths' ! lenght of the domain along each direction
    Vnames( 4) = 'periodic'! periodicity along each direction
    Vnames( 5) = 'spacings'! spacing along each direction
    Vnames( 6) = 'np'      ! number of processes along each direction
    Vnames( 7) = 'reorder' ! reorder the processes to fit cpu architecture
    Vnames( 8) = 'scheme'  ! differential schemes
    Vnames( 9) = 'uvw_BC'      ! boundary conditions
    

    Cmv = 0 ! counter of mandatory variables initialized to 0
    
    ReadFile: DO
    
        SkipComments: DO
            READ (f_in, *, IOSTAT=stat) Vname
            IF (stat.LT.0)          EXIT ReadFile
            IF (Vname (1:1) /= "#") EXIT SkipComments
        END DO SkipComments
    
        BACKSPACE(f_in)
        READ(f_in,*) Vname ! read the name of variable only
    
        ii = 0
    
        DO i = 1, Nv
            IF(Vname .EQ. Vnames(i)) THEN
                ii = i
            END IF
        END DO
    
        SELECT CASE (ii)
    
            CASE(1)
                BACKSPACE(f_in)
                READ(f_in,*) Vname, iscalar
                ndims = iscalar
                ! allocate arrays of dimension ndims
                ALLOCATE(L(ndims), periodic(ndims), Ntot(ndims), spacings(ndims), dims(ndims))
                ALLOCATE(all_schemes(ndims, 0:2, 1:2))
                ! increase the number of mandatory variables to be read
                Nmv = Nmv + 5*ndims
                ALLOCATE(uvw_BC%string(-ndims:+ndims))
                ALLOCATE(uvw_BC%values(-ndims:+ndims, ndims))
                ! TO DO: qua dovrei mettere le assegnazioni di default
                DO i = -ndims, +ndims
                    uvw_BC%string(i) = 'none'
                END DO
                uvw_BC%values = 0
                Nmv = Nmv + 2*ndims
                Cmv = Cmv + 1
    
            CASE(2)
                BACKSPACE(f_in)
                ALLOCATE(iarray(ndims))
                READ(f_in,*) Vname, iarray
                Ntot = iarray
                DEALLOCATE(iarray)
                Cmv  = Cmv + SIZE(Ntot)
    
            CASE(3)
                BACKSPACE(f_in)
                ALLOCATE(rarray(ndims))
                READ(f_in,*) Vname, rarray
                L   = rarray
                DEALLOCATE(rarray)
                Cmv = Cmv + SIZE(L)
    
            CASE(4)
                BACKSPACE(f_in)
                ALLOCATE(rarray(ndims))
                READ(f_in,*) Vname, rarray
                DO i = 1,ndims
                   IF      (rarray(i) == 1.0) THEN
                       periodic(i) = .TRUE.
                   ELSE IF (rarray(i) == 0.0) THEN
                       periodic(i) = .FALSE.
                   END IF
                END DO
                DEALLOCATE(rarray)
                Cmv = Cmv + SIZE(periodic)
    
            CASE(5)
                BACKSPACE(f_in)
                ALLOCATE(iarray(ndims))
                READ(f_in,*) Vname, iarray
                spacings = iarray
                DEALLOCATE(iarray)
                Cmv      = Cmv + SIZE(spacings)
    
            CASE(6)
                BACKSPACE(f_in)
                ALLOCATE(iarray(ndims))
                READ(f_in,*) Vname, iarray
                dims = iarray
                DEALLOCATE(iarray)
                Cmv  = Cmv + SIZE(dims)
    
            CASE(7)
                BACKSPACE(f_in)
                READ(f_in,*) Vname, rscalar
                IF (rscalar == 0.0) reorder = .FALSE.
    
            CASE(8)
                ! Only the master process reads the schemes
                master_process_read_schemes: IF (myid == 0) THEN
                    BACKSPACE(f_in)
                    READ(f_in,*) Vname, dir_der_stag, iscalar
                    ! **************************************************************
                    ! ********** TO DO: questa limitazione dovrebbe essere eliminata
                    IF (MOD(iscalar,2) == 0) THEN
                        WRITE(*,*) 'ERROR (in input_mod.f90). An even number of schemes along&
                                    & a direction is not allowed (not yet).'
                        CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                        STOP
                    END IF
                    ! ***************** questa limitazione dovrebbe essere eliminata
                    ! **************************************************************
                    ALLOCATE(strings((1-iscalar)/2:(iscalar-1)/2))
                    BACKSPACE(f_in)
                    READ(f_in,*) Vname, dir_der_stag, iscalar, str1, strings
                    
                    ! determine the direction
                    READ(dir_der_stag(4:4),*) i
                    ! determine the order of derivation
                    READ(dir_der_stag(9:9),*) j
                    ! determine the staggering
                    IF      ((dir_der_stag(11:13) .EQ. 'c2f') .OR. (dir_der_stag(11:13) .EQ. 'c2c')) THEN
                        k = 1
                    ELSE IF ((dir_der_stag(11:13) .EQ. 'f2c') .OR. (dir_der_stag(11:13) .EQ. 'f2f')) THEN
                        k = 2
                    ELSE
                        WRITE(*,*) 'ERROR (in input_mod.f90). Unrecognized &
                                    & staggering of FD schemes in input file:', &
                                    & dir_der_stag(11:13)
                        CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                        STOP
                    END IF
                    ALLOCATE(all_schemes(i,j,k)%sch((1-iscalar)/2:(iscalar-1)/2))
                    DO jj = (1-iscalar)/2,(iscalar-1)/2
                        all_schemes(i,j,k)%sch(jj) = 'D'//dir_der_stag(9:13)//'_'//strings(jj)
                    END DO
                    DEALLOCATE(strings)

                END IF master_process_read_schemes

            CASE(9)

                BACKSPACE(f_in)
                READ(f_in,*) Vname, str1, str2

                read_face: SELECT CASE (str1) ! idf is the face index, and it is minus the direction

                CASE ('Left')
                    idf = -1
                    
                CASE ('Right')
                    idf = +1
                    
                CASE ('Down')
                    idf = -2
                    
                CASE ('Up')
                    idf = +2
                    
                CASE ('Back')
                    idf = -3
                    
                CASE ('Front')
                    idf = +3

                CASE DEFAULT
                    PRINT *, 'ERROR in input_mod.f90: unrecognized face in &
                            & read_face CASE construct.'
                    CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                    STOP
                    
                END SELECT read_face

                ! assign the string identifying the BC
                uvw_BC%string(idf) = str2

                read_appended_values: SELECT CASE (str2) 

                CASE ('steady_solid_wall')
                    ! no needed values for steady wall
                    uvw_BC%string(idf) = str2
                    uvw_BC%string(0)   = 'do_not_use'
                    uvw_BC%values(idf,:) = 0 ! zero velocity along all directions

                CASE ('sliding_solid_wall')
                    ! need to read tangetial components of velocity
                    ! the tangential components are 2 or 1 in 3D and 2D respectively
                    ALLOCATE(rarray(ndims-1)) 
                    ! read again
                    BACKSPACE(f_in)
                    READ(f_in,*) Vname, str1, str2, rarray

                    uvw_BC%string(0)   = 'do_not_use' ! it's just a remainder

                    uvw_BC%values(idf, ABS(idf)) = 0 ! set to zero the normal component

                    ! the following tesa contains the array [1, ..., ndims] without idf
                    ! (so it can be [1, 2], [1, 3], or [2, 3] in 3D, and signletons
                    ! [1] or [2] in 2D), that is, the two tangential directions.
                    tesa = [(jj, jj = 1, ABS(idf)-1), (jj, jj = ABS(idf)+1, ndims)]
                    DO jj = 1, ndims-1
                        ! insert read values as tangential components of velocity
                        uvw_BC%values(idf, tesa(jj)) = rarray(jj)
                    END DO

                    DEALLOCATE(rarray) 

                CASE DEFAULT
                    PRINT *, 'Error in input_mod.f90: boundary condition not recognized.'
                    CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                    STOP

                END SELECT read_appended_values
                Cmv  = Cmv + 1
    
            CASE DEFAULT
                WRITE(*,*) 'ERROR (in input_mod.f90). Field not identified in input: ', Vname
                CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
                STOP
    
        END SELECT
    
    END DO ReadFile
    
    !1 CLOSE(f_in)
    
    IF (Cmv < Nmv) THEN
        WRITE(*,*) 'ERROR: ', Nmv - Cmv, &
                 & ' mandatory input variable(s) missing (in rdin.f90).'
        CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
        STOP
    END IF
    
    RETURN
    
    END SUBROUTINE readin

END MODULE input
