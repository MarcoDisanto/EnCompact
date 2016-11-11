PROGRAM codice

USE MPI
USE essentials 
USE MPI_module
USE input
USE grids
USE variables
USE compact
USE ic_and_bc 
USE set_pressure 
USE solve_pressure 
USE SPIKE

IMPLICIT NONE

INTEGER :: ierr, i
INTEGER, PARAMETER :: out_unit = 1
INTEGER :: MPI_v, MPI_sub_v, MPI_v_len
CHARACTER(len=140) :: MPI_v_str

CALL MPI_INIT(ierr)

! Current process reads data from file (nothing massive, just an handful
! of  scalars and strings).
CALL readin

! Creation of cartesian communicator and everything else is needed to identify
! the current process (i.e.: is it near a boundary? Which is its neighborhood? Which
! are its cartesian coordinates inside the cartesian communicator? ...)
CALL set_topology

! Determination of the number of cells of the current process along each
! dimension.
CALL distribute_cells

! The master process sets the grids (such grids are allocated only by the master
! process, since they are used only to build the matrices of the compact scheme.
! After this, they can also be deallocated).
CALL set_grids

! The master process builds the matrices of compact schemes, then it scatters them
! among other processes.
CALL set_compacts

! The current process allocates the structures aimed to contain velocity,
! presure and other variables.
CALL set_variables

! The current process sets an initial condition.
CALL set_ic

! The current process sets the boundary condition (if it touches any boundary).
CALL set_bc

! The current process inserts NaNs on values of velocity that should never be used.
CALL set_nans

! Initialize the SPIKE algorithm.
CALL SPIKE_init


! le seguenti assegnazioni devono provenire da una lettura da file
tol=1e-8
It_type='SRJ'
! It_type='J'
IF(myid==0)THEN

    ALLOCATE(pres_schemes(ndims,1:2))  

    Do i=1,ndims     

        ALLOCATE(pres_schemes(i,1)%sch(-1:1), pres_schemes(i,2)%sch(-1:1))
        pres_schemes(i,1)%sch(-1)='zero'
        pres_schemes(i,1)%sch( 0)='D1_c2f_2E_C'
        pres_schemes(i,1)%sch(+1)='zero'

        pres_schemes(i,2)%sch(-1)='D1_f2c_2E_C'
        pres_schemes(i,2)%sch( 0)='D1_f2c_2E_C'
        pres_schemes(i,2)%sch(+1)='D1_f2c_2E_C'

        !ALLOCATE(pres_schemes(i,1)%sch(-3:3),pres_schemes(i,2)%sch(-2:2))
        !pres_schemes(i,1)%sch(-3)='D1_c2f_2E_L'
        !pres_schemes(i,1)%sch(-2)='D1_c2f_3E_L'
        !pres_schemes(i,1)%sch(-1)='D1_c2f_4E_C'
        !pres_schemes(i,1)%sch( 0)='D1_c2f_6E_C'
        !pres_schemes(i,1)%sch(+1)='D1_c2f_4E_C'
        !pres_schemes(i,1)%sch(+2)='D1_c2f_3E_R'
        !pres_schemes(i,1)%sch(+3)='D1_c2f_2E_R'

        !pres_schemes(i,2)%sch(-2)='D1_f2c_3E_L'
        !pres_schemes(i,2)%sch(-1)='D1_f2c_4E_C'
        !pres_schemes(i,2)%sch( 0)='D1_f2c_6E_C'
        !pres_schemes(i,2)%sch(+1)='D1_f2c_4E_C'
        !pres_schemes(i,2)%sch(+2)='D1_f2c_3E_R'

        !pres_schemes(i,1)%sch(-3)='D1_c2f_6E_LLL'
        !pres_schemes(i,1)%sch(-2)='D1_c2f_6E_LL'
        !pres_schemes(i,1)%sch(-1)='D1_c2f_6E_L'
        !pres_schemes(i,1)%sch( 0)='D1_c2f_6E_C'
        !pres_schemes(i,1)%sch(+1)='D1_c2f_6E_R'
        !pres_schemes(i,1)%sch(+2)='D1_c2f_6E_RR'
        !pres_schemes(i,1)%sch(+3)='D1_c2f_6E_RRR'

        !pres_schemes(i,2)%sch(-2)='D1_f2c_6E_LL'
        !pres_schemes(i,2)%sch(-1)='D1_f2c_6E_L'
        !pres_schemes(i,2)%sch( 0)='D1_f2c_6E_C'
        !pres_schemes(i,2)%sch(+1)='D1_f2c_6E_R'
        !pres_schemes(i,2)%sch(+2)='D1_f2c_6E_RR'

    END DO


end if


CALL Set_p


! Allocation of source term in the Poisson equation
ALLOCATE(b(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2)))
! Questo termine deve essere in realtà la divergenza del campo *
b = .0
IF (myid == 13) THEN
    b(1,2,1) = 1.0
    b(1,3,1) = -1.0
END IF

! uvwp(4)%values=100.0
p = 123.0
! CALL RANDOM_NUMBER(p)
! p = myid

! Master process deallocates grids (which have been used to set the compacts),
! and other stuff.
CALL deallocate_grids

! Solve Poisson equation for pressure
!CALL Solve_p

! writing text files for the pressure field's post processing phase
CALL MPI_BARRIER(procs_grid, ierr)
DO i = 0, nprocs-1
    IF (myid == i) THEN
        IF (i == 0) THEN
            OPEN (UNIT = out_unit, FILE = 'domain.txt', ACTION = 'write', STATUS = 'replace')
        ELSE
            OPEN (UNIT = out_unit, FILE = 'domain.txt', ACTION = 'write', ACCESS = 'append', STATUS = 'old')
        END IF
        WRITE(out_unit,*) myid, N(1), N(2), N(3)
        CLOSE(out_unit)
    END IF 
    CALL MPI_BARRIER(procs_grid, ierr)
END DO

DO i = 0, nprocs-1
    IF (myid == i) THEN
        IF (i == 0) THEN
            OPEN (UNIT = out_unit, FILE = 'pressure.txt', ACTION = 'write', STATUS = 'replace')
        ELSE
            OPEN (UNIT = out_unit, FILE = 'pressure.txt', ACTION = 'write', ACCESS = 'append', STATUS = 'old')
        END IF
        WRITE(out_unit,*) p
        ! WRITE(out_unit,*) uvwp(4)%values
        CLOSE(out_unit)
    END IF 
    CALL MPI_BARRIER(procs_grid, ierr)
END DO




CALL SPIKE_exchange


! print MPI library version
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (myid == 0) THEN
!    PRINT *, ''
!    PRINT *, ''
!    PRINT *, ''
!    CALL MPI_GET_LIBRARY_VERSION(MPI_v_str, MPI_v_len, ierr)
!    PRINT *, 'MPI_GET_LIBRARY_VERSION: ', MPI_v_str
!    CALL MPI_GET_VERSION(MPI_v, MPI_sub_v, ierr)
!    PRINT *, 'MPI_GET_VERSION: ', MPI_v, '.',  MPI_sub_v
!END IF

CALL MPI_FINALIZE(ierr)

END PROGRAM codice
