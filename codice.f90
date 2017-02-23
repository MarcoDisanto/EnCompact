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
USE diffusive_term
USE convective_term
USE bandedmatrix
USE time_advancement
USE output
USE Thomas_suite
USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines


IMPLICIT NONE

INTEGER :: ierr, ic, id, istag
INTEGER, PARAMETER :: out_unit = 1
INTEGER :: MPI_v, MPI_sub_v, MPI_v_len
CHARACTER(len=140) :: MPI_v_str
REAL, DIMENSION(3, 3, 3) :: A
INTEGER :: i, j, k, iii
REAL, DIMENSION(:), ALLOCATABLE :: q, psi
REAL, DIMENSION(:), ALLOCATABLE :: coeffs
INTEGER, POINTER :: lid, uid
REAL :: dt, ni, vel_tol
INTEGER :: sz, prova, nun
INTEGER(MPI_ADDRESS_KIND) :: i1, i2, diff
INTEGER :: ext1, ext2
CHARACTER (LEN = 1024) :: strtry
INTEGER :: wrtind, wrtfreq         ! required to output multiple times
REAL               :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r


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
        !IF (periodic(i)) THEN
        !   pres_schemes(i,1)%sch(-1)='periodic'
        !   pres_schemes(i,1)%sch( 0)='D1_c2f_2E_C'
        !   pres_schemes(i,1)%sch(+1)='periodic'

        !   pres_schemes(i,2)%sch(-1)='D1_f2c_2E_C'
        !   pres_schemes(i,2)%sch( 0)='D1_f2c_2E_C'
        !   pres_schemes(i,2)%sch(+1)='D1_f2c_2E_C'
        !ELSE
           pres_schemes(i,1)%sch(-1)='zero'
           pres_schemes(i,1)%sch( 0)='D1_c2f_2E_C'
           pres_schemes(i,1)%sch(+1)='zero'

           pres_schemes(i,2)%sch(-1)='D1_f2c_2E_C'
           pres_schemes(i,2)%sch( 0)='D1_f2c_2E_C'
           pres_schemes(i,2)%sch(+1)='D1_f2c_2E_C'
        !END IF

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

! Master process deallocates grids (which have been used to set the compacts),
! and other stuff.
CALL deallocate_grids

CALL set_diff_bounds
CALL set_conv_bounds
CALL set_conv_bc
CALL set_grad_p
CALL set_output
CALL set_tot_der ! call this subroutine only if multistage methods are used

! Allocation of source term in the Poisson equation
ALLOCATE(b(indv(1,1):indv(1,2),indv(2,1):indv(2,2),indv(3,1):indv(3,2)))
! Questo termine deve essere in realtà la divergenza del campo *

! Decide whether to use pentdag or cypent - basically a periodic pentdag - depending
! on the values of variable "periodic"
DO id = 1, ndims

  IF (periodic(id)) THEN
    pent_point(id)%true_point => cypent
  ELSE
    pent_point(id)%true_point => pentdag
  END IF

END DO

!IF (myid==0) PRINT *, 'n° of processes'
!IF (myid==0) PRINT *, dims

! TODO: inserire dt e ni nel file di input e modificare il modulo di lettura in
! maniera opportuna
dt = 0.0104 ! C = 1
ni = 1e-3
vel_tol = 1e-3 ! criterio di arresto


IF (myid==0) PRINT *, ''
IF (myid==0) PRINT *, ''
IF (myid==0) PRINT *, ''
IF (myid==0) PRINT *, 'Start calculating'

!!!!!!!!!! Time advancement !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i = 0
vel_res = 1
wrtfreq = 500 ! how often do you want to write an output file
wrtind = 1
CALL SPIKE_exchange_uvw ! exchanging velocity values
t_start = MPI_Wtime()

!DO WHILE (vel_res>vel_tol .AND. i<100000)
!
!  i = i+1
!
!  ! upgrading solution
!  CALL RK4(dt, ni)
!  !CALL ExplEuler(dt, ni)
!
!  CALL residual_eval(dt)
!
!  IF (myid==0) PRINT *, i, '     Residual = ', vel_res
!
!  ! Output
!  !IF (MOD(i, wrtfreq)==0) THEN
!  !  CALL raw_out(wrtind)
!  !  IF (myid==0) PRINT *, 'Output written'
!  !  wrtind = wrtind+1
!  !END IF
!
!END DO
t_end = MPI_Wtime()
CALL divergence_calc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculating execution time
t_exec_proc = t_end-t_start
CALL MPI_REDUCE(t_exec_proc, t_exec, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

! Output
!CALL set_output
wrtind = 1
!CALL raw_out(wrtind)
!CALL output_read_aid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
id = 1
istag = 2

iii = 0
i = 1

!PRINT *, 'Process', myid, '     Velocity difference =', vel_diff_proc


IF (myid==0) PRINT *, ''
IF (myid==0) PRINT *, ''
IF (myid==0) PRINT *, ''
IF (myid==0) PRINT *, 'Elapsed time =', t_exec
!IF (myid==iii) CALL printmatrix(uvwp(i)%values)

!IF (myid==iii) print *, 'grad'
!IF (myid==iii) CALL printmatrix(GraDiv(i, 1)%matrix)
!IF (myid==iii) print *, 'div'
!IF (myid==iii) CALL printmatrix(GraDiv(i, 2)%matrix)
!IF (myid==iii) print *, 'lapl'
!IF (myid==iii) CALL printmatrix(Lapl(i)%matrix)



! print MPI library version
CALL MPI_BARRIER(procs_grid, ierr)
IF (myid == 0) THEN
    PRINT *, ''
    PRINT *, ''
    PRINT *, ''
    CALL MPI_GET_LIBRARY_VERSION(MPI_v_str, MPI_v_len, ierr)
    PRINT *, 'MPI_GET_LIBRARY_VERSION: ', MPI_v_str
    CALL MPI_GET_VERSION(MPI_v, MPI_sub_v, ierr)
    PRINT *, 'MPI_GET_VERSION: ', MPI_v, '.',  MPI_sub_v
END IF

CALL MPI_FINALIZE(ierr)

END PROGRAM codice
