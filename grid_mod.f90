MODULE grids
! This module contains variables and functions related to the independent
! physical variables, i.e. the spacial grids (such grids are used to build
! dimensional differential operators).
!
! Procedures:
!    SUBROUTINE set_grids
!    SUBROUTINE deallocate_grids

    IMPLICIT NONE

    REAL,    DIMENSION(:), ALLOCATABLE :: L        ! Length of the domain along the ndims directions
    INTEGER, DIMENSION(:), ALLOCATABLE :: spacings ! integers identifying the spacing along the ndims directions (uniform, Chebyshev, ...)

    TYPE grid1D
    ! This type has just one field in order to emulate MATLAB's cell arrays.
    ! two_grids(2)%g would be the same as MATLAB's two_grids{2}.
        REAL, DIMENSION(:), ALLOCATABLE :: g
    END TYPE grid1D

    TYPE grid3D
    ! A variable of this type is going to contain the ndims meshes along the
    ! ndims spatial directions.
    ! The two fields are needed to distinguish between faces and centers of cells.
        TYPE(grid1D), DIMENSION(:), ALLOCATABLE :: f
        TYPE(grid1D), DIMENSION(:), ALLOCATABLE :: c
    END TYPE grid3D

    TYPE(grid1D), DIMENSION(:, :), ALLOCATABLE :: x_glob ! vectors of global 1D internal coordinates for every variable;
    TYPE(grid1D), DIMENSION(:, :), ALLOCATABLE :: x_loc  ! vectors of local 1D internal coordinates for every variable;
                                                         ! NOTE: only process 0 allocates them
                                                         ! first index : component
                                                         ! second index: direction

    TYPE(grid3D) :: m3Dfc_tot ! stands for 3 meshes (or 1 or 2, depending on
                              ! ndims), for faces and centers of cells.
    ! At this point, m3Dfc_tot%f(2)%g would access to the y-coordinates of the
    ! faces (obviously those perpendicular to the y-direction).

    ! The variable all_grids is going to contain the meshes relative to
    ! 1st index: direction
    ! 2nd   "  : order of derivation
    ! 3rd   "  : staggering (1: c2c or c2f, 2: f2c or f2f)
    ! 4th   "  : sides (1: LHS, 2: RHS)
    TYPE(grid1D), DIMENSION(:,:,:,:), ALLOCATABLE :: all_grids

CONTAINS

    SUBROUTINE set_grids

        USE essentials, ONLY : log2int => logical2integer, linspace
        USE MPI_module, ONLY : myid, Ntot, ndims, periodic

        IMPLICIT NONE

        INTEGER         :: i ! generic index
        REAL, PARAMETER :: pi = ACOS(-1.0) ! pi greco

        master_process: IF (myid == 0) THEN

            ALLOCATE(all_grids(ndims,0:2,2,2))
            ! 1st index: direction
            ! 2nd index: order of derivation
            ! 3rd index: staggering (1: c2c or c2f, 2: f2c or f2f)
            ! 4th index: sides (1: LHS, 2: RHS)

            ALLOCATE(m3Dfc_tot%f(ndims))
            ALLOCATE(m3Dfc_tot%c(ndims))

            for_each_dimension: DO i = 1,ndims ! For each dimension (1, 2 or 3)...

                ALLOCATE(m3Dfc_tot%f(i)%g(Ntot(i)+1))
                ALLOCATE(m3Dfc_tot%c(i)%g(Ntot(i)))
                ! Uniform mesh from 0 to 1 (other spacings are based on this)
                m3Dfc_tot%f(i)%g = linspace(0.0,1.0,Ntot(i)+1)

                mesh_spacing: SELECT CASE (spacings(i))

                    CASE (1) ! uniform (from 0 to L)
                        m3Dfc_tot%f(i)%g = m3Dfc_tot%f(i)%g*L(i)

                    CASE (2) ! Chebyshev (from 0 to L)
                        m3Dfc_tot%f(i)%g = (1 + COS(pi*(1-m3Dfc_tot%f(i)%g)))*L(i)/2

                    CASE DEFAULT
                        ! TODO: implement error control
                        PRINT *, 'ERROR in grid_mod.f90: no such spacing is coded'

                END SELECT mesh_spacing

                ! NOTE: in the periodic case the last face is equal to the first
                ! one and should not be considered as an unknown; however the
                ! following data is used to calculate the coefficients of the
                ! schemes, thus it needs to be included in order to avoid strange
                ! results due to out-of-bound.

                ! The grid of the centers has to be built from the grid of the faces.
                m3Dfc_tot%c(i)%g = (m3Dfc_tot%f(i)%g(1:Ntot(i)) + m3Dfc_tot%f(i)%g(2:Ntot(i)+1))/2

                ! Grids of centers and faces are used to assemble the grids of
                ! right- and left-hand side of the various compact schemes used.
                ALLOCATE(all_grids(i,0,1,1)%g((2 - log2int(periodic(i))):Ntot(i)))
                all_grids(i,0,1,1)%g = m3Dfc_tot%f(i)%g((2 - log2int(periodic(i))):Ntot(i))
                ALLOCATE(all_grids(i,0,1,2)%g(log2int(periodic(i)):(Ntot(i)+1-log2int(periodic(i)))))
                IF (.NOT. periodic(i)) THEN
                    all_grids(i,0,1,2)%g = [m3Dfc_tot%f(i)%g(1),&
                                          & m3Dfc_tot%c(i)%g,&
                                          & m3Dfc_tot%f(i)%g(Ntot(i)+1)]
                ELSE
                    all_grids(i,0,1,2)%g = m3Dfc_tot%c(i)%g
                END IF

                ALLOCATE(all_grids(i,0,2,1)%g(1:Ntot(i)))
                all_grids(i,0,2,1)%g = m3Dfc_tot%c(i)%g
                ALLOCATE(all_grids(i,0,2,2)%g(1:Ntot(i)+1-log2int(periodic(i))))
                all_grids(i,0,2,2)%g = m3Dfc_tot%f(i)%g(1:Ntot(i)+1-log2int(periodic(i)))

                ! incognite anche le facce di bordo, poi scartate (questo è ANOMALO [cerca questa parola nei file per individuare le dipendenze])
                ALLOCATE(all_grids(i,1,1,1)%g(1:Ntot(i)+1-log2int(periodic(i))))
                all_grids(i,1,1,1)%g = m3Dfc_tot%f(i)%g(1:Ntot(i)+1-log2int(periodic(i)))
                ALLOCATE(all_grids(i,1,1,2)%g(log2int(periodic(i)):Ntot(i)+1-log2int(periodic(i))))
                IF (.NOT. periodic(i)) THEN
                    all_grids(i,1,1,2)%g = [m3Dfc_tot%f(i)%g(1),&
                                          & m3Dfc_tot%c(i)%g,&
                                          & m3Dfc_tot%f(i)%g(Ntot(i)+1)]
                ELSE
                    all_grids(i,1,1,2)%g = m3Dfc_tot%c(i)%g
                END IF

                ALLOCATE(all_grids(i,1,2,1)%g(1:Ntot(i)))
                all_grids(i,1,2,1)%g = m3Dfc_tot%c(i)%g
                ALLOCATE(all_grids(i,1,2,2)%g(1:Ntot(i)+1-log2int(periodic(i))))
                all_grids(i,1,2,2)%g = m3Dfc_tot%f(i)%g(1:Ntot(i)+1-log2int(periodic(i)))

                ALLOCATE(all_grids(i,2,1,1)%g(1:Ntot(i)))
                all_grids(i,2,1,1)%g = m3Dfc_tot%c(i)%g
                ALLOCATE(all_grids(i,2,1,2)%g(log2int(periodic(i)):Ntot(i)+1-log2int(periodic(i))))
                IF (.NOT. periodic(i)) THEN
                    all_grids(i,2,1,2)%g = [m3Dfc_tot%f(i)%g(1),&
                                          & m3Dfc_tot%c(i)%g,&
                                          & m3Dfc_tot%f(i)%g(Ntot(i)+1)]
                ELSE
                    all_grids(i,2,1,2)%g = m3Dfc_tot%c(i)%g
                END IF

                ALLOCATE(all_grids(i,2,2,1)%g(2-log2int(periodic(i)):Ntot(i)))
                all_grids(i,2,2,1)%g = m3Dfc_tot%f(i)%g(2-log2int(periodic(i)):Ntot(i))
                ALLOCATE(all_grids(i,2,2,2)%g(1:Ntot(i)+1-log2int(periodic(i))))
                all_grids(i,2,2,2)%g = m3Dfc_tot%f(i)%g(1:Ntot(i)+1-log2int(periodic(i)))

            END DO for_each_dimension

        END IF master_process

    END SUBROUTINE set_grids


    SUBROUTINE deallocate_grids
    ! This subroutine deallocates the allocatable variables of this module

        USE MPI_module, ONLY : myid

        IMPLICIT NONE

        IF (myid == 0) THEN
            DEALLOCATE(m3Dfc_tot%f)
            DEALLOCATE(m3Dfc_tot%c)
            !DEALLOCATE(L)    ! NOTE: mi serviva
            DEALLOCATE(spacings)
            DEALLOCATE(all_grids)
        END IF

    END SUBROUTINE deallocate_grids


    SUBROUTINE set_glob_coordinates
      ! This subroutine defines global internal coordinates for each component
      ! and broadcast the

      USE MPI, ONLY: MPI_DOUBLE_PRECISION, MPI_INTEGER
      USE MPI_module, ONLY: myid, ndims, periodic, slab_comm, mycoords
      USE essentials, ONLY: log2int=>logical2integer
      !USE variables, ONLY: block, uvwp

      INTEGER :: i, j
      INTEGER :: ierr
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: dmn ! this one is used to allocate x_glob%g

      ! Pencil master receive info on coordinate array size
      ALLOCATE(dmn(ndims, ndims))
      dmn = 0
      DO i = 1, ndims
        DO j = 1, ndims
          IF (myid==0) THEN
            IF (j == i) THEN
              dmn(i, j) = SIZE(m3Dfc_tot%f(j)%g)-1-log2int(.NOT.periodic(j))
            ELSE
              dmn(i, j) = SIZE(m3Dfc_tot%c(j)%g)
            END IF
          END IF
          CALL MPI_BCAST(dmn(i, j), 1, MPI_INTEGER, 0, slab_comm(j), ierr)
        END DO
      END DO


      ! allocate array of 1D vectors containing global internal coordinates
      ALLOCATE(x_glob(ndims, ndims))
      DO i = 1, ndims
        DO j = 1, ndims

          ! Pencil masters allocate x_glob%g
          IF (mycoords(j)==0) THEN
            ALLOCATE(x_glob(i, j)%g(dmn(i, j)))
          END IF

          ! Process 0 calculates coordinate arrays
          IF (myid==0) THEN
            IF (j == i) THEN
              x_glob(i, j)%g = m3Dfc_tot%f(j)%g(2-log2int(periodic(j)) : SIZE(m3Dfc_tot%f(j)%g)-1)
            ELSE
              x_glob(i, j)%g = m3Dfc_tot%c(j)%g
            END IF
          END IF

        END DO
      END DO

      ! Global coordinate x_glob(i, j) is broadcasted to pencil masters
      DO i = 1, ndims
        DO j = 1, ndims
          IF (mycoords(j)==0) THEN
            CALL MPI_BCAST(x_glob(i, j)%g, SIZE(x_glob(i, j)%g), MPI_DOUBLE_PRECISION, 0, slab_comm(j), ierr)
          END IF
        END DO
      END DO

      ! allocate array of 1D vectors containing local internal coordinates
      ALLOCATE(x_loc(ndims, ndims))

    END SUBROUTINE set_glob_coordinates



    SUBROUTINE set_loc_coordinates
      ! This subroutine would be more fit to lie in grids module. However, in order
      ! to avoid loops, it had to be placed here. It calculates local internal
      ! coordinates for each variable.

      USE MPI_module, ONLY: N, ndims, Neq, dims, mycoords, pencil_comm, myid, idm
      USE MPI, ONLY: MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_PROC_NULL
      USE essentials, ONLY: log2int=>logical2integer, KronDelta

      IMPLICIT NONE

      INTEGER :: ic, id, ip, ierr
      INTEGER, DIMENSION(:), ALLOCATABLE :: nd ! number of internal variables along a given direction
                                               ! lenght equal to n° of processes along direction
                                               ! only pencil masters allocate it
      INTEGER :: nc ! number of internal variables for a specific component
      INTEGER, DIMENSION(:), ALLOCATABLE :: displ ! displacements required by SCATTERV

      ! each local coordinates vector is allocated using bound taken from uvwp so
      ! as to avoid any recalculation
      DO ic = 1, ndims
        DO id = 1, ndims
          ALLOCATE(x_loc(ic, id)%g(1 + KronDelta(id,ic)*log2int(idm(id) == MPI_PROC_NULL) : N(id)))
        END DO
      END DO

      directions: DO id = 1, ndims

        ALLOCATE(nd(0 : dims(id)-1))
        ALLOCATE(displ(0 : dims(id)-1))
        nd = 0
        displ = 0

        components: DO ic = 1, ndims

          nc = Neq(id, 2, 1 + log2int(ic==id))
          CALL MPI_GATHER(nc, 1, MPI_INTEGER, nd, 1, MPI_INTEGER, 0, pencil_comm(id), ierr)

          IF (mycoords(id) == 0) THEN
            DO ip = 1, dims(id)-1
              displ(ip) = displ(ip-1) + nd(ip-1)
            END DO
          END IF
          CALL MPI_SCATTERV(x_glob(ic, id)%g, nd, displ, MPI_DOUBLE_PRECISION, &
                            x_loc(ic, id)%g,  nc, MPI_DOUBLE_PRECISION, 0, &
                            pencil_comm(id), ierr)


        END DO components

        DEALLOCATE(nd)
        DEALLOCATE(displ)
      END DO directions

    END SUBROUTINE set_loc_coordinates


END MODULE grids
