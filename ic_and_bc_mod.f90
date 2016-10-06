MODULE IC_and_BC 
! This module contains the following procedures
!    SUBROUTINE set_bc
!    SUBROUTINE assign_bc_i_3D(blocks, faceid)
!    SUBROUTINE set_face_condition(ndims_array, idf)
!    SUBROUTINE set_ic
        
    USE MPI, ONLY : MPI_COMM_WORLD

    IMPLICIT NONE

    ! A variable of the following derived data type is meant to contain
    ! the boundary conditions of velocity. It is allocated at read-time,
    ! based on ndims, by the subroutine readin in the input_mod.f90 file.
    TYPE boundary_conditions
        CHARACTER(LEN = 20), DIMENSION(:),   ALLOCATABLE :: string
        REAL,                DIMENSION(:,:), ALLOCATABLE :: values
    END TYPE boundary_conditions
    ! Commit boundary conditions array
    TYPE(boundary_conditions) :: uvw_BC

    ! Interface to put subroutines relative to 2D and 3D cases under the same name
    INTERFACE assign_bc_i
        ! MODULE PROCEDURE assign_bc_i_2D
        MODULE PROCEDURE assign_bc_i_3D
    END INTERFACE assign_bc_i

    ! Make public only a few stuff
    PRIVATE
    PUBLIC :: set_bc, set_ic, uvw_BC

CONTAINS

    SUBROUTINE set_bc
    ! This subroutine sets the boundary condition along the boundaries of the
    ! domain. It actually calls another subroutine to set the boundary condition
    ! along just one direction and one side.

        USE MPI,        ONLY : MPI_PROC_NULL
        USE MPI_module, ONLY : ndims, idm, idp
        USE variables,  ONLY : uvwp

        IMPLICIT NONE

        ! internal variables
        INTEGER :: id ! index of direction

        for_each_direction: DO id = 1, ndims

            if_touches_minus_side: IF (idm(id) == MPI_PROC_NULL) THEN

                ! Here the process touches the minus boundary along the
                ! id-th direction; assign the boundary condition, to the
                ! velocity (the 1st, 2nd and 3rd components of uvwp).
                CALL assign_bc_i(uvwp(1:ndims), -id)
                
            END IF if_touches_minus_side

            if_touches_plus_side: IF (idp(id) == MPI_PROC_NULL) THEN

                ! Here the process touches the plus boundary along the
                ! id-th direction; assign the boundary condition, ...
                CALL assign_bc_i(uvwp(1:ndims), +id)

            END IF if_touches_plus_side

        END DO for_each_direction

    END SUBROUTINE set_bc


    ! TO DO: aggiungere la subroutine assign_bc_i_2D per il caso 2D. Si tratta
    ! di un copia e incolla della assign_bc_i_3D con modifiche abbastanza immediate.
    SUBROUTINE assign_bc_i_3D(blocks, faceid)
    ! This is one of two subroutines (the other one is assign_bc_i_2D) which are
    ! called under the same name, assign_bc_i. The present should be called when
    ! ndims == 3, the other one when ndims == 2.
    ! Both subroutines assign the boundary condition at the boundary (second
    ! argument), of a domain (first argument).

        USE MPI_module, ONLY : ndims
        USE variables,  ONLY : block, flat

        IMPLICIT NONE

        ! IN/OUT/INOUT variables
        TYPE(block), DIMENSION(ndims), INTENT(INOUT), TARGET :: blocks
        INTEGER,                       INTENT(IN)            :: faceid
        ! faceid can be    -3 -2 -1 +1 +2 +3
        ! corresponding to  B  D  L  R  U  F

        ! internal variables
        TYPE(flat), DIMENSION(ndims) :: flats
        INTEGER :: ierr, errorcode
        INTEGER :: ic

        ! The following case construct makes 'flats' elements point
        ! to the 'faceid' face of the 'blocks' array elements
        select_face: SELECT CASE (faceid)

        CASE (-1) ! Left face
            DO ic = 1, ndims
                flats(ic)%values => blocks(ic)%values(blocks(ic)%b_bc(1,1), &
                                                    & blocks(ic)%b_bc(2,1):blocks(ic)%b_bc(2,2), &
                                                    & blocks(ic)%b_bc(3,1):blocks(ic)%b_bc(3,2))
            END DO

        CASE (+1) ! Right face
            DO ic = 1, ndims
                flats(ic)%values => blocks(ic)%values(                     blocks(ic)%b_bc(1,2), &
                                                    & blocks(ic)%b_bc(2,1):blocks(ic)%b_bc(2,2), &
                                                    & blocks(ic)%b_bc(3,1):blocks(ic)%b_bc(3,2))
            END DO

        CASE (-2) ! Down face
            DO ic = 1, ndims
                flats(ic)%values => blocks(ic)%values(blocks(ic)%b_bc(1,1):blocks(ic)%b_bc(1,2), & 
                                                    & blocks(ic)%b_bc(2,1), &
                                                    & blocks(ic)%b_bc(3,1):blocks(ic)%b_bc(3,2))
            END DO

        CASE (+2) ! Up face
            DO ic = 1, ndims
                flats(ic)%values => blocks(ic)%values(blocks(ic)%b_bc(1,1):blocks(ic)%b_bc(1,2), & 
                                                    &                      blocks(ic)%b_bc(2,2), &
                                                    & blocks(ic)%b_bc(3,1):blocks(ic)%b_bc(3,2))
            END DO

        CASE (-3) ! Back face
            DO ic = 1, ndims
                flats(ic)%values => blocks(ic)%values(blocks(ic)%b_bc(1,1):blocks(ic)%b_bc(1,2), & 
                                                    & blocks(ic)%b_bc(2,1):blocks(ic)%b_bc(2,2), &
                                                    & blocks(ic)%b_bc(3,1))
            END DO

        CASE (+3) ! Front face
            DO ic = 1, ndims
                flats(ic)%values => blocks(ic)%values(blocks(ic)%b_bc(1,1):blocks(ic)%b_bc(1,2), & 
                                                    & blocks(ic)%b_bc(2,1):blocks(ic)%b_bc(2,2), &
                                                    &                      blocks(ic)%b_bc(3,2))
            END DO

        CASE DEFAULT
            PRINT *, 'ERROR in assign_bc_i_3D (ic_and_bc_mod.f90): How is it possible &
                    & to end up here?!'
            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            STOP

        END SELECT select_face

        CALL set_face_condition(flats, faceid)
        
    END SUBROUTINE assign_bc_i_3D


    SUBROUTINE set_face_condition(ndims_array, idf)
    ! This subroutine assigns the values to 'array', on the face 'idf',
    ! based on the array of strings uvw_BC, which is in input from file
    ! and is indicative of the boundary condition.

        USE MPI_module, ONLY : ndims
        USE variables,  ONLY : flat
        USE essentials

        IMPLICIT NONE

        ! IN/OUT/INOUT variables
        TYPE(flat), DIMENSION(ndims), INTENT(INOUT) :: ndims_array
        INTEGER,                      INTENT(IN)    :: idf

        ! internal varibles
        INTEGER :: ic

        for_each_component: DO ic = 1, ndims
            ndims_array(ic)%values = uvw_BC%values(idf, ic)
        END DO for_each_component

    END SUBROUTINE set_face_condition


    SUBROUTINE set_ic
    ! This subroutine sets the initial condition in terms of initial velocity
    ! (the initial pressure is non-sense in the context of a pressure
    ! projection method).
    ! TO DO: in a future, I could implement a conditional to retrieve the
    ! initial conditions from file.

        USE MPI,       ONLY : MPI_PROC_NULL
        USE variables, ONLY : uvwp, n_flow_variables
        USE MPI_module, ONLY: myid, N, idm, idp
        USE essentials, ONLY : log2int => logical2integer, KronDelta
        USE, INTRINSIC :: IEEE_ARITHMETIC ! to use IEEE routines
        ! such as IEEE_QUIET_NAN, IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF
    
        IMPLICIT NONE

        ! internal variables
        INTEGER :: iv ! generic index

        ! IEEE
        REAL :: r = 0.0, NaN ! r is a dummy real used to define a NaN of type real
        NaN = IEEE_VALUE(r, IEEE_QUIET_NAN) ! NaN of the same type as r

        zero_ic: DO iv = 1, n_flow_variables-1

            uvwp(iv)%values = NaN
            uvwp(iv)%values(uvwp(iv)%b(1,1):uvwp(iv)%b(1,2), &
                          & uvwp(iv)%b(2,1):uvwp(iv)%b(2,2), &
                          & uvwp(iv)%b(3,1):uvwp(iv)%b(3,2)) = myid

        END DO zero_ic

    END SUBROUTINE set_ic

END MODULE IC_and_BC
