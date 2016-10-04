MODULE library

IMPLICIT NONE

CONTAINS

    FUNCTION detBband(schemes)
    ! This routine determines the lower band, detBband(1), and upper band, detBband(2),
    ! of the matrix B, based on the variables schemes, which is supposed to contain the
    ! strings relative to the boundary and central schemes.
    ! NB: detBband(1) is minus the number of lower diagonals, whereas detBband(2) is the
    ! number of upper diagonals.

        USE mpi, ONLY : MPI_COMM_WORLD

        IMPLICIT NONE

        ! input variable
        CHARACTER(len = 20), DIMENSION (:), INTENT(IN) :: schemes

        ! output variable
        INTEGER,             DIMENSION(2)              :: detBband

        ! internal variable
        INTEGER :: i, ierr, errorcode

        detBband = 0
        for_each_scheme: DO i = lbound(schemes,1), ubound(schemes,1)

            scheme: SELECT CASE (schemes(i))

                CASE ('D0_c2f_6C_C')
                    ! 6-th order compact interpolation scheme from centers to faces
                    ! The array of the centers is supposed to be completed with the
                    ! first and last faces.
                    ! The array of faces doesn't contain the first and last faces.
                    !   o   o       o       o                 o       o   o
                    !   |___|_______|_______|__  .....  ______|_______|___|
                    !           |       |                 |       |
                    !           o       o                 o       o
                    detBband = [min(detBband(1), -2), max(detBband(2), +1)]

                CASE ('D0_f2c_6C_C')
                    ! 6-th order compact interpolation scheme from faces to centers
                    !   o       o       o                o       o       o
                    !   |_______|_______|_  .....  ______|_______|_______|
                    !       |       |                |       |       |
                    !       o       o                o       o       o
                    detBband = [min(detBband(1), -1), max(detBband(2), +2)]

                CASE ('D0_c2f_4C_L')
                    !   o   o       o
                    !   |___|_______|____
                    !           |       |
                    !           o       o
                    detBband = [min(detBband(1), -2), max(detBband(2), 0)]

                CASE ('D0_f2c_5C_L')
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !       |       |
                    !       o       o
                    detBband = [min(detBband(1), 0), max(detBband(2), +3)]

                CASE ('D0_c2f_4C_R')
                    !       o       o   o
                    !   ____|_______|___|
                    !   |       |
                    !   o       o
                    detBband = [min(detBband(1), -1), max(detBband(2), +1)]

                CASE ('D0_f2c_5C_R')
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !               |       |
                    !               o       o
                    detBband = [min(detBband(1), -2), max(detBband(2), +1)]

                CASE ('D1_f2c_2E_C')
                    !              o       o
                    !              |_______|
                    !                  |
                    !                  \

                    detBband = [min(detBband(1), 0), max(detBband(2), +1)]

                CASE ('D1_f2c_4E_C')
                    !              o       o       o       o
                    !              |_______|_______|_______|
                    !                          |
                    !                          \

                    detBband = [min(detBband(1), -1), max(detBband(2), +2)]

                CASE ('D1_f2c_6E_C')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                          |
                    !                          \

                    detBband = [min(detBband(1), -2), max(detBband(2), +3)]

                CASE ('D1_f2c_4E_L')
                    ! 4-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o
                    !      |_______|_______|_______|_______|
                    !          |
                    !          \

                    detBband = [min(detBband(1), 0), max(detBband(2), +4)]

                CASE ('D1_f2c_5E_L')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                  |
                    !                  \

                    detBband = [min(detBband(1), -1), max(detBband(2), +4)]

                CASE ('D1_f2c_5E_LL')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !          |
                    !          \

                    detBband = [min(detBband(1), 0), max(detBband(2), +5)]

                CASE ('D1_f2c_6E_L')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                  |
                    !                  \

                    detBband = [min(detBband(1), -1), max(detBband(2), +5)]

                CASE ('D1_f2c_6E_LL')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !          |
                    !          \

                    detBband = [min(detBband(1), 0), max(detBband(2), +6)]

                CASE ('D1_f2c_4E_R')
                    ! 4-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o
                    !      |_______|_______|_______|_______|
                    !                                  |
                    !                                  \

                    detBband = [min(detBband(1), -3), max(detBband(2), +1)]

                CASE ('D1_f2c_5E_R')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                                  |
                    !                                  \

                    detBband = [min(detBband(1), -3), max(detBband(2), +2)]

                CASE ('D1_f2c_5E_RR')
                    ! 5-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|
                    !                                          |
                    !                                          \

                    detBband = [min(detBband(1), -4), max(detBband(2), +1)]

                CASE ('D1_f2c_6E_R')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                                          |
                    !                                          \

                    detBband = [min(detBband(1), -4), max(detBband(2), +2)]

                CASE ('D1_f2c_6E_RR')
                    ! 6-th order explicit first derivative scheme from faces to centers
                    !      o       o       o       o       o       o       o
                    !      |_______|_______|_______|_______|_______|_______|
                    !                                                  |
                    !                                                  /

                    detBband = [min(detBband(1), -5), max(detBband(2), 1)]

                CASE ('D1_c2f_2E_C')
                    !              o       o
                    !              |_______|
                    !                  |
                    !                  \

                    detBband = [min(detBband(1), -1), max(detBband(2), 0)]

                CASE ('D1_c2f_4E_C')
                    !              o       o       o       o
                    !              |_______|_______|_______|
                    !                          |
                    !                          \

                    detBband = [min(detBband(1), -2), max(detBband(2), +1)]

                CASE ('D1_c2f_6E_C')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                               |
                    !                               \

                    detBband = [min(detBband(1), -3), max(detBband(2), +2)]
                
                CASE ('D1_c2f_1E_L')
                    ! 1-th order explicit first derivative scheme from faces to centers
                    !           o       o
                    !        ___|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +1)]
                
                CASE ('D1_c2f_2E_L')
                    ! 2-th order explicit first derivative scheme from faces to centers
                    !           o       o       o
                    !        ___|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +2)]
                
                CASE ('D1_c2f_3E_L')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), -1), max(detBband(2), +2)]

                CASE ('D1_c2f_3E_LL')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o
                    !        ___|_______|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +3)]

                CASE ('D1_c2f_4E_L')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !           |_______|_______|_______|_______|
                    !               |
                    !               \

                    detBband = [min(detBband(1), -1), max(detBband(2), +3)]

                CASE ('D1_c2f_4E_LL')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !        ___|_______|_______|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +4)]

                CASE ('D1_c2f_5E_L')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                       |
                    !                       \

                    detBband = [min(detBband(1), -2), max(detBband(2), +3)]

                CASE ('D1_c2f_5E_LL')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !               |
                    !               \

                    detBband = [min(detBband(1), -1), max(detBband(2), +4)]

                CASE ('D1_c2f_5E_LLL')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !        ___|_______|_______|_______|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +5)]

                CASE ('D1_c2f_6E_L')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|_______|
                    !                       |
                    !                       \

                    detBband = [min(detBband(1), -2), max(detBband(2), +4)]

                CASE ('D1_c2f_6E_LL')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|_______|
                    !               |
                    !               \

                    detBband = [min(detBband(1), -1), max(detBband(2), +5)]

                CASE ('D1_c2f_6E_LLL')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o       o
                    !        ___|_______|_______|_______|_______|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +6)]
                
                CASE ('D1_c2f_1E_R')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !                   o       o
                    !                   |_______|___
                    !                               |
                    !                               \

                    detBband = [min(detBband(1), -2), max(detBband(2), -1)]
                
                CASE ('D1_c2f_2E_R')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !           o       o       o
                    !           |_______|_______|___
                    !                               |
                    !                               \

                    detBband = [min(detBband(1), -3), max(detBband(2), -1)]
                
                CASE ('D1_c2f_3E_R')
                    ! 3-th order explicit first derivative scheme from faces to centers
                    !                   o       o       o       o
                    !                   |_______|_______|_______|
                    !                                       |
                    !                                       \

                    detBband = [min(detBband(1), -3), max(detBband(2), 0)]

                CASE ('D1_c2f_3E_RR')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !                   o       o       o       o
                    !                   |_______|_______|_______|___
                    !                                               |
                    !                                               \

                    detBband = [min(detBband(1), -4), max(detBband(2), -1)]

                CASE ('D1_c2f_4E_R')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !           |_______|_______|_______|_______|
                    !                                       |
                    !                                       \

                    detBband = [min(detBband(1), -4), max(detBband(2), 0)]

                CASE ('D1_c2f_4E_RR')
                    ! 4-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o
                    !           |_______|_______|_______|_______|___
                    !                                               |
                    !                                               \

                    detBband = [min(detBband(1), -5), max(detBband(2), -1)]

                CASE ('D1_c2f_5E_R')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                                       |
                    !                                       \

                    detBband = [min(detBband(1), -4), max(detBband(2), +1)]

                CASE ('D1_c2f_5E_RR')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|
                    !                                               |
                    !                                               \

                    detBband = [min(detBband(1), -5), max(detBband(2), 0)]

                CASE ('D1_c2f_5E_RRR')
                    ! 5-th order explicit first derivative scheme from centers to faces
                    !           o       o       o       o       o       o
                    !           |_______|_______|_______|_______|_______|___
                    !                                                       |
                    !                                                       \

                    detBband = [min(detBband(1), -6), max(detBband(2), -1)]

                CASE ('D1_c2f_6E_R')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !    o       o       o       o       o       o       o
                    !    |_______|_______|_______|_______|_______|_______|
                    !                                        |
                    !                                        \

                    detBband = [min(detBband(1), -5), max(detBband(2), +1)]

                CASE ('D1_c2f_6E_RR')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !    o       o       o       o       o       o       o
                    !    |_______|_______|_______|_______|_______|_______|
                    !                                                |
                    !                                                \

                    detBband = [min(detBband(1), -6), max(detBband(2), +0)]

                CASE ('D1_c2f_6E_RRR')
                    ! 6-th order explicit first derivative scheme from centers to faces
                    !    o       o       o       o       o       o       o
                    !    |_______|_______|_______|_______|_______|_______|___
                    !                                                        |
                    !                                                        \

                    detBband = [min(detBband(1), -7), max(detBband(2), -1)]

                CASE ('D1_f2c_6C_C')
                    ! 6-th order compact first derivative scheme from faces to centers
                    !   o       o       o                o       o       o
                    !   |_______|_______|_  .....  ______|_______|_______|
                    !       |       |                |       |       |
                    !       \       \                \       \       \

                    detBband = [min(detBband(1), -1), max(detBband(2), +2)]

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

                    detBband = [min(detBband(1), -2), max(detBband(2), +1)]

                CASE ('D2_c2c_6C_C','D2_f2f_6C_C')
                    ! 6-th order compact second derivative collocated scheme (from
                    ! faces to faces)
                    !   o       o       o                   o       o       o
                    !   |_______|_______|____  .....  ______|_______|_______|
                    !           |       |                   |       |
                    !           \\      \\                  \\      \\

                    detBband = [min(detBband(1), -2), max(detBband(2), +2)]

                CASE ('D1_f2c_3E_L')
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !       |
                    !       \

                    detBband = [min(detBband(1), 0), max(detBband(2), +3)]

                CASE ('D1_c2f_4C_L')
                    !       o       o
                    !   ____|_______|____
                    !   |       |       |
                    !   \       \       \

                    detBband = [min(detBband(1), -1), max(detBband(2), 0)]

                CASE ('D1_c2f_4C_LL')
                    !   o   o       o       o
                    !   |___|_______|_______|
                    !   |       |
                    !   \       \

                    detBband = [min(detBband(1), -1), max(detBband(2), 2)]

                CASE ('D2_c2c_5C_L','D2_f2f_5C_L')
                    !   o       o       o       o       o       o
                    !   |_______|_______|_______|_______|_______|
                    !           |       |
                    !           \\      \\

                    detBband = [min(detBband(1), -1), max(detBband(2), +4)]


                CASE ('D1_f2c_3E_R')
                    !   o       o       o       o
                    !   |_______|_______|_______|
                    !                       |
                    !                       \

                    detBband = [min(detBband(1), -2), max(detBband(2), +1)]


                CASE ('D1_c2f_4C_R')
                    !       o       o
                    !   ____|_______|____
                    !   |       |       |
                    !   \       \       \

                    detBband = [min(detBband(1), -1), max(detBband(2), 0)]


                CASE ('D1_c2f_4C_RR')
                    !   o       o       o   o
                    !   |_______|_______|___|
                    !               |       |
                    !               \       \

                    detBband = [min(detBband(1), -3), max(detBband(2), 0)]


                CASE ('D2_c2c_5C_R','D2_f2f_5C_R')
                    !   o       o       o       o       o       o
                    !   |_______|_______|_______|_______|_______|
                    !                           |       |
                    !                           \\      \\

                    detBband = [min(detBband(1), -4), max(detBband(2), +1)]

                CASE ('zero')
                    ! No boundary scheme

                CASE DEFAULT

                    PRINT *, 'Error in compact_mod.f90: no such scheme defined.'
                    PRINT *, schemes(i)
                    CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                    STOP

            END SELECT scheme

        END DO for_each_scheme

    END FUNCTION detBband

END MODULE library
