MODULE bandedmatrix

    IMPLICIT NONE

    TYPE CDS ! CDS (Compressed Diagonal Storage) type
        ! TODO: improve the following comment!!!
        ! Lower bounds of both directions are allowed to be not 1. For what concerns the 2nd direction
        ! the reason is that rectangular matrices
        ! can come from the partitioning of (rectangular as well as square) matrices; a consequence of the
        ! partitioning is that the vector which is to be multiplied by the partition of the matrix has to
        ! be provided with overlap values; in the same way, the matrix has to contain the coefficients that
        ! are to be multiplied by these values.
        ! Cyclic matrices can be thought of as rectangular matrices, as long as the vector which is going
        ! to be multiplied by them is long enough and filled with overlap values (received by neighboring
        ! processes, in a cyclic sense).
        REAL,    DIMENSION(:,:), ALLOCATABLE :: matrix ! sparse banded matrix (diagonals written as columns)
        INTEGER, DIMENSION(2)                :: lb, ub ! matrix's lower and upper bounds along the two dimensions
        INTEGER                              :: ld, ud ! number of lower diagonals, number of upper diagonals
        INTEGER                              :: lid, uid ! as ld and ud, but considering only the inner rows
        ! Note that  lbound(matrix,1) and ubound(matrix,1) should be equal to lb(1) and ub(1) respectively.
        ! Similarly, lbound(matrix,2) and ubound(matrix,2) should be equal to -ld and +ud.
    END TYPE CDS

    INTERFACE OPERATOR(*)
        MODULE PROCEDURE CDS_mat_x_full_vec
        MODULE PROCEDURE CDS_mat_x_CDS_mat
    END INTERFACE

    PRIVATE ! hides items not listed on the following public statement
    PUBLIC :: CDS, &
              OPERATOR (*)

CONTAINS

    FUNCTION CDS_mat_x_CDS_mat(A, B)
    ! This function performs the matrix priduct between the full counterparts
    ! of A and B.
    ! IMPORTANT: the upper elements of lower diagonals and, if A represent a square
    ! matrix, the lower elements of upper diagonals are ignored. Such elements could be
    ! used to perform the product in the case of A representing a cyclic banded matrix.
    ! In this sense either a new function associated with a new derived data type could be
    ! written, or this function could be modified and the CDS type provided with a
    ! logical field identifying the represented matrix to be cyclic or not.

        USE MPI

        IMPLICIT NONE

        ! IN/OUT variables
        TYPE(CDS), INTENT(IN) :: A, B
        TYPE(CDS)             :: CDS_mat_x_CDS_mat

        ! internal variables
        INTEGER :: errorcode, ierr, i, j, k

        ! check dimensions
        IF (A%lb(2) /= B%lb(1) .OR. A%ub(2) /= B%ub(1)) THEN
            PRINT *, 'ERROR in bandedmatrix_mod.f90: matrix dimensions must agree in CDS_mat_x_CDS_mat!'
            PRINT *, 'bounds of A are [', A%lb(1), A%ub(1), ']x[', A%lb(2), A%ub(2), ']'
            PRINT *, 'bounds of B are [', B%lb(1), B%ub(1), ']x[', B%lb(2), B%ub(2), ']'
            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            STOP
        END IF

        CDS_mat_x_CDS_mat%lb(1) = A%lb(1)
        CDS_mat_x_CDS_mat%ub(1) = A%ub(1)
        CDS_mat_x_CDS_mat%lb(2) = B%lb(2)
        CDS_mat_x_CDS_mat%ub(2) = B%ub(2)
        CDS_mat_x_CDS_mat%ld = A%ld + B%ld
        CDS_mat_x_CDS_mat%ud = A%ud + B%ud
        ALLOCATE(CDS_mat_x_CDS_mat%matrix(CDS_mat_x_CDS_mat%lb(1):CDS_mat_x_CDS_mat%ub(1), &
                                          & -CDS_mat_x_CDS_mat%ld:+CDS_mat_x_CDS_mat%ud))

        CDS_mat_x_CDS_mat%matrix = 0
        ! TODO: il primo ciclo può essere fatto su un insieme di indici più
        ! stretto, se per caso la prima matrice ha le prime o ultime colonne
        ! nulle
        for_each_row: DO i = CDS_mat_x_CDS_mat%lb(1), CDS_mat_x_CDS_mat%ub(1)

            for_each_column: DO j = MAX(-CDS_mat_x_CDS_mat%ld, CDS_mat_x_CDS_mat%lb(2)-i), &
                                  & MIN(+CDS_mat_x_CDS_mat%ud, CDS_mat_x_CDS_mat%ub(2)-i)

                for_each_non_null_common_term: DO k = max(i-a%ld, j-b%ud+i, a%lb(2)), min(i+a%ud, j+b%ld+i, a%ub(2))

                    CDS_mat_x_CDS_mat%matrix(i, j) = CDS_mat_x_CDS_mat%matrix(i, j) + A%matrix(i, k-i)*B%matrix(k, j-k+i)

                END DO for_each_non_null_common_term

            END DO for_each_column

        END DO for_each_row

    END FUNCTION CDS_mat_x_CDS_mat


    FUNCTION CDS_mat_x_full_vec(A, v)
    ! This function performs the matrix priduct between the full counterpart of A
    ! and the full vector v.
    ! IMPORTANT: the upper elements of lower diagonals and, if A represent a square
    ! matrix, the lower elements of upper diagonals are ignored. Such elements could be
    ! used to perform the product in the case of A representing a cyclic banded matrix.
    ! In this sense either a new function associated with a new derived data type could be
    ! written, or this function could be modified and the CDS type provided with a
    ! logical field identifying the represented matrix to be cyclic or not.

        USE MPI, ONLY : MPI_COMM_WORLD

        ! input arguments
        TYPE(CDS),                       INTENT(IN) :: A ! matrix in its sparse-banded form
        REAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: v ! full array multiplying A

        ! output argument
        REAL, DIMENSION(A%lb(1):A%ub(1)) :: CDS_mat_x_full_vec ! result

        ! internal arguments
        INTEGER :: lband, uband ! number of lower diagonals, number of upper diagonals (both supposed >0)
        INTEGER :: lv, uv       ! lower and upper indices of the array v (supposed to be lv<=0, uv>0)
        INTEGER :: i            ! generic index

        INTEGER :: errorcode, ierr

        uband  = A%ud ! it's positive when super-diagonals are present
        lband  = A%ld ! it's positive when   sub-diagonals are present
        uv  = ubound(v,1)
        lv  = lbound(v,1)

        ! check dimensions
        IF (ANY([lv, uv] /= [A%lb(2), A%ub(2)])) THEN
            PRINT *, 'ERROR in bandedmatrix_mod.f90: matrix dimensions must agree in CDS_mat_x_full_vec!'
            CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            STOP
        END IF

        ! product
        DO i = A%lb(1), A%ub(1)
            CDS_mat_x_full_vec(i) = DOT_PRODUCT(A%matrix(i,max(-lband,lv-i):min(uband,uv-i)), &
                                                       & v(max(i-lband,lv):min(i+uband,uv)))
        END DO

    END FUNCTION CDS_mat_x_full_vec

END MODULE bandedmatrix
