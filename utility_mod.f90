MODULE utility

IMPLICIT NONE

INTERFACE row_reorder
  MODULE PROCEDURE row_reorder_REAL
  MODULE PROCEDURE row_reorder_INT
END INTERFACE row_reorder


PRIVATE
PUBLIC :: sort, count_zeros, row_reorder

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION count_zeros(A) RESULT(nz)

    IMPLICIT NONE

    REAL, DIMENSION(:, :), INTENT(IN)    :: A
    INTEGER, DIMENSION(:), ALLOCATABLE   :: nz
    INTEGER                              :: i, j
    INTEGER, DIMENSION(2)                :: dmn

    dmn = SHAPE(A)

    ALLOCATE(nz(dmn(1)))

    nz = 0
    rows : DO i = 1, dmn(1)
      columns : DO j = 1, dmn(2)
        IF (A(i, j)/=0) THEN
          nz(i) = nz(i) + 1
        END IF
      END DO columns
    END DO rows

  END FUNCTION count_zeros
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE sort(a, ord)

    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN OUT)    :: a
    INTEGER, DIMENSION(SIZE(a)), INTENT(OUT) :: ord
    INTEGER                                  :: i, j
    INTEGER                                  :: dum

    ord = [(i, i=1, SIZE(ord))]

    DO i = 1, SIZE(a)
      DO j = i+1, SIZE(a)
        IF (a(j)>a(i)) THEN
          dum = a(i)
          a(i) = a(j)
          a(j) = dum

          dum = ord(i)
          ord(i) = ord(j)
          ord(j) = dum
        END IF
      END DO
    END DO

  END SUBROUTINE sort
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION row_reorder_REAL(A, ord) RESULT(B)

    IMPLICIT NONE

    REAL, DIMENSION(:, :), INTENT(IN)          :: A
    INTEGER, DIMENSION(:), INTENT(IN)          :: ord
    REAL, DIMENSION(SIZE(A, 1), SIZE(A, 2))    :: B
    INTEGER                                    :: i

    IF (SIZE(ord)/=SIZE(A, 1)) THEN
      PRINT *, 'dimensions not coherent'
    ELSE

      DO i = 1, SIZE(A, 1)
        B(i, :) = A(ord(i), :)
      END DO

    END IF

  END FUNCTION row_reorder_REAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION row_reorder_INT(A, ord) RESULT(B)

    IMPLICIT NONE

    INTEGER, DIMENSION(:, :), INTENT(IN)        :: A
    INTEGER, DIMENSION(:), INTENT(IN)           :: ord
    INTEGER, DIMENSION(SIZE(A, 1), SIZE(A, 2))  :: B
    INTEGER                                     :: i

    IF (SIZE(ord)/=SIZE(A, 1)) THEN
      PRINT *, 'dimensions not coherent'
    ELSE

      DO i = 1, SIZE(A, 1)
        B(i, :) = A(ord(i), :)
      END DO

    END IF

  END FUNCTION row_reorder_INT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



END MODULE utility
