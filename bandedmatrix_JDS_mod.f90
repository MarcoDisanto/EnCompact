MODULE bandedmatrix_JDS

!Modified version of bandedmatrix_mod.f90 module designed to work with matrices stored by means of the JDS format
!More informations can be found at http://netlib.org/linalg/html_templates/node95.html

USE utility
USE bandedmatrix

IMPLICIT NONE

TYPE JDS
  INTEGER                            :: n_row, n_col   !total number of rows and columns
  INTEGER, DIMENSION(:), ALLOCATABLE :: perm           !original row ordering
  REAL,    DIMENSION(:), ALLOCATABLE :: jdiag          !vector of jagged diagonals packed column-wise
  INTEGER, DIMENSION(:), ALLOCATABLE :: col_ind        !same size of jdiag; memorizes the column index for each non-null element
  INTEGER, DIMENSION(:), ALLOCATABLE :: jptr           !vector of length equal to the max number of non-null elements on a row;
                                                       !contains the indices of the first elements of each column of the jagged matrix,
                                                       !thus letting distinguish different columns even after the "packaging" operation.
                                                       !ATTENTION : originally it was supposed to be a pointer
END TYPE JDS

PRIVATE
PUBLIC :: JDS, mat2JDS, CDS2JDS, JDS_mat_x_vec

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION mat2JDS(A)
    !This function coverts a regularly stored matrix into JDS format

    IMPLICIT NONE

    REAL, DIMENSION(:, :), INTENT(IN)          :: A
    REAL, DIMENSION(SIZE(A, 1), SIZE(A, 2))    :: Adum
    INTEGER, DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: C, Cdum
    TYPE(JDS)                                  :: mat2JDS
    INTEGER                                    :: i
    INTEGER, DIMENSION(SIZE(A, 1))             :: nz

    !dimensioni matrice
    mat2JDS%n_row = SIZE(A, 1)
    mat2JDS%n_col = SIZE(A, 2)

    !numero di elementi non nulli in ogni riga
    nz = count_zeros(A)
    !riordinamento di nz
    ALLOCATE(mat2JDS%perm(mat2JDS%n_row))
    CALL sort(nz, mat2JDS%perm)
    !creazione matrice con indici di colonna
    C = 0
    WHERE (A/=0) C = 1
    DO i = 1, SIZE(C, 2)
      C(:, i) = i*C(:, i)
    END DO

    !!!!!creazione vettore elementi e indici di colonna!!!!!!!!!!!!!!!!!!!!!!!!!
    Adum = 0
    Cdum = 0
    DO i = 1, mat2JDS%n_row
      Adum(i, 1:nz(i)) = PACK(A(i, :), A(i, :)/=0)
      Cdum(i, 1:nz(i)) = PACK(C(i, :), C(i, :)/=0)
    END DO
    !riordinamento righe
    Adum = row_reorder(Adum, mat2JDS%perm)
    Cdum = row_reorder(Cdum, mat2JDS%perm)
    !vettore elementi
    ALLOCATE(mat2JDS%jdiag(SUM(nz)))
    mat2JDS%jdiag = PACK(Adum, Adum/=0)
    !vettore indici di colonna
    ALLOCATE(mat2JDS%col_ind(SUM(nz)))
    mat2JDS%col_ind = PACK(Cdum, Cdum/=0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!creazione vettore inizio colonna!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(mat2JDS%jptr(MAXVAL(nz)))

    i = 1
    DO WHILE (Adum(1, i)/=0)
      mat2JDS%jptr(i) = COUNT(Adum(:, 1:i-1)/=0) + 1
      i = i+1
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END FUNCTION mat2JDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CDS2JDS(A)
    !This function coverts a CDS format into a JDS one
    !ATTENTION : CDS supposed to store diagonals as columns!!!
    !CDS type needs an additional parameter to distinguish between column-wise
    !and row-wise storage (tip)

    IMPLICIT NONE

    TYPE(CDS), INTENT(IN)                                    :: A
    TYPE(JDS)                                                :: CDS2JDS
    INTEGER                                                  :: i
    INTEGER, DIMENSION(SIZE(A%matrix, 1))                    :: nz
    INTEGER, DIMENSION(SIZE(A%matrix, 1), -A%ld:A%ud)        :: C
    REAL, DIMENSION(SIZE(A%matrix, 1), SIZE(A%matrix, 2))    :: Adum
    INTEGER, DIMENSION(SIZE(A%matrix, 1), SIZE(A%matrix, 2)) :: Cdum


    CDS2JDS%n_row = (A%ub(1)-A%lb(1))+1
    CDS2JDS%n_col = (A%ub(2)-A%lb(2))+1

    !numero di elementi non nulli in ogni riga
    nz = count_zeros(A%matrix)
    !riordinamento di nz
    ALLOCATE(CDS2JDS%perm(CDS2JDS%n_row))
    CALL sort(nz, CDS2JDS%perm)

    !creazione matrice con indici di colonna
    C = 0
    !WHERE (A%matrix/=0) C = 1
    C(:, 0) = [(i, i=1,SIZE(C, 1))]
    DO i = -A%ld, A%ud
      C(:, i) = C(:, 0)+i
    END DO
    WHERE (A%matrix==0) C = 0                           !prima fase di eliminazione elementi superflui
    WHERE ((C<=0) .OR. (C>(A%ub(2)-A%lb(2)+1))) C = 0   !seconda fase di eliminazione elementi superflui

    !!!!!creazione vettore elementi e indici di colonna!!!!!!!!!!!!!!!!!!!!!!!!!
    Adum = 0
    Cdum = 0
    DO i = 1, CDS2JDS%n_row
      Adum(i, 1:nz(i)) = PACK(A%matrix(i, :), A%matrix(i, :)/=0)
      Cdum(i, 1:nz(i)) = PACK(C(i, :), C(i, :)/=0)
    END DO
    !riordinamento righe
    Adum = row_reorder(Adum, CDS2JDS%perm)
    Cdum = row_reorder(Cdum, CDS2JDS%perm)
    !vettore elementi
    ALLOCATE(CDS2JDS%jdiag(SUM(nz)))
    CDS2JDS%jdiag = PACK(Adum, Adum/=0)
    !vettore indici di colonna
    ALLOCATE(CDS2JDS%col_ind(SUM(nz)))
    CDS2JDS%col_ind = PACK(Cdum, Cdum/=0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!creazione vettore inizio colonna!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(CDS2JDS%jptr(MAXVAL(nz)))

    i = 1
    DO WHILE (Adum(1, i)/=0)
      CDS2JDS%jptr(i) = COUNT(Adum(:, 1:i-1)/=0) + 1
      i = i+1
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END FUNCTION CDS2JDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION JDS_mat_x_vec(A, v)
  !function performing the matrix product between a JDS stored matrix and a vector

  IMPLICIT NONE

  TYPE(JDS), INTENT(IN)              :: A
  REAL, DIMENSION(:), INTENT(IN)     :: v
  REAL, DIMENSION(A%n_row)           :: JDS_mat_x_vec
  INTEGER                            :: i, j, row=0
  !INTEGER, DIMENSION(:), ALLOCATABLE :: row  !vettore con gli indici di riga

  ! check dimensions
  IF (A%n_col/=SIZE(v)) THEN
      PRINT *, 'ERROR in bandedmatrix_mod_JDS.f90: matrix dimensions must agree in JDS_mat_x_vec!'
      STOP
  END IF

  !vettore indici di riga
  !ALLOCATE(row(SIZE(A%jdiag)))
  !DO i = 1, SIZE(A%jptr)
  !  row(A%jptr(i)) = 1
  !  DO j = A%jptr(i)+1, A%jptr(min(i+1, size(A%jptr)))-1
  !    row(j) = row(j-1) + 1
  !  END DO
  !END DO

  !prodotto
  JDS_mat_x_vec = 0
  DO i = 1, SIZE(A%jdiag)

    DO j = 1, SIZE(A%jptr)
      if (A%jptr(j)==i) row = 1
    END DO

    JDS_mat_x_vec(A%perm(row)) = JDS_mat_x_vec(A%perm(row)) + A%jdiag(i)*v(A%col_ind(i))

    row = row + 1
  END DO

  END FUNCTION JDS_mat_x_vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




END MODULE bandedmatrix_JDS
