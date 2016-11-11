MODULE essentials

    IMPLICIT NONE

    INTERFACE lininterpmid1D
        MODULE PROCEDURE midpoint_1D
        MODULE PROCEDURE midpoint_2D
        MODULE PROCEDURE midpoint_3D
        MODULE PROCEDURE midpoint_4D
        !MODULE PROCEDURE midpoint_5D
        !MODULE PROCEDURE midpoint_6D
        !MODULE PROCEDURE midpoint_7D
    END INTERFACE lininterpmid1D

    INTERFACE printmatrix
        MODULE PROCEDURE print2DmatrixINT
        MODULE PROCEDURE print3DmatrixINT
        MODULE PROCEDURE print1DmatrixREAL
        MODULE PROCEDURE print2DmatrixREAL
        MODULE PROCEDURE print3DmatrixREAL
    END INTERFACE printmatrix

    PRIVATE ! hides items not listed in the public statement 
    PUBLIC :: lininterpmid1D, &
              sub2ind, ind2sub, &
              linspace, logspace, &
              logical2integer, logical2real, &
              printmatrix, fliplr, extract_elements, KronDelta
              
CONTAINS

    INTEGER FUNCTION KronDelta(i,j) 
    ! Kronecker delta function between integers

        IMPLICIT NONE

        INTEGER :: i, j 
        KronDelta = 0
        IF (i == j) KronDelta = 1 

    END FUNCTION KronDelta


    ! FIXME: this function is maybe useless, since Fortran should already support vector-valued indices,
    ! as it is said here -> http://www.moreisdifferent.com/2015/07/16/why-physicsts-still-use-fortran/
    FUNCTION extract_elements(array, indices) RESULT(extracted)

        IMPLICIT NONE

        ! IN/OUT variables
        REAL,    DIMENSION(:), INTENT(IN),  ALLOCATABLE :: array ! deferred shape array
        INTEGER, DIMENSION(:), INTENT(IN)               :: indices
        REAL,    DIMENSION(SIZE(indices))               :: extracted

        ! internal variables
        INTEGER :: i

        DO i = 1, SIZE(indices)
            extracted(i) = array(indices(i))
        END DO

    END FUNCTION extract_elements


    ELEMENTAL INTEGER FUNCTION logical2integer(a)
        LOGICAL, INTENT(IN) :: a
        IF (a) THEN
          logical2integer = 1
        ELSE
          logical2integer = 0
        END IF
    END FUNCTION logical2integer
    

    ELEMENTAL REAL FUNCTION logical2real(a)
        LOGICAL, INTENT(IN) :: a
        IF (a) THEN
          logical2real = 1.0
        ELSE
          logical2real = 0.0
        END IF
    END FUNCTION logical2real
    

    FUNCTION linspace(L1,L2,N) RESULT(V)

        ! variables in argument list
        REAL,    INTENT(IN)              :: L1,L2
        INTEGER, INTENT(IN)              :: N
        REAL,               DIMENSION(N) :: V
        
        ! internal variables
        INTEGER :: i
        
        V = L1 + [(i, i = 0,N-1)]*(L2-L1)/(N-1)

        RETURN

    END FUNCTION linspace


    FUNCTION logspace(logL1,logL2,N)

        ! variables in argument list
        
        REAL,    INTENT(IN) :: logL1,logL2
        INTEGER, INTENT(IN) :: N
        REAL,  DIMENSION(N) :: logspace
        
        logspace = 10**linspace(logL1,logL2,N)

        RETURN

    END FUNCTION logspace


    FUNCTION sub2ind(shape,subscripts)
        ! Linear index from multiple subscripts.
        ! sub2ind is used to determine the equivalent single index
        ! corresponding to a given set of subscript values.

        IMPLICIT NONE

        ! output variables (the function itself)
        INTEGER :: sub2ind
        
        ! input variables
        INTEGER, DIMENSION(:), INTENT(IN) :: subscripts, shape
        
        ! internal variables
        INTEGER :: r, i, index
        
        ! rank of the array whose shape is shape
        r = size(shape)

        ! initialize 
        index = 1

        ! cycle
        DO i = 1,r
            index = index + (subscripts(i) - 1)*product(shape(1:i-1))
        END DO

        ! assign the output to the function itself
        sub2ind = index

    END FUNCTION sub2ind


    FUNCTION ind2sub(shape,index)
        ! Multiple subscripts from linear index.
        ! ind2sub is used to determine the equivalent subscript values
        ! corresponding to a given single index into an array.

        IMPLICIT NONE

        ! input variables
        INTEGER, DIMENSION(:), INTENT(IN) :: shape
        INTEGER, INTENT(IN) :: index

        ! output variables (the function itself)
        INTEGER, DIMENSION(size(shape)) :: ind2sub

        ! internal variables
        INTEGER :: r, i, dividend, divisor, quotient, remainder

        ! rank of the array whose shape is shape
        r = size(shape)

        ! initialize
        ind2sub = shape
        dividend = index

        ! cycle
        DO i = r,1,-1
            divisor = product(shape(1:i-1))
            quotient = dividend/divisor
            remainder = dividend - quotient*divisor
            IF (remainder == 0) THEN
                ind2sub(i) = quotient
                RETURN
            ELSE
                ind2sub(i) = quotient + 1
                dividend = remainder
            END IF
        END DO

    END FUNCTION ind2sub
    

    SUBROUTINE print2DmatrixINT(A)

        IMPLICIT NONE

        INTEGER, DIMENSION(:,:) :: A
        INTEGER :: ix

        DO ix = 1,size(A,1)
            WRITE(*,*) A(ix,:)
        END DO

    END SUBROUTINE print2DmatrixINT


    SUBROUTINE print3DmatrixINT(A)

        IMPLICIT NONE

        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: A
        INTEGER, DIMENSION(3) :: lb, ub
        INTEGER :: ix, iz
        CHARACTER(len=3), DIMENSION(0:4) :: ordinal

        ordinal(0) = '-th'
        ordinal(1) = '-st'
        ordinal(2) = '-nd'
        ordinal(3) = '-rd'
        ordinal(4) = '-th'

        lb = LBOUND(A)
        ub = UBOUND(A)
        DO iz = lb(3), ub(3)
            WRITE(*,*) iz, ordinal(MIN(ABS(iz),4)), " Z slice (colums indexed from",lb(2),"to",ub(2),")"
            DO ix = lb(1), ub(1)
                IF (ix == 1) THEN
                    WRITE(*,*) '#1 ->', A(ix,:,iz)
                ELSE
                    WRITE(*,*) '     ', A(ix,:,iz)
                END IF
            END DO
        END DO

    END SUBROUTINE print3DmatrixINT


    SUBROUTINE print1DmatrixREAL(v)

        IMPLICIT NONE

        REAL, DIMENSION(:) :: v
        INTEGER :: ix

        DO ix = 1,size(v)
            WRITE(*,*) v(ix)
        END DO

    END SUBROUTINE print1DmatrixREAL


    SUBROUTINE print2DmatrixREAL(A)

        IMPLICIT NONE

        REAL, DIMENSION(:,:) :: A
        INTEGER :: ix

        DO ix = 1,size(A,1)
            WRITE(*,*) A(ix,:)
        END DO

    END SUBROUTINE print2DmatrixREAL


    SUBROUTINE print3DmatrixREAL(A)

        IMPLICIT NONE

        REAL, DIMENSION(:,:,:), ALLOCATABLE :: A
        INTEGER, DIMENSION(3) :: lb, ub
        INTEGER :: ix, iz
        CHARACTER(len=3), DIMENSION(0:4) :: ordinal

        ordinal(0) = '-th'
        ordinal(1) = '-st'
        ordinal(2) = '-nd'
        ordinal(3) = '-rd'
        ordinal(4) = '-th'

        lb = LBOUND(A)
        ub = UBOUND(A)
        DO iz = lb(3), ub(3)
            WRITE(*,*) iz, ordinal(MIN(ABS(iz),4)), " Z slice (colums indexed from",lb(2),"to",ub(2),")"
            DO ix = lb(1), ub(1)
                IF (ix == 1) THEN
                    WRITE(*,*) '#1 ->', A(ix,:,iz)
                ELSE
                    WRITE(*,*) '     ', A(ix,:,iz)
                END IF
            END DO
        END DO

    END SUBROUTINE print3DmatrixREAL


    FUNCTION midpoint_1D(A,dir)
        IMPLICIT NONE
        ! input variables
        REAL, DIMENSION(:), INTENT(IN) :: A
        INTEGER :: dir

        ! output variable (the function itself)
        REAL, DIMENSION(size(A)-1) :: midpoint_1D

        ! internal variables
        INTEGER :: siz, i1


        IF (dir/=1) THEN
            WRITE(*,*) "codice incompleto"
            ! E' evidente che per array 1D la funzione può agire in una sola direzione.
            ! Tuttavia si potrebbe implementare la stessa logica di MATLAB, nel quale,
            ! se dir è una direzione lungo la quale l'estensione è 1, l'output è un
            ! vettore vuoto.
            STOP
        END IF

        siz = size(A)
        DO i1 = 1, siz-1
            midpoint_1D(i1) = .5*(A(i1) + A(i1+1))
        END DO
    END FUNCTION midpoint_1D


    FUNCTION midpoint_2D(A,dim)

        IMPLICIT NONE
        ! input variables
        REAL, DIMENSION(:,:), INTENT(IN) :: A
        INTEGER :: dim

        ! output variable (the function itself)
        REAL, DIMENSION(:,:), ALLOCATABLE :: midpoint_2D

        ! internal variables
        INTEGER, DIMENSION(2) :: shap
        INTEGER :: i1, i2


        shap = shape(A)
        
        SELECT CASE(dim)

            CASE(1)
                ALLOCATE(midpoint_2D(shap(1)-1,shap(2)))
                DO i2 = 1,shap(2)
                    DO i1 = 1,shap(1)-1
                        midpoint_2D(i1,i2) = .5*(A(i1+1,i2) + A(i1,i2))
                    END DO
                END DO

            CASE(2)
                ALLOCATE(midpoint_2D(shap(1),shap(2)-1))
                DO i1 = 1,shap(1)
                    DO i2 = 1,shap(2)-1
                        midpoint_2D(i1,i2) = .5*(A(i1,i2) + A(i1,i2+1))
                    END DO
                END DO

            CASE DEFAULT
                WRITE(*,*) "codice incompleto"
            ! E' evidente che per array 2D la funzione può agire in due direzioni.
            ! Tuttavia si potrebbe implementare la stessa logica di MATLAB, nel quale,
            ! se dim è una direzione lungo la quale l'estensione è 1, l'output è un
            ! vettore vuoto.
                STOP

        END SELECT
    END FUNCTION midpoint_2D


    FUNCTION midpoint_3D(A,dim)

        IMPLICIT NONE
        ! input variables
        REAL, DIMENSION(:,:,:), INTENT(IN) :: A
        INTEGER :: dim

        ! output variable (the function itself)
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: midpoint_3D

        ! internal variables
        INTEGER, DIMENSION(3) :: shap
        INTEGER :: i1, i2, i3


        shap = shape(A)
        
        SELECT CASE(dim)

            CASE(1)
                ALLOCATE(midpoint_3D(shap(1)-1,shap(2),shap(3)))
                DO i2 = 1,shap(2)
                    DO i3 = 1,shap(3)
                        DO i1 = 1,shap(1)-1
                            midpoint_3D(i1,i2,i3) = .5*(A(i1+1,i2,i3) + A(i1,i2,i3))
                        END DO
                    END DO
                END DO

            CASE(2)
                ALLOCATE(midpoint_3D(shap(1),shap(2)-1,shap(3)))
                DO i1 = 1,shap(1)
                    DO i3 = 1,shap(3)
                        DO i2 = 1,shap(2)-1
                            midpoint_3D(i1,i2,i3) = .5*(A(i1,i2+1,i3) + A(i1,i2,i3))
                        END DO
                    END DO
                END DO

            CASE(3)
                ALLOCATE(midpoint_3D(shap(1),shap(2),shap(3)-1))
                DO i1 = 1,shap(1)
                    DO i2 = 1,shap(2)
                        DO i3 = 1,shap(3)-1
                            midpoint_3D(i1,i2,i3) = .5*(A(i1,i2,i3+1) + A(i1,i2,i3))
                        END DO
                    END DO
                END DO

            CASE DEFAULT
                WRITE(*,*) "codice incompleto"
            ! E' evidente che per array 2D la funzione può agire in due direzioni.
            ! Tuttavia si potrebbe implementare la stessa logica di MATLAB, nel quale,
            ! se dim è una direzione lungo la quale l'estensione è 1, l'output è un
            ! vettore vuoto.
                STOP

        END SELECT
    END FUNCTION midpoint_3D


    FUNCTION midpoint_4D(A,dim)

        IMPLICIT NONE
        ! input variables
        REAL, DIMENSION(:,:,:,:), INTENT(IN) :: A
        INTEGER :: dim

        ! output variable (the function itself)
        REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: midpoint_4D

        ! internal variables
        INTEGER, DIMENSION(4) :: shap
        INTEGER :: i1, i2, i3, i4


        shap = shape(A)
        
        SELECT CASE(dim)

            CASE(1)
                ALLOCATE(midpoint_4D(shap(1)-1,shap(2),shap(3),shap(4)))
                DO i2 = 1,shap(2)
                    DO i3 = 1,shap(3)
                        DO i4 = 1,shap(4)
                            DO i1 = 1,shap(1)-1
                                midpoint_4D(i1,i2,i3,i4) = .5*(A(i1+1,i2,i3,i4) + A(i1,i2,i3,i4))
                            END DO
                        END DO
                    END DO
                END DO

            CASE(2)
                ALLOCATE(midpoint_4D(shap(1),shap(2)-1,shap(3),shap(4)))
                DO i1 = 1,shap(1)
                    DO i3 = 1,shap(3)
                        DO i4 = 1,shap(4)
                            DO i2 = 1,shap(2)-1
                                midpoint_4D(i1,i2,i3,i4) = .5*(A(i1,i2+1,i3,i4) + A(i1,i2,i3,i4))
                            END DO
                        END DO
                    END DO
                END DO

            CASE(3)
                ALLOCATE(midpoint_4D(shap(1),shap(2),shap(3)-1,shap(4)))
                DO i1 = 1,shap(1)
                    DO i2 = 1,shap(2)
                        DO i4 = 1,shap(4)
                            DO i3 = 1,shap(3)-1
                                midpoint_4D(i1,i2,i3,i4) = .5*(A(i1,i2,i3+1,i4) + A(i1,i2,i3,i4))
                            END DO
                        END DO
                    END DO
                END DO

            CASE(4)
                ALLOCATE(midpoint_4D(shap(1),shap(2),shap(3),shap(4)-1))
                DO i1 = 1,shap(1)
                    DO i2 = 1,shap(2)
                        DO i3 = 1,shap(3)
                            DO i4 = 1,shap(4)-1
                                midpoint_4D(i1,i2,i3,i4) = .5*(A(i1,i2,i3,i4+1) + A(i1,i2,i3,i4))
                            END DO
                        END DO
                    END DO
                END DO

            CASE DEFAULT
                WRITE(*,*) "codice incompleto"
            ! E' evidente che per array 2D la funzione può agire in due direzioni.
            ! Tuttavia si potrebbe implementare la stessa logica di MATLAB, nel quale,
            ! se dim è una direzione lungo la quale l'estensione è 1, l'output è un
            ! vettore vuoto.
                STOP

        END SELECT
    END FUNCTION midpoint_4D


    FUNCTION fliplr(v_dummy)

        IMPLICIT NONE
        REAL,ALLOCATABLE, DIMENSION(:) :: fliplr
        REAL,  DIMENSION(:):: v_dummy
        INTEGER:: siz

        ! Only for one-dimensional vectors!!!!!!!!
        siz=SIZE(v_dummy)

        ALLOCATE(fliplr(siz))

        fliplr(siz:1:-1)=v_dummy(1:siz)

    END FUNCTION fliplr

END MODULE essentials
