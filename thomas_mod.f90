MODULE Thomas_suite

    PRIVATE
    PUBLIC :: Thomas, cypent, pentdag

CONTAINS

SUBROUTINE Thomas(As,b,x)
! Algoritmo di Thomas (propriamente detto).
!
! Il sistema tridiagonale
!                           A x = b
! viene passato alla subroutine mediante
!   -   la matrice As ('s' sta per sparsa), che ha dimensioni (N,3) se
!       A ha dimensioni (N,N) ed costituita dalle sole diagonali A
!   -   il termine noto b
!
! Gli elementi di As corrispondono a quelli non nulli di A, nel modo seguente
!
!    As(1,0)  As(1,+1)     0        0        0      0    ...   0           0          0
!    As(2,-1) As(2,0)  As(2,+1)     0        0      0    ...   0           0          0
!       0     As(3,-1) As(3,0)  As(3,+1)     0      0    ...   0           0          0
!       0       0      As(3,-1) As(3,0)  As(3,+1)   0    ...   0           0          0
!       :       :
!       :       :
!       0       :                                         As(N-1,-1)   As(N-1,0)   As(N-1,+1)
!       0       0                                                      As(N,-1)    As(N,0)
!
! Si osserva esplicitamente che, avendo le subdiagonali
! un numero di N-1 elementi (N è il numero di elementi della
! diagonale principale), i valori degli elementi di tali subdiagonali
! sono memorizzati in modo che As(1,-1) e As(N,+1) non trovino posto nella
! matrice e quindi non devono essere utilizzati nelle elaborazioni.
! Tale scelta ha lo scopo di semplificare le notazioni relative
! appunto agli elementi appartenenti ad una stessa riga.


REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)  :: As ! deferred-shape array (so that the bounds are automatically passed)
REAL, DIMENSION(:),                INTENT(IN)  :: b
REAL, DIMENSION(:),                INTENT(OUT) :: x
INTEGER                                        :: lb1, ub1  ! row bounds are needed to generalise the algorithm

INTEGER                           :: N, i
REAL, DIMENSION(:,:), ALLOCATABLE :: LU
REAL, DIMENSION(:),   ALLOCATABLE :: y

! -----------------Allocazioni---------------------------------------------
N = SIZE(As,1)
ALLOCATE(LU(N,-1:1))
ALLOCATE(y(N))

! -----------------Bounds--------------------------------------------------
lb1 = LBOUND(As, 1)
ub1 = UBOUND(As, 1)

! -----------------Algoritmo di Thomas (decomposizione LU)-----------------
! La matrice LU potrebbe anche essere proprio la As, dichiarata con
! INTENT(INOUT) e modificata di volta in volta (il che rappresenterebbe
! anche un bel risparmio di memoria!).
LU(1,0)  = As(lb1,0)
LU(1,+1) = As(lb1,+1)
y(1) = b(1)
DO i = 2,N
    LU(i,-1) = As(lb1 + i-1,-1)/LU(i-1,0)
    LU(i,0)  = As(lb1 + i-1,0) - LU(i,-1)*As(lb1 + i-2,+1)
    LU(i,+1) = As(lb1 + i-1,+1)
    y(i) = b(i)-LU(i,-1)*y(i-1)
END DO
! -----------------Algoritmo di Thomas (back substitution)-----------------
x(N) = y(N)/LU(N,0)
DO i = N-1,1,-1
    x(i) = (y(i)-LU(i,1)*x(i+1))/LU(i,0)
END DO

END SUBROUTINE Thomas




SUBROUTINE ThomasC(As,b,x)
! Algoritmo di Thomas per il caso circolante.
!
! Il sistema tridiagonale circolante
!                           A x = b
! viene passato alla subroutine mediante
!   -   la matrice As ('s' sta per sparsa), che ha dimensioni (N,3) se
!       A ha dimensioni (N,N) ed costituita dalle sole diagonali A
!   -   il termine noto b
!
! Gli elementi di As corrispondono a quelli non nulli di A, nel modo seguente
!
!    As(1,0)  As(1,+1)     0        0        0      0    ...   0           0       As(1,-1)
!    As(2,-1) As(2,0)  As(2,+1)     0        0      0    ...   0           0          0
!       0     As(3,-1) As(3,0)  As(3,+1)     0      0    ...   0           0          0
!       0       0      As(3,-1) As(3,0)  As(3,+1)   0    ...   0           0          0
!       :       :
!       :       :
!       0       :                                         As(N-1,-1)   As(N-1,0)   As(N-1,+1)
!    As(N,+1)   0                                                      As(N,-1)    As(N,0)
!
! Si osserva esplicitamente che, le subdiagonali hanno un numero N
! di elementi, al pari della diagonale principale.



REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)  :: As ! deferred-shape array (so that the bounds are automatically passed)
REAL, DIMENSION(:),                INTENT(IN)  :: b
REAL, DIMENSION(:),                INTENT(OUT) :: x

INTEGER                                 :: N
REAL, DIMENSION(:,:), ALLOCATABLE :: Ac
REAL, DIMENSION(:),   ALLOCATABLE :: uvw_BC, delta, xc1, xc2
!REAL                                          :: xN

! -----------------Allocazioni---------------------------------------------
N = SIZE(As,1)

ALLOCATE(Ac(N-1,-1:1))
ALLOCATE(uvw_BC(N-1),delta(N-1),xc1(N-1),xc2(N-1))
delta = 0

! -----------------Riduzione del sistema-----------------------------------
Ac           = As(1:N-1,:) ! Scarto l'ultima equazione...
uvw_BC         = b(1:N-1)    ! ...e coerentemente anche il termine noto
Ac(1,-1)   = 0           ! Elimino i termini che moltiplicano...
Ac(N-1,+1) = 0           ! ...l'ultima incognita (considerata un parametro)
delta(1)   = As(1,-1)    ! Costruisto il vettore che moltiplica...
delta(N-1) = As(N-1,+1)  ! ...l'ultima incognita

! -----------------Soluzione dei due sistemi-------------------------------
! La soluzione del sistema ridotto è combinazione lineare delle soluzioni
! xc1, relativa al termine noto uvw_BC,
CALL Thomas(Ac,uvw_BC,xc1)
! e xc2, relativa al termine noto delta,
CALL Thomas(Ac,delta,xc2)
! tramite coefficienti 1 e -x(N), ancora incognita. Troviamo allora x(N)
! dall'equazione prima scartata.
x(N) = (b(N) - As(N,+1)*xc1(1) - As(N,-1)*xc1(N-1))/(As(N,0) - As(N,+1)*xc2(1) - As(N,-1)*xc2(N-1))
! Calcoliamo la soluzione finale
x(1:N-1) = xc1 - x(N)*xc2

END SUBROUTINE ThomasC


SUBROUTINE pentdag(a,b,c,d,e,f,u,n)
      IMPLICIT NONE
      SAVE

! solves for a vector u of length n in the pentadiagonal linear system
!  a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i
! input are the a, b, c, d, e, and f and they are not modified

! in its clearest incarnation, this algorithm uses three storage arrays
! called p, q and r. here, the solution vector u is used for r, cutting
! the extra storage down to two arrays.

! declare the pass
      INTEGER :: n
      REAL :: a(n),b(n),c(n),d(n),e(n),f(n),u(n)

! local variables
      INTEGER :: i
      INTEGER, PARAMETER :: nmax = 500
      REAL :: p(nmax),q(nmax),bet,den


! initialize elimination and backsubstitution arrays
      IF (c(1) .EQ. 0.0)  STOP 'eliminate u2 trivially'
      bet  = 1.0d0/c(1)
      p(1) = -d(1) * bet
      q(1) = -e(1) * bet
      u(1) = f(1)  * bet

      bet = c(2) + b(2)*p(1)
      IF (bet .EQ. 0.0) STOP 'singular 1 in pentdag'
      bet = -1.0d0/bet
      p(2) = (d(2) + b(2)*q(1)) * bet
      q(2) = e(2) * bet
      u(2) = (b(2)*u(1) - f(2)) * bet


! reduce to upper triangular
      DO i=3,n
       bet = b(i) + a(i) * p(i-2)
       den = c(i) + a(i)*q(i-2) + bet*p(i-1)
       IF (den .EQ. 0.0) STOP 'singular 2 in pentdag'
       den = -1.0d0/den
       p(i) = (d(i) + bet*q(i-1)) * den
       q(i) = e(i) * den
       u(i) = (a(i)*u(i-2) + bet*u(i-1) - f(i)) * den
      END DO

! backsubstitution
      u(n-1) = u(n-1) + p(n-1) * u(n)
      DO i=n-2,1,-1
       u(i) = u(i) + p(i) * u(i+1) + q(i) * u(i+2)
      END DO
      RETURN

END SUBROUTINE pentdag


SUBROUTINE cypent(a,b,c,d,e,f,cp1,cp2,cp3,cp4,cp5,cp6,x,n)
      IMPLICIT NONE
      SAVE

! solves for x(1:n) a pentadiagonal matrix with nonzero entries in the
! lower left and upper right corners of the matrix:

!   x    x    x    0    0  ... 0  cp1  cp2
!   x    x    x    x    0  ... 0  0    cp3
!   .    .    .    .    .  .   .  .    .
!   cp4  0    0    ...  0  0   x  x    x
!   cp5  cp6  0    ...  0  x   x  x    x

! the woodbury formula is applied to the pentadiagonal matrix.


! declare the pass
      INTEGER :: n
      REAL    :: a(n),b(n),c(n),d(n),e(n),f(n),x(n)
      REAL    :: cp1,cp2,cp3,cp4,cp5,cp6


! local variables
      INTEGER :: i,j,k
      INTEGER, PARAMETER :: nmax = 500
      INTEGER :: indx(nmax)
      REAL    :: u(nmax,4),v(nmax,4),z(nmax,4)
      REAL    :: r(nmax),s(nmax),y(nmax)
      REAL    :: h(4,4),p(4,4),sum


! popular formats
 !01   FORMAT(1x,1p7e10.2)



! initialize
      IF (n .LE. 2) STOP 'n < 2 in routine cypent'
      IF (n .GT. nmax) STOP 'n > nmax in routine cypent'

      DO j=1,4
       DO i=1,n
        u(i,j) = 0.0d0
        v(i,j) = 0.0d0
        z(i,j) = 0.0d0
       END DO
      END DO

      u(1,1)   = 1.0d0
      u(2,2)   = 1.0d0
      u(n-1,3) = 1.0d0
      u(n,4)   = 1.0d0

      v(n-1,1) = cp1
      v(n,1)   = cp2
      v(n,2)   = cp3
      v(1,3)   = cp4
      v(1,4)   = cp5
      v(2,4)   = cp6


! solve the auxillary systems
! recipies equation 2.7.17 and 2.7.20
      CALL pentdag(a,b,c,d,e,u(1,1),z(1,1),n)
      CALL pentdag(a,b,c,d,e,u(1,2),z(1,2),n)
      CALL pentdag(a,b,c,d,e,u(1,3),z(1,3),n)
      CALL pentdag(a,b,c,d,e,u(1,4),z(1,4),n)
      CALL pentdag(a,b,c,d,e,f,y,n)


! form the 4x4 matrix h
! recipies equation 2.7.19

      DO j=1,4
       DO i=1,4
        sum = 0.0d0
        DO k=1,n
         sum = sum + v(k,j) * z(k,i)
        end DO
        p(j,i) = sum
       end DO
      END DO
      DO i=1,4
       p(i,i) = p(i,i) + 1.0d0
      END DO
      CALL luinv(p,4,4,indx,h)

! form the solution
! recipe equation 2.7.21

      DO j=1,4
       r(j) = 0.0d0
       DO k=1,n
        r(j) = r(j) + v(k,j) * y(k)
       END DO
      END DO

      DO j=1,4
       s(j) = 0.0d0
       DO k=1,4
        s(j) = s(j) + h(j,k) * r(k)
       END DO
      END DO

      DO j=1,n
       sum = 0.0d0
       DO k=1,4
        sum = sum + z(j,k) * s(k)
       END DO
       x(j) = y(j) - sum
      END DO

      RETURN
END SUBROUTINE cypent



SUBROUTINE luinv(a,n,np,indx,y)
      IMPLICIT NONE
      SAVE

! this routine takes as input the n by n matrix a, of physical dimension
! np by np and on output fills y with the inverse of a
!
! declare
      INTEGER :: n,np,i,j,indx(np)
      REAL    :: a(np,np),y(np,np),d

! set y to the identity matrix
      DO j=1,n
       DO i=1,n
        y(i,j) = 0.0d0
       END DO
      END DO
      DO i=1,n
       y(i,i) = 1.0d0
      END DO

! decomp and backsubstitute each column
      CALL ludcmp(a,n,np,indx,d)
      DO j=1,n
       CALL lubksb(a,n,np,indx,y(1,j))
      END DO
      RETURN
END SUBROUTINE luinv



SUBROUTINE lubksb(a,n,np,indx,b)
      IMPLICIT NONE
      SAVE

! solves a set of n linear equations ax=b. a is input in its lu decomposition
! form, determined by the routine above ludcmp. indx is input as the
! permutation vector also returned by ludcmp. b is input as the right hand
! side vector and returns with the solution vector x.
! a,n ans np are not modified by this routine and thus can be left in place
! for successive calls (i.e matrix inversion)
!
! declare
      INTEGER :: n,np,indx(np),i,ii,j,ll
      REAL    :: a(np,np),b(np),sum

! when ii is > 0, ii becomes the index of the first nonzero element of b
! this is forward substitution of equation 2.3.6, and unscamble in place
      ii = 0
      DO i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       IF (ii .NE. 0) THEN
        DO j=ii,i-1
         sum = sum - a(i,j) * b(j)
        ENDDO

! nonzero element was found, so dos the sums in the loop above
       ELSE IF (sum .NE. 0.0) THEN
        ii  = i
       END IF
       b(i) = sum
      ENDDO

! back substitution equation 2.3.7
      DO i = n,1,-1
       sum = b(i)
       IF (i .LT. n) THEN
        DO j=i+1,n
         sum = sum - a(i,j) * b(j)
        ENDDO
       END IF
       b(i) = sum/a(i,i)
      ENDDO
      RETURN
END SUBROUTINE lubksb


subroutine ludcmp(a,n,np,indx,d)
      implicit none
      save

! given the matrix a(n,n), with physical dimsnsions a(np,ap) this routine
! replaces a by the lu decompostion of a row-wise permutation of itself.
! input are a,n,np. output is a, indx which records the row
! permutations effected by the partial pivoting, and d which is 1 if
! the number of interchanges is even, -1 if odd.
! use routine lubksb to solve a system of linear equations.
!
! nmax is the largest expected value of n
!
! declare
      INTEGER          n,np,indx(np),i,j,k,imax
      INTEGER, PARAMETER :: nmax = 500
      REAL :: a(np,np),d,vv(nmax),aamax,sum,dum
      REAL, PARAMETER :: tiny=1.0!d-20


! vv stores the implicit scaling of each row
! loop over the rows to get the scaling information
      d = 1.0d0
      DO i=1,n
       aamax = 0.0d0
       DO j=1,n
        IF (abs(a(i,j)) .GT. aamax) aamax = abs(a(i,j))
       ENDDO
       IF (aamax .EQ. 0.0) STOP 'singular matrix in ludcmp'
       vv(i) = 1.0d0/aamax
      ENDDO

! for each column apply crouts method; see equation 2.3.12
      DO j=1,n
       DO i=1,j-1
        sum = a(i,j)
        DO k=1,i-1
         sum = sum - a(i,k)*a(k,j)
        ENDDO
        a(i,j) = sum
       ENDDO

! find the largest pivot element
       aamax = 0.0d0
       DO i=j,n
        sum = a(i,j)
        DO k=1,j-1
         sum = sum - a(i,k)*a(k,j)
        ENDDO
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        IF (dum .GE. aamax) THEN
         imax  = i
         aamax = dum
        END IF
       ENDDO

! IF we need to interchange rows
       IF (j .NE. imax) THEN
        DO k=1,n
         dum       = a(imax,k)
         a(imax,k) = a(j,k)
         a(j,k)    = dum
        ENDDO
        d          = -d
        vv(imax)   = vv(j)
       END IF

! divide by the pivot element
       indx(j) = imax
       IF (a(j,j) .EQ. 0.0) a(j,j) = tiny
       IF (j .NE. n) THEN
        dum = 1.0d0/a(j,j)
        DO i=j+1,n
         a(i,j) = a(i,j)*dum
        ENDDO
       END IF

! and go back for another column of crouts method
      ENDDO
      RETURN
END SUBROUTINE ludcmp



END MODULE Thomas_suite
