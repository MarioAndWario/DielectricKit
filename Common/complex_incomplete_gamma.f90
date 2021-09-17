MODULE complex_incomplete_gamma_m
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2002-03-31  Time: 00:10:14

  ! --- Written By Eric Kostlan & Dmitry Gokhman
  ! --- March  1986
  ! --- For documentation, see:
  ! http://www.math.utsa.edu/~gokhman/papers/igf.html

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
  PRIVATE
  PUBLIC  :: cdig

CONTAINS

  FUNCTION cdig(alpha, x) RESULT(fn_val)

    COMPLEX (dp), INTENT(IN)  :: alpha
    COMPLEX (dp), INTENT(IN)  :: x
    COMPLEX (dp)              :: fn_val

    COMPLEX (dp)  :: p, q
    INTEGER       :: i, ilim
    REAL (dp), PARAMETER     :: zero = 0.0_dp, xlim = 1.0_dp
    COMPLEX (dp), PARAMETER  :: cone = (1.0_dp, 0.0_dp)
    REAL (dp), PARAMETER     :: re = (0.36787944117144232_dp, 0.0_dp)
    INTEGER, PARAMETER       :: ibuf = 34

    ! --- If x is near the negative real axis, then shift to x=1.
    IF (dnrm(x) < xlim .OR. REAL(x, KIND=dp) < zero .AND. ABS(AIMAG(x)) < xlim) THEN
       fn_val = re / cdh(alpha, cone)
       ilim = REAL(x/re, KIND=dp)
       DO  i = 0, ibuf - ilim
          CALL term(alpha, x, i, p, q)
          fn_val = fn_val + p * q
       END DO
    ELSE
       fn_val = EXP(-x + alpha*LOG(x)) / cdh(alpha, x)
    END IF
    RETURN
  END FUNCTION cdig

  FUNCTION cdh(alpha, x) RESULT(fn_val)
    ! --- Written By Eric Kostlan & Dmitry Gokhman
    ! --- March  1986

    COMPLEX (dp), INTENT(IN)  :: alpha
    COMPLEX (dp), INTENT(IN)  :: x
    COMPLEX (dp)              :: fn_val

    COMPLEX (dp)  :: term, sum, cn, alpha1
    INTEGER       :: i, n
    REAL (dp), PARAMETER  :: one = 1.0_dp

    ! --- If Re(alpha-x) is too big, shift alpha.
    n = REAL(alpha-x, KIND=dp)
    IF (n > 0) THEN
       cn = n
       alpha1 = alpha - cn
       term = one / x
       sum = term
       DO  i = 1, n - 1
          cn = n - i
          term = term * (alpha1 + cn) / x
          sum = term + sum
       END DO
       sum = sum + term * alpha1 / cdhs(alpha1, x)
       fn_val = one / sum
    ELSE
       fn_val = cdhs(alpha, x)
    END IF
    RETURN
  END FUNCTION cdh

  FUNCTION cdhs(alpha, x) RESULT(fn_val)
    ! --- Written By Eric Kostlan & Dmitry Gokhman
    ! --- March  1986

    COMPLEX (dp), INTENT(IN)  :: alpha
    COMPLEX (dp), INTENT(IN)  :: x
    COMPLEX (dp)              :: fn_val

    COMPLEX (dp)  :: p0, q0, p1, q1, r0, r1, ci, factor
    INTEGER       :: i
    REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp
    REAL (dp), PARAMETER  :: tol1 = 1.0D+10, tol2 = 1.0D-10, error = 5.D-18
    INTEGER, PARAMETER    :: ilim = 100000

    q0 = one
    q1 = one
    p0 = x
    p1 = x + one - alpha
    DO  i = 1, ilim
       ci = i
       IF (p0 /= zero .AND. q0 /= zero .AND. q1 /= zero) THEN
          r0 = p0 / q0
          r1 = p1 / q1
          IF (dnrm(r0-r1) <= dnrm(r1)*error) THEN
             fn_val = r1
             RETURN
          END IF
          ! --------- Occasionally renormalize the sequences to avoid over(under)flow.
          IF (dnrm(p0) > tol1 .OR. dnrm(p0) < tol2 .OR. dnrm(q0) > tol1  &
               .OR. dnrm(q0) < tol2) THEN
             factor = p0 * q0
             p0 = p0 / factor
             q0 = q0 / factor
             p1 = p1 / factor
             q1 = q1 / factor
          END IF
       END IF
       p0 = x * p1 + ci * p0
       q0 = x * q1 + ci * q0
       p1 = p0 + (ci+one-alpha) * p1
       q1 = q0 + (ci+one-alpha) * q1
    END DO
    ! --- If the peripheral routines are written correctly,
    ! --- the following four statements should never be executed.
    WRITE(*, *) 'cdhs:  *** Warning: i >', ilim
    WRITE(*, *) 'cdhs:  *** r0,r1= ', r0, r1
    fn_val = half * (r0+r1)
    RETURN
  END FUNCTION cdhs

  SUBROUTINE term(alpha, x, i, p, q)
    ! --- Calculate p*q = -1**i(1 - x**(alpha+i))/(alpha+i)i ! carefully.

    COMPLEX (dp), INTENT(IN)   :: alpha
    COMPLEX (dp), INTENT(IN)   :: x
    INTEGER, INTENT(IN)        :: i
    COMPLEX (dp), INTENT(OUT)  :: p
    COMPLEX (dp), INTENT(OUT)  :: q

    COMPLEX (dp)  :: ci, cdlx, alphai
    REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
    REAL (dp), PARAMETER  :: tol = 3.0D-7, xlim = 39.0_dp

    IF (i == 0) q = one
    ci = i
    alphai = alpha + ci
    IF (x == zero) THEN
       p = one / alphai
       IF (i /= 0) q = -q / ci
       RETURN
    END IF
    cdlx = LOG(x)

    ! --- If (1 - x**alphai) = -x**alphai on the computer,
    ! --- then change the inductive scheme to avoid overflow.
    IF (REAL(alphai*cdlx, KIND=dp) > xlim .AND. i /= 0) THEN
       p = p * (alphai - one) / alphai
       q = -q * x / ci
       RETURN
    END IF
    IF (dnrm(alphai) > tol) THEN
       p = (one - x**alphai) / alphai
    ELSE
       p = -cdlx * (one + cdlx*alphai/two)
    END IF
    IF (i /= 0) q = -q / ci
    RETURN
  END SUBROUTINE term

  FUNCTION dnrm(z) RESULT(fn_val)
    COMPLEX (dp), INTENT(IN)  :: z
    REAL (dp)                 :: fn_val

    fn_val = ABS(REAL(z, KIND=dp)) + ABS(AIMAG(z))
    RETURN
  END FUNCTION dnrm

END MODULE complex_incomplete_gamma_m
