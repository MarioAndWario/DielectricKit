#ifdef HAS_GVEC
#define is_greater(ii,jj) \
  (AA(ii)-AA(jj)>TOL.or.(AA(ii)-AA(jj)>-TOL.and.GK(ii)>GK(jj)))
#else
#define is_greater(ii,jj) \
  (AA(ii)>AA(jj))
#endif

  ! This file is based on the work of Michel Olagnon.
  ! The original code for the MRGRNK subroutine is available at:
  ! http://fortran-2000.com/rank/mrgrnk.f90

  ! MRGRNK - Copyright (c) Michel Olagnon
  ! Copying and distribution of this file, with or without modification,
  ! are permitted in any medium without royalty provided the copyright
  ! notice and this notice are preserved.  This file is offered as-is,
  ! without any warranty.


  ! FHJ: WARNING - make sure you don`t change the following lines too much,
  ! otherwise they will be longer than 120 characters after the preprocessors kicks in.
  ! Note that, if there the extra "gvec" argument, we use a tolerance to figure
  ! out if the two items AA(ii) and AA(jj) are degenerate.

  ! FHJ: Note - we need to nest the "JOIN" macro otherwise the symbol LABEL
  ! doesn`t get expanded by the C preprocessor.
#define JOIN2(x,y) x ## _ ## y
#define JOIN(x,y) JOIN2(x,y)
#define LABEL_INSERTSORT JOIN(LABEL,insertsort)

  !> Sorts (actually, ranks) a small array AA using the insert sort algorithm.
  subroutine LABEL_INSERTSORT(NVAL, AA, ord &
#ifdef HAS_GVEC
    , GK &
#endif
    )
    integer, intent(in) :: NVAL
    DTYPE, intent(in) :: AA(NVAL)
    integer, intent(inout) :: ord(NVAL)
#ifdef HAS_GVEC
    integer, intent(in) :: GK(NVAL)
    DTYPE, parameter :: TOL=TOL_ZERO
#endif

    integer :: ii, jj, tord

    PUSH_SUB(LABEL_INSERTSORT)

    do ii = 2, NVAL
       tord = ord(ii)
       jj = ii - 1
       do while (jj>0)
          if (.not.is_greater(ord(jj),tord)) exit
          ord(jj+1) = ord(jj)
          jj = jj - 1
       enddo
       ord(jj+1) = tord
    enddo

    POP_SUB(LABEL_INSERTSORT)

  end subroutine LABEL_INSERTSORT

  !> Sorts (actually, ranks) a real/integer array AA.
  !! The rank is written to the output array ord.
  !! This subroutine is based on the routine MRGRNK by Michel Olagnon, which
  !! uses the merge sort algorithm.
  subroutine LABEL(NVAL, AA, ord &
#ifdef HAS_GVEC
    , gvec &
#endif
    )
    integer, intent(in) :: NVAL
    DTYPE, intent(in) :: AA(NVAL)
    integer, intent(out) :: ord(NVAL)
#ifdef HAS_GVEC
    integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data

    integer :: GK(NVAL), fftgrid(3)
    DTYPE, parameter :: TOL=TOL_ZERO
#endif
    !
    integer :: JT(NVAL)
    integer :: LMTNA, LMTNC, IRNG1, IRNG2
    integer :: IIND, ID, IWRK, IWRKF, JINDA, IA, IB

    PUSH_SUB(LABEL)

#ifdef HAS_GVEC
    fftgrid(1:3) = maxval(gvec(1:3,1:NVAL), 2) - minval(gvec(1:3,1:NVAL), 2) + 1
    do IIND=1,NVAL
       GK(IIND) = gvec(3,IIND) + fftgrid(3)*(gvec(2,IIND) + fftgrid(2)*gvec(1,IIND))
    enddo
#endif
    !
    !  Fill-in the index array, creating ordered couples
    !
    Do IIND = 2, NVAL, 2
       If (&
            is_greater(IIND,IIND-1)&
            ) Then
          ord (IIND-1) = IIND - 1
          ord (IIND) = IIND
       Else
          ord (IIND-1) = IIND
          ord (IIND) = IIND - 1
       End If
    End Do
    If (Modulo(NVAL, 2) /= 0) Then
       ord (NVAL) = NVAL
    End If

    ! FHJ - shortcut if the array is small enough
    if (NVAL<16) then
       call LABEL_INSERTSORT(NVAL, AA, ord &
#ifdef HAS_GVEC
       , GK &
#endif
       )
       return
    endif

    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into    C  -  C  - ...
    !
    LMTNA = 2
    LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    Do
       If (NVAL <= 2) Exit
       !
       !  Loop on merges of A and B into C
       !
       Do ID = 0, NVAL - 1, 4
          If ((ID+4) > NVAL) Then
             If ((ID+2) >= NVAL) Exit
             !
             !  1 2 3
             !
             If (&
                  is_greater(ord(ID+3),ord(ID+2))&
                  ) Exit
             !
             !  1 3 2
             !
             If (&
                  is_greater(ord(ID+3),ord(ID+1))&
                  ) Then
                IRNG2 = ord (ID+2)
                ord (ID+2) = ord (ID+3)
                ord (ID+3) = IRNG2
                !
                !  3 1 2
                !
             Else
                IRNG1 = ord (ID+1)
                ord (ID+1) = ord (ID+3)
                ord (ID+3) = ord (ID+2)
                ord (ID+2) = IRNG1
             End If
             Exit
          End If
          !
          !  1 2 3 4
          !
          If (&
               is_greater(ord(ID+3),ord(ID+2))&
               ) Cycle
          !
          !  1 3 x x
          !
          If (&
               is_greater(ord(ID+3),ord(ID+1))&
               ) Then
             IRNG2 = ord (ID+2)
             ord (ID+2) = ord (ID+3)
             If (&
                  is_greater(ord(ID+4),IRNG2)&
                  ) Then
                !  1 3 2 4
                ord (ID+3) = IRNG2
             Else
                !  1 3 4 2
                ord (ID+3) = ord (ID+4)
                ord (ID+4) = IRNG2
             End If
             !
             !  3 x x x
             !
          Else
             IRNG1 = ord (ID+1)
             IRNG2 = ord (ID+2)
             ord (ID+1) = ord (ID+3)
             If (&
                  is_greater(ord(ID+4),IRNG1)&
                  ) Then
                ord (ID+2) = IRNG1
                If (&
                     is_greater(ord(ID+4),IRNG2)&
                     ) Then
                   !  3 1 2 4
                   ord (ID+3) = IRNG2
                Else
                   !  3 1 4 2
                   ord (ID+3) = ord (ID+4)
                   ord (ID+4) = IRNG2
                End If
             Else
                !  3 4 1 2
                ord (ID+2) = ord (ID+4)
                ord (ID+3) = IRNG1
                ord (ID+4) = IRNG2
             End If
          End If
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 4
       Exit
    End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    Do
       If (LMTNA >= NVAL) Exit
       IWRKF = 0
       LMTNC = 2 * LMTNC
       !
       !  Loop on merges of A and B into C
       !
       Do
          IWRK = IWRKF
          ID = IWRKF + 1
          JINDA = IWRKF + LMTNA
          IWRKF = IWRKF + LMTNC
          If (IWRKF >= NVAL) Then
             If (JINDA >= NVAL) Exit
             IWRKF = NVAL
          End If
          IA = 1
          IB = JINDA + 1
          !
          !  Shortcut for the case when the max of A is smaller
          !  than the min of B. This line may be activated when the
          !  initial set is already close to sorted.
          !
          IF (&
               is_greater(ord(IB),ord(JINDA))&
               ) CYCLE
          !
          !  One steps in the C subset, that we build in the final rank array
          !
          !  Make a copy of the rank array for the merge iteration
          !
          JT (1:LMTNA) = ord (ID:JINDA)
          !
          Do
             IWRK = IWRK + 1
             !
             !  We still have unprocessed values in both A and B
             !
             If (&
                  is_greater(JT(IA),ord(IB))&
                  ) Then
                ord (IWRK) = ord (IB)
                IB = IB + 1
                If (IB > IWRKF) Then
                   !  Only A still with unprocessed values
                   ord (IWRK+1:IWRKF) = JT (IA:LMTNA)
                   Exit
                End If
             Else
                ord (IWRK) = JT (IA)
                IA = IA + 1
                If (IA > LMTNA) Exit! Only B still with unprocessed values
             End If
             !
          End Do
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 2 * LMTNA
    End Do
    !
    POP_SUB(LABEL)
    !
  End Subroutine LABEL

#undef is_greater
