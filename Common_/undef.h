! use between inclusions of f_defs.h in template modules
! list here everything defined differently by flavor in f_defs.h
! these undefs prevent lots of warnings from cpp
#undef MYFLAVOR
#undef SCALAR
#undef MPI_SCALAR
#undef MYCONJG
#undef MYREAL
#undef X
#undef pX
#undef ZERO
#undef ONE
#undef SCALARIFY
#undef SCALARIFY2
#undef COMPLEXIFY
#undef SCALARSIZE
#undef ONLYIFCPLX
#undef MYABS2

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
