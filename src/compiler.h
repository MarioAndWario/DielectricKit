! Created Sept 2011 by DAS. 
! Define characteristics of various compilers, via compiler symbols (e.g. -DGNU) 
! to be used directly from the arch.mk files, and then defining what we need to do 
! for that compiler via the symbols for various properties (e.g. NOSIZEOF). 
! Ideally, to support a new compiler, one need only change this file, adding a 
! new block to define what -DNEWCOMPILER would mean. 
! NOTE: of course, Makefile-level issues still need to be handled in common-rules.mk 

#ifdef PGI
#  define MCLOCKINT
#  define CENTISEC
#  define HOSTNAMEINT
#  define FTELLINT
#  define HAS_SLEEP
#  define COMPILER_STR "PGI"
#endif

#ifdef INTEL
#  define SIZEOF64
#  define HOSTNAMENOE
#  define HOSTNAMEMOD
#  define MILLISEC
#  define FTELLINT
#  define HAS_SLEEP
#  define COMPILER_STR "INTEL"
#endif

! very ancient version may require NOSIZEOF 
#ifdef GNU
#  define SIZEOF64
#  define MICROSEC
#  define NOIARGCINT
#  define FTELLROUTINE
#  define HAS_SLEEP
#  define COMPILER_STR "GNU"
#endif

#ifdef G95
#  define CPUTIME
#  define NOFSEEK
#  define NOOPENMP
#  define HAS_SLEEP
#  define COMPILER_STR "G95"
#endif

! open64 is very similar to path, it is an open-sourced version of it
! omp_lib.f90 needed to do OpenMP, see common-rules.mk.
#ifdef OPEN64
#  define NOSIZEOF
#  define CALLFLUSH
#  define HOSTNAMEINT
#  define NOFSEEK ! compiles, but does not seem to work
! #  define FTELLROUTINE
#  define NOOPENMP
#  define HAS_SLEEP
#  define COMPILER_STR "OPEN64"
#endif

! path before 4.0.9 lacks SIZEOF
#ifdef PATH
#  define CALLFLUSH
#  define HOSTNAMEINT
#  define NOFSEEK ! compiles, but does not seem to work
! #  define FTELLROUTINE
#  define HAS_SLEEP
#  define COMPILER_STR "PATH"
#endif

! both open64 and path die on fseek with:
!lib-5002 : UNRECOVERABLE library error 
!  This FFIO request is not supported.
!
!Encountered during a GETPOS on unit 8

#ifdef XLF
#  define NOFLUSH
#  define HOSTNAMEUNDERSCORE
#  define MCLOCKINT
#  define CENTISEC
#  define HOSTNAMEINT
#  define NOFSEEK ! compiles, but does not seem to work
! #  define FTELLINT
#  define COMPILER_STR "XLF" 
#endif 
 
#ifdef SUN 
#  define NOSIZEOF 
#  define CPUTIME 
#  define HOSTNAMEINT 
#  define NOFSEEK ! compiles, but does not seem to work
! #  define FTELLINT
#  define HAS_SLEEP 
#  define COMPILER_STR "SUN" 
#endif 
 
#ifdef NAG 
#  define NOSIZEOF 
#  define CPUTIME 
#  define NOOPENMP 
#  define GETHOSTNAME 
#  define NOIARGCINT 
#  define NOFSEEK 
#  define SYSTEMMOD_NAG 
! sleep is defined in f90_unix_env
#  define HAS_SLEEP 
#  define COMPILER_STR "NAG" 
#endif 
 
#ifdef ABSOFT 
#  define CALLFLUSH 
#  define CPUTIME 
#  define NOIARGCINT 
#  define SYSTEMFUNCTION 
#  define FSEEKFUNCTION 
#  define SYSTEMMOD_AB 
! sleep is defined in unix_library
#  define NOFSEEK ! compiles, but does not seem to work
! #  define FTELLINT sleep is defined in unix_library
#  define HAS_SLEEP
#  define COMPILER_STR "ABSOFT"
#endif

! cce 7.4.4 and before support sizeof for intrinsic types, but need NOSIZEOF_TYPE
! cce 8.0.0 and later do not allow sizeof for multidimensional arrays, requiring us
! to turn sizeof off everywhere. Why would Cray do this?
#ifdef CRAY
#  define HOSTNAMEINT
#  define NOFSEEK
#  define CPUTIME
#  define NOSIZEOF
#  define HAS_SLEEP
#  define COMPILER_STR "CRAY"
#endif

#ifndef COMPILER_STR
#  define COMPILER_STR "UNKNOWN"
#endif

! It is considered a bug in OPEN64 that sizeof will not work in our code. 
#if defined NOSIZEOF
#  define SIZEOF(x) 1
#elif defined SIZEOFTRICK
#  define SIZEOF(x) int(size(x),8)*int(sizeof(reshape(x,(/ 1 /))),8)
#else
#  define SIZEOF(x) sizeof(x)
#endif

#if defined CALLFLUSH
! XLF will not accept this version 
#  define FLUSH(x) call flush(x)
#elif defined NOFLUSH
#  define FLUSH(x)
#else
! Fortran 2003 prefers this as a statement, not an intrinsic 
#  define FLUSH(x) flush(x)
#endif

! on some platforms there is a different return value for sizeof if build is 64-bit 
#if defined SIZEOF64 && defined _LP64 
#  define INTSIZEOF integer*8
#elif defined SIZEOFTRICK
#  define INTSIZEOF integer*8
#else
#  define INTSIZEOF integer
#endif

! name of routine to get name of host program is running on 
#if defined HOSTNAMEUNDERSCORE
#  define HOSTNAME hostnm_
#elif defined HOSTNAMENOE
#  define HOSTNAME hostnam
#else
#  define HOSTNAME hostnm
#endif

#ifdef GETHOSTNAME
#  define HOSTNAMECALL(name, stat) call gethostname(name)
#else
#  define HOSTNAMECALL(name, stat) stat = HOSTNAME(name)
#endif

! module required for HOSTNAME routine to be usable 
#if defined HOSTNAMEMOD
#  define USEHOSTNAMEMOD use ifport, only : hostnam
#else
#  define USEHOSTNAMEMOD
#endif

! note: interfaces are split into two lines because ifort for some reason decrees of
! 'end interface' : "END statement must be only statement on line". 

! HOSTNAMEINT enables interface for HOSTNAME routine in intrinsics_m 
#if defined HOSTNAMEINT
#  define HOSTNAME_INTERFACE \
  interface; \
  integer function HOSTNAME(name); \
    character(len=*), intent(out) :: name; \
  end function HOSTNAME
#  define HOSTNAME_INTERFACE2 end interface
#else
#  define HOSTNAME_INTERFACE
#  define HOSTNAME_INTERFACE2
#endif

! how to get the cpu time in seconds 
#if defined CPUTIME
#  define TIMGET(cpu) call cpu_time(cpu)
#elif defined MICROSEC
#  define TIMGET(cpu) cpu=mclock()*1.0d-6
#elif defined MILLISEC
#  define TIMGET(cpu) cpu=mclock()*1.0d-3
#elif defined CENTISEC
#  define TIMGET(cpu) cpu=mclock()*1.0d-2
#else
#  define TIMGET(cpu) cpu=0.0d0
#endif

! interface required for mclock routine (timing) to be usable 
#if defined MCLOCKINT
#  define MCLOCK_INTERFACE  interface; integer function mclock(); end function mclock
#  define MCLOCK_INTERFACE2 end interface
#else
#  define MCLOCK_INTERFACE
#  define MCLOCK_INTERFACE2
#endif

! interface required for iargc routine (command-line arguments) to be usable 
#ifdef NOIARGCINT
#  define IARGC_INTERFACE
#  define IARGC_INTERFACE2
#else
#  define IARGC_INTERFACE  interface; integer function iargc(); end function iargc
#  define IARGC_INTERFACE2 end interface
#endif

! ftell gives you the current location in a file, to fseek back to it 
#ifdef FTELLINT
#  define FTELL_INTERFACE \
  interface; \
    integer function ftell(unit); \
      integer, intent(in) :: unit; \
    end function ftell
#  define FTELL_INTERFACE2 end interface
#else
#  define FTELL_INTERFACE
#  define FTELL_INTERFACE2
#endif

! if no fseek, ftell is useless 
#ifdef FTELLROUTINE
#  define FTELLCALL(unit, loc) call ftell(unit, loc)
#elif defined NOFSEEK
#  define FTELLCALL(unit, loc)
#else
#  define FTELLCALL(unit, loc) loc = ftell(unit)
#endif

! #warning This compiler does not support fseek. 
! fseek returns to a location in a file, bookmarked by ftell. G95 lacks it 
#ifdef NOFSEEK
#  define FSEEKCALL(unit, loc, whence) call die("This compiler does not support fseek.")
#elif defined FSEEKFUNCTION
#  define FSEEKCALL(unit, loc, whence) alc = fseek(unit, loc, whence)
#else
#  define FSEEKCALL(unit, loc, whence) call fseek(unit, loc, whence)
#endif

! intrinsic module for OpenMP. external for Open64 (see common-rules.mk). NAG and G95 do not support OpenMP 
#ifdef NOOPENMP
#  define USEOMPLIB
#else
#  define USEOMPLIB use omp_lib
#endif

! using a global var here to avoid need for conditional local declaration 
#ifdef NOSYSTEM
#  define SYSTEMCALL(command) write(0,*) 'WARNING: system call not available.'
#elif defined SYSTEMFUNCTION
#  define SYSTEMCALL(command) alc = system(command)
#else
#  define SYSTEMCALL(command) call system(command)
#endif

#ifdef SYSTEMMOD_AB
#  define SYSTEMMOD use unix_library
#elif defined SYSTEMMOD_NAG
#  define SYSTEMMOD use f90_unix_env; use f90_unix_proc
#else
#  define SYSTEMMOD
#endif

#ifdef HAS_SLEEP
#  define MYSLEEP(x) call sleep(x)
#else
#  define MYSLEEP(x)
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
