AC_PREREQ([2.69])
AC_INIT([saber], [0.1.1], [raphael.mabit@gmail.com])
AC_CONFIG_SRCDIR([src/saber_shim_stan.cpp])
AC_CONFIG_HEADERS([config.h])

# find a C compiler
AC_PROG_CC

# find pkg-config
PKG_PROG_PKG_CONFIG

# Check for saber-lib ≥ 0.1.0
# – on success: stash flags into SABER_CFLAGS/SABER_LIBS
# – on failure: abort with a clear error
PKG_CHECK_MODULES([SABER], [saber-lib >= 0.1.0],
  [
    AC_MSG_NOTICE([found saber-lib via pkg-config])
    AC_SUBST([SABER_CFLAGS])
    AC_SUBST([SABER_LIBS])
  ],
  [
    AC_MSG_ERROR([
*** ERROR: pkg-config could not find saber-lib ≥ 0.1.0
    Please install saber-lib and make sure its .pc file is on PKG_CONFIG_PATH.
])
  ]
)

# Check for STAN
AC_MSG_CHECKING([for CmdStan directory])
AC_ARG_WITH([cmdstan],
            AS_HELP_STRING([--with-cmdstan=DIR], [path to CmdStan root]))
AC_ARG_VAR([CMDSTAN], [Path to CmdStan root directory])

# 1) user‑supplied value wins
if test -n "$with_cmdstan"; then
  CMDSTAN="$with_cmdstan"
fi

# 2) scan a few typical install roots
FALLBACK=""

# User’s home?
if test -n "$HOME" && test -d "$HOME"; then
  FALLBACK="$FALLBACK $HOME/cmdstan $HOME/.cmdstan"
else
  AC_MSG_ERROR([HOME is not set or not a directory – skipping home‑relative guesses])
fi

# System‑wide guesses (always safe to try)
FALLBACK="$FALLBACK /opt/cmdstan /usr/local/cmdstan"

# Iterate with per‑path feedback
for d in $CMDSTAN $FALLBACK ; do
  test -z "$d" && continue
  AC_MSG_CHECKING([whether CmdStan is in $d])
  if test -d "$d"; then
    CMDSTAN=$d
    AC_MSG_RESULT([yes])
    break
  else
    AC_MSG_RESULT([no])
  fi
done

# 3) final sanity checks
if test -z "$CMDSTAN" || test ! -d "$CMDSTAN"; then
  AC_MSG_ERROR([*** CmdStan not found.  Set CMDSTAN or use --with-cmdstan=... ])
fi
AC_MSG_RESULT([found CmdStan in $CMDSTAN])

STAN_MATH_INCDIR="$CMDSTAN/stan/lib/stan_math"
test -d "$STAN_MATH_INCDIR" || AC_MSG_ERROR([Stan‑Math headers missing])

# Find Eigen bundeld with Stan-Math
EIGEN_DIR=`ls -1d "$STAN_MATH_INCDIR"/lib/eigen_* 2>/dev/null | head -n 1`
if test -z "$EIGEN_DIR" || test ! -d "$EIGEN_DIR"; then
  AC_MSG_ERROR([bundled Eigen headers not found under $STAN_MATH_INCDIR/lib])
fi

# Find Boost
BOOST_DIR=`ls -1d "$STAN_MATH_INCDIR"/lib/boost_* 2>/dev/null | head -n 1`
if test -z "$BOOST_DIR" || test ! -d "$BOOST_DIR"; then
  AC_MSG_ERROR([bundled Eigen headers not found under $STAN_MATH_INCDIR/lib])
fi

# Find tbb
TBB_ROOT=`ls -1d "$STAN_MATH_INCDIR"/lib/tbb_* 2>/dev/null | head -n 1`
if test -z "$TBB_ROOT" || test ! -d "$TBB_ROOT/include/tbb"; then
  AC_MSG_ERROR([bundled TBB headers not found under $STAN_MATH_INCDIR/lib])
fi
TBB_INC="$TBB_ROOT/include"

TBB_DIR=`ls -1d "$STAN_MATH_INCDIR"/lib/tbb 2>/dev/null | head -n 1`
test -f "$TBB_DIR/libtbb.so" || AC_MSG_ERROR([libtbb.so not found in $TBB_DIR])
AC_SUBST([STAN_TBB_LIBS], ["-L$TBB_DIR -ltbb"])

# Find sundials
SUNDIALS_ROOT=`ls -1d "$STAN_MATH_INCDIR"/lib/sundials_* 2>/dev/null | head -n 1`
if test -z "$SUNDIALS_ROOT" || test ! -d "$SUNDIALS_ROOT/include/sundials"; then
  AC_MSG_ERROR([bundled sundials headers not found under $STAN_MATH_INCDIR/lib])
fi
SUNDIALS_INC="$SUNDIALS_ROOT/include"

AC_SUBST([STAN_CPPFLAGS],
         ["-I$STAN_MATH_INCDIR \
           -isystem $EIGEN_DIR \
           -isystem $BOOST_DIR \
           -isystem $TBB_INC \
           -isystem $SUNDIALS_INC"])

# Threads
AC_SEARCH_LIBS([pthread_create], [pthread])
AC_MSG_CHECKING([for compiler flag to enable POSIX threads])
AX_PTHREAD
AC_MSG_RESULT([$PTHREAD_CFLAGS])

# Export both compile‑time and link‑time flags
AC_SUBST([STAN_PTHREAD_CFLAGS], [$PTHREAD_CFLAGS])
AC_SUBST([STAN_PTHREAD_LIBS],   [$PTHREAD_LIBS])

# Generate src/Makevars from src/Makevars.in
AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
