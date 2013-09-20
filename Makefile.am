# options passed to aclocal, which is a tool for making macroes visible to
# autoconf. We use -I to tell aclocal where we put the local macros.
ACLOCAL_AMFLAGS = -I build/autotools/m4

# Options passed to the C and C++ PreProcessor (CPP) and compiler
AM_CPPFLAGS = -I${top_srcdir}/ -I$(top_srcdir)/src/ $(DEPS_CFLAGS)

# tell Libtool what the name of the library is.
lib_LTLIBRARIES = libmathicgb.la

# libraries that are needed by this library
libmathicgb_la_LIBADD= $(DEPS_LIBS)

# the sources that are built to make libmathicgb. Listing the headers in
# sources ensure that those files are included in distributions.
libmathicgb_la_SOURCES = src/mathicgb/Range.hpp							\
  src/mathicgb/ReducerPack.hpp src/mathicgb/ReducerPack.cpp				\
  src/mathicgb/ClassicGBAlg.cpp src/mathicgb/ClassicGBAlg.hpp			\
  src/mathicgb/MonoLookup.hpp src/mathicgb/MonoLookup.cpp				\
  src/mathicgb/StaticMonoMap.hpp src/mathicgb/SigPolyBasis.cpp			\
  src/mathicgb/SigPolyBasis.hpp src/mathicgb/Basis.cpp					\
  src/mathicgb/Basis.hpp src/mathicgb/io-util.cpp						\
  src/mathicgb/io-util.hpp src/mathicgb/KoszulQueue.hpp					\
  src/mathicgb/ModuleMonoSet.cpp src/mathicgb/ModuleMonoSet.hpp			\
  src/mathicgb/Poly.cpp src/mathicgb/Poly.hpp							\
  src/mathicgb/PolyBasis.cpp src/mathicgb/PolyBasis.hpp					\
  src/mathicgb/PolyHashTable.cpp src/mathicgb/PolyHashTable.hpp			\
  src/mathicgb/PolyRing.cpp src/mathicgb/PolyRing.hpp					\
  src/mathicgb/Reducer.cpp src/mathicgb/Reducer.hpp						\
  src/mathicgb/ReducerDedup.hpp src/mathicgb/ReducerDedup.cpp			\
  src/mathicgb/ReducerHash.hpp src/mathicgb/ReducerHash.cpp				\
  src/mathicgb/ReducerHashPack.hpp src/mathicgb/ReducerHashPack.cpp		\
  src/mathicgb/ReducerHelper.hpp src/mathicgb/ReducerNoDedup.hpp		\
  src/mathicgb/ReducerNoDedup.cpp src/mathicgb/ReducerPackDedup.hpp		\
  src/mathicgb/ReducerPackDedup.cpp src/mathicgb/SignatureGB.cpp		\
  src/mathicgb/SignatureGB.hpp src/mathicgb/SigSPairs.cpp				\
  src/mathicgb/SigSPairs.hpp src/mathicgb/SPairs.cpp					\
  src/mathicgb/SPairs.hpp src/mathicgb/stdinc.h							\
  src/mathicgb/SigSPairQueue.hpp src/mathicgb/SigSPairQueue.cpp			\
  src/mathicgb/SparseMatrix.hpp src/mathicgb/SparseMatrix.cpp			\
  src/mathicgb/QuadMatrixBuilder.hpp									\
  src/mathicgb/QuadMatrixBuilder.cpp src/mathicgb/TypicalReducer.cpp	\
  src/mathicgb/TypicalReducer.hpp src/mathicgb/F4Reducer.hpp			\
  src/mathicgb/F4Reducer.cpp src/mathicgb/F4MatrixBuilder.hpp			\
  src/mathicgb/F4MatrixBuilder.cpp src/mathicgb/QuadMatrix.hpp			\
  src/mathicgb/QuadMatrix.cpp src/mathicgb/F4MatrixReducer.cpp			\
  src/mathicgb/F4MatrixReducer.hpp src/mathicgb/MonomialMap.hpp			\
  src/mathicgb/RawVector.hpp src/mathicgb/Atomic.hpp					\
  src/mathicgb/FixedSizeMonomialMap.hpp src/mathicgb/CFile.hpp			\
  src/mathicgb/CFile.cpp src/mathicgb/LogDomain.hpp						\
  src/mathicgb/LogDomain.cpp src/mathicgb/LogDomainSet.hpp				\
  src/mathicgb/F4MatrixBuilder2.hpp src/mathicgb/F4MatrixBuilder2.cpp	\
  src/mathicgb/LogDomainSet.cpp src/mathicgb/F4ProtoMatrix.hpp			\
  src/mathicgb/F4ProtoMatrix.cpp src/mathicgb/F4MatrixProject.hpp		\
  src/mathicgb/F4MatrixProjection.cpp src/mathicgb/ScopeExit.hpp		\
  src/mathicgb.cpp src/mathicgb.h src/mathicgb/mtbb.hpp					\
  src/mathicgb/PrimeField.hpp src/mathicgb/MonoMonoid.hpp				\
  src/mathicgb/MonoProcessor.hpp src/mathicgb/MonoOrder.hpp				\
  src/mathicgb/Scanner.hpp src/mathicgb/Scanner.cpp						\
  src/mathicgb/Unchar.hpp src/mathicgb/MathicIO.hpp						\
  src/mathicgb/NonCopyable.hpp


# The headers that libmathicgb installs.
mathicgbA_include_HEADERS = src/mathicgb.h
mathicgbA_includedir = $(includedir)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = build/autotools/mathicgb.pc



# When making a distribution file, Automake knows to include all files
# that are necessary to build the project. EXTRA_DIST specifies files
# to include beyond those used in the build process.
EXTRA_DIST = autogen.sh

bin_PROGRAMS = mgb

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
mgb_SOURCES = src/cli/GBMain.cpp src/cli/CommonParams.hpp		\
  src/cli/CommonParams.cpp src/cli/GBAction.hpp src/cli/GBAction.cpp	\
  src/cli/GBCommonParams.hpp src/cli/GBCommonParams.cpp			\
  src/cli/MatrixAction.cpp src/cli/MatrixAction.hpp			\
  src/cli/SigGBAction.hpp src/cli/SigGBAction.cpp			\
  src/cli/HelpAction.hpp src/cli/HelpAction.cpp
mgb_LDADD = $(top_builddir)/libmathicgb.la $(DEPS_LIBS)

# set up tests to run on "make check"
if with_gtest

TESTS=unittest
check_PROGRAMS=$(TESTS)

unittest_LDADD = $(DEPS_LIBS) $(top_builddir)/libmathicgb.la

test_LIBS=
unittest_SOURCES=src/test/Range.cpp src/test/gtestInclude.cpp			\
  src/test/testMain.cpp src/test/gb-test.cpp src/test/ideals.cpp		\
  src/test/poly-test.cpp src/test/ideals.hpp src/test/SparseMatrix.cpp	\
  src/test/QuadMatrixBuilder.cpp src/test/F4MatrixBuilder.cpp			\
  src/test/F4MatrixReducer.cpp src/test/mathicgb.cpp					\
  src/test/PrimeField.cpp src/test/MonoMonoid.cpp src/test/Scanner.cpp	\
  src/test/MathicIO.cpp

else

check:
	@echo
	@echo "Configure did not locate gtest, so unittests cannot be run."

endif
