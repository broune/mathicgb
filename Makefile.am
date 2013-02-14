# options passed to aclocal, which is a tool for making macroes visible to
# autoconf. We use -I to tell aclocal where we put the local macros.
ACLOCAL_AMFLAGS = -I build/autotools/m4

# Options passed to the C and C++ PreProcessor (CPP) and compiler
AM_CPPFLAGS = -I${top_srcdir}/
libmathicgb_la_CPPFLAGS = $(DEPS_CFLAGS)

# tell Libtool what the name of the library is.
noinst_LTLIBRARIES = libmathicgb.la

# set the C++ compiler to include src/
AM_CXXFLAGS=-I$(top_srcdir)/src/ -std=gnu++0x
AM_LDFLAGS=

# libraries that are needed by this library
libmathicgb_la_LIBADD= $(DEPS_LIBS)

# the sources that are built to make libmathicgb. Listing the headers in
# sources ensure that those files are included in distributions.
libmathicgb_la_SOURCES = src/mathicgb/BjarkeGeobucket2.cpp		\
  src/mathicgb/BjarkeGeobucket2.hpp src/mathicgb/BjarkeGeobucket.cpp	\
  src/mathicgb/BjarkeGeobucket.hpp src/mathicgb/BuchbergerAlg.cpp	\
  src/mathicgb/BuchbergerAlg.hpp src/mathicgb/ChainedHashTable.cpp	\
  src/mathicgb/ChainedHashTable.hpp src/mathicgb/DivisorLookup.hpp	\
  src/mathicgb/DivisorLookup.cpp src/mathicgb/DivLookup.hpp		\
  src/mathicgb/FreeModuleOrder.cpp src/mathicgb/FreeModuleOrder.hpp	\
  src/mathicgb/GroebnerBasis.cpp src/mathicgb/GroebnerBasis.hpp	\
  src/mathicgb/HashTourReducer.cpp src/mathicgb/HashTourReducer.hpp	\
  src/mathicgb/Ideal.cpp src/mathicgb/Ideal.hpp			\
  src/mathicgb/io-util.cpp src/mathicgb/io-util.hpp			\
  src/mathicgb/KoszulQueue.cpp src/mathicgb/KoszulQueue.hpp		\
  src/mathicgb/MonomialHashTable.hpp					\
  src/mathicgb/MonTableDivList.hpp src/mathicgb/MonTableKDTree.hpp	\
  src/mathicgb/MonTableNaive.cpp src/mathicgb/MonTableNaive.hpp		\
  src/mathicgb/MTArray.cpp src/mathicgb/MTArray.hpp			\
  src/mathicgb/PairTriangle.cpp src/mathicgb/PairTriangle.hpp		\
  src/mathicgb/Poly.cpp src/mathicgb/Poly.hpp				\
  src/mathicgb/PolyBasis.cpp src/mathicgb/PolyBasis.hpp			\
  src/mathicgb/PolyGeoBucket.cpp src/mathicgb/PolyGeoBucket.hpp		\
  src/mathicgb/PolyHashReducer.cpp src/mathicgb/PolyHashReducer.hpp	\
  src/mathicgb/PolyHashTable.cpp src/mathicgb/PolyHashTable.hpp		\
  src/mathicgb/PolyHeap.cpp src/mathicgb/PolyHeap.hpp			\
  src/mathicgb/PolyReducer.cpp src/mathicgb/PolyReducer.hpp		\
  src/mathicgb/PolyRing.cpp src/mathicgb/PolyRing.hpp			\
  src/mathicgb/Reducer.cpp src/mathicgb/Reducer.hpp			\
  src/mathicgb/ReducerDedup.hpp src/mathicgb/ReducerHash.hpp		\
  src/mathicgb/ReducerHashPack.hpp src/mathicgb/ReducerHelper.hpp	\
  src/mathicgb/ReducerNoDedup.hpp src/mathicgb/ReducerPack.hpp		\
  src/mathicgb/ReducerPackDedup.hpp src/mathicgb/SignatureGB.cpp	\
  src/mathicgb/SignatureGB.hpp src/mathicgb/SigSPairs.cpp		\
  src/mathicgb/SigSPairs.hpp src/mathicgb/SPairs.cpp			\
  src/mathicgb/SPairs.hpp src/mathicgb/stdinc.h				\
  src/mathicgb/TournamentReducer.cpp					\
  src/mathicgb/TournamentReducer.hpp src/mathicgb/SigSPairQueue.hpp	\
  src/mathicgb/SigSPairQueue.cpp src/mathicgb/SparseMatrix.hpp		\
  src/mathicgb/SparseMatrix.cpp src/mathicgb/QuadMatrixBuilder.hpp	\
  src/mathicgb/QuadMatrixBuilder.cpp src/mathicgb/TypicalReducer.cpp\
  src/mathicgb/TypicalReducer.cpp src/mathicgb/F4Reducer.hpp		\
  src/mathicgb/F4Reducer.cpp src/mathicgb/F4MatrixBuilder.hpp		\
  src/mathicgb/F4MatrixBuilder.cpp src/mathicgb/QuadMatrix.hpp		\
  src/mathicgb/QuadMatrix.cpp src/mathicgb/F4MatrixReducer.cpp		\
  src/mathicgb/F4MatrixReducer.hpp src/mathicgb/MonomialMap.hpp	\
  src/mathicgb/RawVector.hpp src/mathicgb/Atomic.hpp			\
  src/mathicgb/FixedSizeMonomialMap.hpp src/mathicgb/CFile.hpp		\
  src/mathicgb/CFile.cpp src/mathicgb/LogDomain.hpp			\
  src/mathicgb/LogDomain.cpp src/mathicgb/LogDomainSet.hpp \
  src/mathicgb/F4MatrixBuilder2.hpp src/mathicgb/F4MatrixBuilder2.cpp \
  src/mathicgb/LogDomainSet.cpp src/mathicgb/ProtoMatrix.hpp \
  src/mathicgb/ProtoMatrix.cpp src/mathicgb/F4MatrixProject.hpp \
  src/mathicgb/F4MatrixProjection.cpp

# When making a distribution file, Automake knows to include all files
# that are necessary to build the project. EXTRA_DIST specifies files
# to include beyond those used in the build process.
EXTRA_DIST = autogen.sh libs/gtest

bin_PROGRAMS = mgb

# set up the tests. Listing the headers in sources ensure that those
# files are included in distributions.
mgb_CPPFLAGS = $(DEPS_CFLAGS)
mgb_SOURCES = src/cli/GBMain.cpp src/cli/CommonParams.hpp				\
  src/cli/CommonParams.cpp src/cli/GBAction.hpp src/cli/GBAction.cpp	\
  src/cli/GBCommonParams.hpp src/cli/GBCommonParams.cpp					\
  src/cli/MatrixAction.cpp src/cli/MatrixAction.hpp						\
  src/cli/SigGBAction.hpp src/cli/SigGBAction.cpp						\
  src/cli/HelpAction.hpp src/cli/HelpAction.cpp
mgb_LDADD = $(top_builddir)/libmathicgb.la

# set up tests to run on "make check"
TESTS=unittest
check_PROGRAMS=$(TESTS)

unittest_CPPFLAGS = $(DEPS_CFLAGS)
unittest_CXXFLAGS=\
  -I$(top_srcdir)/libs/gtest/include\
  -I$(top_srcdir)/libs/gtest/\
  -I$(top_srcdir)/src/ -std=gnu++0x
unittest_LDADD = $(DEPS_LIBS) libmathicgb.la

test_LIBS=
unittest_SOURCES=src/test/FreeModuleOrderTest.cpp						\
  src/test/gtestInclude.cpp src/test/testMain.cpp src/test/gb-test.cpp	\
  src/test/ideals.cpp src/test/poly-test.cpp src/test/ideals.hpp		\
  src/test/SparseMatrix.cpp src/test/QuadMatrixBuilder.cpp				\
  src/test/F4MatrixBuilder.cpp src/test/F4MatrixReducer.cpp
