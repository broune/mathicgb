// MathicGB copyright 2012 all rights reserved. MathicGB comes with ABSOLUTELY
// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.
#include "mathicgb/stdinc.h"

// Includes a file from gtest that pulls in all of the implementation
// of gtest. The gtest docs recommend building gtest individually for
// each program rather than using an installed gtest and this is as
// easy a way of doing it as any. Especially because it guarantees that
// the compiler flags are the same, which is the whole point of the
// recommendation to build gtest for each program.

namespace mgb {}
using namespace mgb;

// the .. goes back from the include/ directory of gtest so we can
// enter the src directory.
#include <../src/gtest-all.cc>
