#include "stdinc.h"
#include "SparseMatrix.hpp"

std::ostream& operator<<(std::ostream& out, SparseMatrix& matrix) {
  matrix.print(out);
  return out;
}
