#include "stdinc.h"
#include "SparseMatrix.hpp"

std::ostream& operator<<(std::ostream& out, const SparseMatrix& matrix) {
  matrix.print(out);
  return out;
}
