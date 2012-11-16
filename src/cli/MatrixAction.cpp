#include "mathicgb/stdinc.h"
#include "MatrixAction.hpp"

#include "mathicgb/QuadMatrix.hpp"
#include "mathicgb/SparseMatrix.hpp"
#include <fstream>
#include <iostream>

MatrixAction::MatrixAction() {
  mParams.registerFileNameExtension(".mat");
}

void MatrixAction::directOptions(
  std::vector<std::string> tokens,
  mic::CliParser& parser
) {
  mParams.directOptions(tokens, parser);
}

void MatrixAction::performAction() {
  mParams.perform();

  QuadMatrix matrix;

  // @todo: fix leak of file on exception

  // read matrix
  SparseMatrix::Scalar modulus;
  {
    FILE* file = fopen((mParams.inputFileNameStem() + ".mat").c_str(), "rb");
    modulus = matrix.read(file);
    fclose(file);
  }

  // write matrix
  {
    FILE* file = fopen((mParams.inputFileNameStem() + ".out").c_str(), "wb");
    matrix.write(modulus, file);
    fclose(file);
  }
}

const char* MatrixAction::staticName() {
  return "matrix";
}

const char* MatrixAction::name() const {
  return staticName();
}

const char* MatrixAction::description() const {
  return "Perform matrix computations. "
    "The name of the matrix file is an optional direct parameter.";
}

const char* MatrixAction::shortDescription() const {
  return "Perform matrix computations.";
}

void MatrixAction::pushBackParameters(
  std::vector<mic::CliParameter*>& parameters
) {
  mParams.pushBackParameters(parameters);
}
