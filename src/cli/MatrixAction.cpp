#include "mathicgb/stdinc.h"
#include "MatrixAction.hpp"

#include "mathicgb/QuadMatrix.hpp"
#include "mathicgb/SparseMatrix.hpp"
#include "mathicgb/F4MatrixReducer.hpp"
#include "mathicgb/SparseMatrix.hpp"
#include "mathicgb/CFile.hpp"
#include <mathic.h>
#include <limits>
#include <fstream>
#include <iostream>

namespace {
  static const char* QuadMatrixExtension = ".qmat";
  static const char* LowerRightMatrixExtension = ".brmat";
  static const char* ReducedLowerRightMatrixExtension = ".rbrmat";

  /// Returns true if the file exists - or more precisely if it can be opened
  /// for reading. Using this function to create a file only if it does not
  /// exist implies a race condition in that the file could have been crated
  /// after checking if it exists and before (re)creating it. Do not use this
  /// approach if that is not acceptable. The advantage here is that this is
  /// portable. There could be a solution with freopen, but unfortunately
  /// freopen is allowed to fail on any change to the mode so it is not
  /// a portable solution.
  bool fileExists(const std::string fileName) {
    return CFile(fileName, "r", CFile::NoThrowTag()).hasFile();
  }
}

MatrixAction::MatrixAction():
  mParams(1, std::numeric_limits<size_t>::max()) {
  mParams.registerFileNameExtension(QuadMatrixExtension);
  mParams.registerFileNameExtension(LowerRightMatrixExtension);
  mParams.registerFileNameExtension(ReducedLowerRightMatrixExtension);
  mParams.registerFileNameExtension(".");
}

void MatrixAction::directOptions(
  std::vector<std::string> tokens,
  mic::CliParser& parser
) {
  mParams.directOptions(tokens, parser);
}

void MatrixAction::performAction() {
  mParams.perform();
  for (size_t i = 0; i < mParams.inputFileCount(); ++i) {
    const auto fileNameStem = mParams.inputFileNameStem(i);
    const auto extension = mParams.inputFileNameExtension(i);
    const auto quadFileName = fileNameStem + QuadMatrixExtension;
    const auto lowerRightFileName = fileNameStem + LowerRightMatrixExtension;
    const auto reducedLowerRightFileName =
      fileNameStem + ReducedLowerRightMatrixExtension;
    std::string inputFileName;

    SparseMatrix lowerRightMatrix;
    SparseMatrix::Scalar modulus;
    if (
      extension == QuadMatrixExtension ||
      extension == "." ||
      extension == ""
    ) {
      inputFileName = quadFileName;
      CFile file(quadFileName, "rb");
      QuadMatrix matrix;
      modulus = matrix.read(file.handle());
      fclose(file.handle());
      lowerRightMatrix = F4MatrixReducer(modulus).reduceToBottomRight(matrix);

      if (!fileExists(lowerRightFileName)) {
        CFile file(lowerRightFileName, "wb");
        lowerRightMatrix.write(modulus, file.handle());
      }
    } else if (extension == LowerRightMatrixExtension) {
      inputFileName = lowerRightFileName;
      CFile file(lowerRightFileName, "rb");
      modulus = lowerRightMatrix.read(file.handle());
    } else {
      mathic::reportError
        ("Unknown input file extension of " + mParams.inputFileName(i));
    }

    lowerRightMatrix = F4MatrixReducer(modulus).
      reducedRowEchelonForm(lowerRightMatrix);
    lowerRightMatrix.sortRowsByIncreasingPivots();

    if (!fileExists(reducedLowerRightFileName)) {
      CFile file(reducedLowerRightFileName.c_str(), "wb");
      lowerRightMatrix.write(modulus, file.handle());
    } else {
      SparseMatrix referenceMatrix;
      CFile file(reducedLowerRightFileName.c_str(), "rb");
      referenceMatrix.read(file.handle());
      if (lowerRightMatrix != referenceMatrix) {
        std::cerr << "Reducing " << inputFileName
          << " does not yield the matrix "
          << reducedLowerRightFileName << ".\n";
      } else if (tracingLevel > 0) {
        std::cerr << "Match for " << inputFileName 
          << " -> " << ReducedLowerRightMatrixExtension << ".\n";
      }
    }
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
