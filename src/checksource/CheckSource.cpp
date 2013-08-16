// Program for checking format of code, with respect to proper
// copyright header, no tabs, order of includes and so on.
//
// This code is not pretty, especially not the duplication of Scanner,
// but this makes this whole thing easy to compile and this is VERY
// quick-and-dirty code anyway. Anyone who wants to clean up thing up
// with proper integration into the build system is more than welcome
// to do so.
#include "Scanner.cpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <string>

void error(std::string str) {
  throw std::runtime_error(str);
}

bool endsWith(const std::string& a, const std::string& b) {
  return a.size() >= b.size() && a.substr(a.size() - b.size()) == b;
}

bool matchInclude(Scanner& in, bool& system) {
  if (!in.match("#include "))
    return false;

  char delim = '\"';
  system = false;
  if (!in.match('\"')) {
    if (in.match('<')) {
      system = true;
      delim = '>';
    } else
      in.reportError("Expected < or \" after #include.");
  }

  in.readUntil(delim);
  in.expect(delim);
  return true;
}

std::string stripEnding(const std::string& str) {
  const auto index = str.find_last_of('.');
  if (index == std::string::npos)
    return str;
  else
    return str.substr(0, index);
}



void checkCopyrightHeader(Scanner& in) {
  const char* const l1 = 
    "// MathicGB copyright 2012 all rights reserved. "
    "MathicGB comes with ABSOLUTELY\n";
  const char* const l2 =
    "// NO WARRANTY and is licensed as GPL v2.0 or later - see LICENSE.txt.\n";
  in.expect(l1);
  in.expect(l2);
}

void checkStdInc(Scanner& in) {
  in.eatWhite();
  if (!in.match("#include \"mathicgb/stdinc.h\"\n"))
      in.expect("#include \"stdinc.h\"\n");
}

void checkIncludes(Scanner& in) {
  bool sawSystem = false;
  bool system;
  while (matchInclude(in, system)) {
    if (sawSystem && !system)
      in.reportError("#include < > must come after all #include \" \".");
    sawSystem = system;
  }
}

void checkInclusionGuard(Scanner& in, const std::string& filename) {
  auto normalize = [](const std::string& str) {
    std::string norm;
    for (size_t i = 0; i < str.size(); ++i)
      if (str[i] != '_')
        norm += std::toupper(str[i]);
    return norm;
  };

  in.expect("#ifndef ");
  const auto macroName = in.readString();
  const auto fileNorm = normalize
    ("MATHICGB_" + stripEnding(filename) + "_GUARD");

  if (fileNorm != normalize(macroName)) {
    std::ostringstream err;
    err << "Inclusion guard name does not match file name.\n";
    err << "Normalization removes _ and converts to upper case.\n";
    err << "  Filename normalizes to: " << fileNorm << '\n';
    err << "macro name normalizes to: " << normalize(macroName) << '\n';
    in.reportError(err.str());
  }

  in.expect("#define ");
  std::string macro2;
  macro2 = in.readString();
  if (macroName != macro2) {
    std::ostringstream out;
    out << "Inclusion guard #ifndef and #define macro name mismatch.\n";
    out << "#ifndef macro name: " << macroName << '\n';
    out << "#define macro name: " << macro2 << '\n';
    in.reportError(out.str());
  }
}

void checkOther(const std::string& filename) {
  std::ifstream file(filename.c_str(), std::ios_base::binary);
  file.peek();
  if (!file)
    error("could not open file");
  Scanner in(file);
  bool mgbNamespace = false;
  while (!in.matchEOF()) {
    if (in.peek() == '\r')
      in.reportError("File contains dos/windows line ending character \\r.");
    if (in.peek() == '\t')
      in.reportError("File contains a tab.");
    if (
      in.match("Local Variables:") ||
      in.match("compile-command:") ||
      in.match("indent-tabs-mode:")
    )
      in.reportError("File contains emacs-specific command comments.");
    if (in.match("namespace mgb") || in.match("MATHICGB_NAMESPACE_BEGIN"))
      mgbNamespace = true;
    else
      in.get();
  }
  if (!mgbNamespace)
    in.reportError("MATHICGB_NAMESPACE_BEGIN does not appear in file");
}

void checkFile(std::string filename) {
  try {
    std::cout << "Checking file " << filename << std::endl;
    const bool hpp = endsWith(filename, ".hpp");
    const bool cpp = endsWith(filename, ".cpp");
    if (!hpp && !cpp)
      return;
    checkOther(filename);

    std::ifstream file(filename.c_str());
    if (!file)
      error("could not open file");
    Scanner in(file);

    if (in.peekWhite())
      in.reportError
        ("File should start with copyright header, not whitespace.");
      
    checkCopyrightHeader(in);
    if (in.peekWhite())
      in.reportError
        ("There should not be whitespace after the copyright header.");

    if (cpp)
      checkStdInc(in);
    else
      checkInclusionGuard(in, filename);
    checkIncludes(in);
  } catch (std::exception& e) {
    std::cout << "*** ERROR in file " << filename << " ***\n"
      << e.what() << std::endl;
  }
}

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i)
    checkFile(argv[i]);
  return 0;
}
