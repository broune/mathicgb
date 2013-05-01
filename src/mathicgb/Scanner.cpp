#include "stdinc.h"
#include "Scanner.hpp"

#include <mathic.h>
#include <limits>
#include <sstream>
#include <cstring>

void reportSyntaxError(std::string s) {
  mathic::reportError(s);
}

void Scanner::reportError(std::string msg) const {
  reportSyntaxError(msg);
}

static const size_t BufferSize = 10 * 1024;

Scanner::Scanner(FILE* input):
  mFile(input),
  mStream(0),
  mLineCount(1),
  mChar(' '),
  mBuffer(BufferSize),
  mBufferPos(mBuffer.end())
{
  get();
}

Scanner::Scanner(std::istream& input):
  mFile(0),
  mStream(&input),
  mLineCount(1),
  mChar(' '),
  mBuffer(BufferSize),
  mBufferPos(mBuffer.end())
{
  get();
}

Scanner::Scanner(const char* const input):
  mFile(0),
  mStream(0),
  mLineCount(1),
  mChar(' '),
  mBuffer(input, input + std::strlen(input)),
  mBufferPos(mBuffer.end())
{
  get();
}

void Scanner::expect(const char* str) {
  MATHICGB_ASSERT(str != 0);

  eatWhite();

  const char* it = str;
  while (*it != '\0') {
    int character = get();
    if (*it == character) {
      ++it;
      continue;
    }

    // Read the rest of what is there to improve error message.
    // TODO: read at least one char in total even if not alnum.
    std::ostringstream got;
    if (character == EOF && it == str)
      got << "no more input";
    else {
      got << '\"' << std::string(str, it);
      if (isalnum(character))
        got << static_cast<char>(character);
      while (isalnum(peek()))
        got << static_cast<char>(get());
      got << '\"';
    }

    reportErrorUnexpectedToken(str, got.str());
  }
}

void Scanner::expectEOF() {
  eatWhite();
  if (get() != EOF)
    reportErrorUnexpectedToken("no more input", "");
}

void Scanner::errorExpectTwo(char a, char b, int got) {
  MATHICGB_ASSERT(a != got && b != got);
  std::ostringstream err;
  err << '\'' << a << "' or '" << b << '\'';
  reportErrorUnexpectedToken(err.str(), got);
}

void Scanner::errorExpectOne(char expected, int got) {
  MATHICGB_ASSERT(expected != got);
  std::ostringstream err;
  err << '\'' << expected << '\'';
  reportErrorUnexpectedToken(err.str(), got);
}

void Scanner::reportErrorUnexpectedToken(const std::string& expected, int got) {
  std::ostringstream gotDescription;
  if (got == EOF)
    gotDescription << "no more input";
  else
    gotDescription << '\'' << static_cast<char>(got) << '\'';
  reportErrorUnexpectedToken(expected, gotDescription.str());
}

void Scanner::reportErrorUnexpectedToken(
  const std::string& expected,
  const std::string& got
) {
  std::ostringstream errorMsg;
  errorMsg << "Expected " << expected;
  if (got != "")
    errorMsg << ", but got " << got;
  errorMsg << '.';
  reportSyntaxError(errorMsg.str());
}

int Scanner::readBuffer() {
  size_t read;
  if (mFile != 0) {
    if (mBuffer.size() < mBuffer.capacity() && (feof(mFile) || ferror(mFile)))
      return EOF;
    mBuffer.resize(mBuffer.capacity());
    read = fread(&mBuffer[0], 1, mBuffer.capacity(), mFile);
  } else if (mStream != 0) {
    MATHICGB_ASSERT(mStream != 0);
    if (mBuffer.size() < mBuffer.capacity() && !mStream->good())
      return EOF;
    mBuffer.resize(mBuffer.capacity());
    mStream->read(reinterpret_cast<char*>(mBuffer.data()), mBuffer.size());
    read = mStream->gcount();
  } else
    return EOF;
  mBuffer.resize(read);
  mBufferPos = mBuffer.begin();
  if (read == 0)
    return EOF;
  const char c = *mBufferPos;
  ++mBufferPos;
  return c;
}
