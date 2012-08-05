#!/bin/env bash

allCPPFLAGS="-Wall -Wextra -Wno-uninitialized -Wno-unused-parameter"

# optimized build
relName="rel";
relCPPFLAGS="-O2 -DNDEBUG";
relMakeArgs="";

# optimized build with asserts
relAssertName="relass";
relAssertCPPFLAGS="-O2 -DMEMTAILOR_DEBUG -DMATHIC_DEBUG -DMATHICGB_DEBUG";
relAssertMakeArgs="";

# debug build with asserts
debName="deb";
debCPPFLAGS="-g -DMEMTAILOR_DEBUG -DMATHIC_DEBUG -DMATHICGB_DEBUG";
debMakeArgs="";

# debug build without asserts
debNoAssertName="debnoass";
debNoAssertCPPFLAGS="-g -DNDEBUG";
debNoAssertMakeArgs="";

# profile build
proName="pro";
proCPPFLAGS="-g -pg -DNDEBUG -O2";
proMakeArgs="";

# analyze build: a larger amount of warnings turned on
anaName="ana";
anaCPPFLAGS="-fsyntax-only -O1 -Wfloat-equal -Wundef\
  -Wno-endif-labels -Wshadow -Wlarger-than-1000 -Wpointer-arith \
  -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wsign-compare \
  -Waggregate-return -Wmissing-noreturn -Wmissing-format-attribute \
  -Wno-multichar -Wno-deprecated-declarations -Wpacked \
  -Wno-redundant-decls -Wunreachable-code -Winline \
  -Wno-invalid-offsetof -Winvalid-pch -Wlong-long \
  -Wdisabled-optimization -D DEBUG -Werror"
anaMakeArgs="";

function makeTarget {
  projectName="$1";
  targetName="$2";
  makeArgs="$3";
  configureArgs="$4";
  if [ "$5" = "" ]; then
    dependency="";
  else
    dependency="$5$targetName";
  fi

  targetDir="$projectName/$targetName";
  prefix="\${PWD}/installed/$targetName";

  echo "$projectName$targetName: ${projectName}basic $dependency"
  echo $'\t'"rm -rf \"$targetDir\""
  echo $'\t'"mkdir -p \"$targetDir\" \"$prefix/lib/pkgconfig\";"
  echo $'\t'"( \\"
  echo $'\t'"  cd \"$targetDir\"; \\"
  echo $'\t'"  export PKG_CONFIG_PATH=\"$prefix/lib/pkgconfig\"; \\";
  echo $'\t'"  ../configure --prefix=\"$prefix\" CXXFLAGS=\"\" CPPFLAGS=\"$configureArgs\"; \\"
  echo $'\t'"  make $makeArgs install; \\"
  echo $'\t'");"
  echo "$targetName: $projectName$targetName";
}

function makeProject {
  name="$1";
  gitUrl="$2";
  dep="$3";
  echo "${name}basic:"
  echo $'\t'"if [ ! -e \"$name/\" ]; then git clone $gitUrl; fi;"
  echo $'\t'"if [ ! -e \"$name/configure\" ]; then (cd $name/; ./autogen.sh;); fi;"
  # It is intentional that ana is not built by default as it is
  # intended to check for warnings and those warnings are treated as
  # errors.  It should be possible to get a successful build even if
  # some of the many very strict warnings that are enabled for that
  # build target are not silenced currently.
  echo "$name: $name$relName"\
       "$name$relAssertName $name$debName"\
       "$name$debNoAssertName $name$proName" # $name$anaName
  makeTarget "$name" "$relName" "$relMakeArgs" "$relCPPFLAGS" "$dep";
  makeTarget "$name" "$relAssertName" "$relAssertMakeArgs" "$relAssertCPPFLAGS" "$dep";
  makeTarget "$name" "$debName" "$debMakeArgs" "$debCPPFLAGS" "$dep";
  makeTarget "$name" "$debNoAssertName" "$debNoAssertMakeArgs" "$debNoAssertCPPFLAGS" "$dep";
  makeTarget "$name" "$proName" "$proMakeArgs" "$proCPPFLAGS" "$dep";
  makeTarget "$name" "$anaName" "$anaMakeArgs" "$anaCPPFLAGS" "$dep";
}

# Causes 8 parallel tasks including within called makefiles. User setting
# overwrites this if the user has specified -jX on the command line.
echo "MAKEFLAGS += -j8 V=0"
echo "all: memtailor mathic mathicgb"
makeProject "memtailor" "https://github.com/broune/memtailor.git" "";
makeProject "mathic" "https://github.com/broune/mathic.git" "memtailor";
makeProject "mathicgb" "https://github.com/broune/mathicgb.git" "mathic";
