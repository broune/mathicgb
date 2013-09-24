#!/usr/bin/env bash

mildWarn="" #-Wall -Wextra"

#############################################
################## projects #################
#############################################

# memtailor
memtailorIndex=${#projectsName[@]}
projectsName+=("memtailor");
projectsGitUrl+=("https://github.com/broune/memtailor.git");
projectsDependencies+=("");

# mathic
mathicIndex=${#projectsName[@]}
projectsName+=("mathic");
projectsGitUrl+=("https://github.com/broune/mathic.git");
projectsDependencies+=("memtailor");

# mathicgb
mathicgbIndex=${#projectsName[@]}
projectsName+=("mathicgb");
projectsGitUrl+=("https://github.com/broune/mathicgb.git");
projectsDependencies+=("memtailor mathic");

#############################################
################## targets ##################
#############################################

# release
relIndex=${#targetsName[@]};
targetsName+=("rel");
targetsDescription+=("Release build. Optimized, no debug symbols, no asserts.");
targetsCPPFLAGS+=("-O2");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("");
targetsMakeArgs+=("");
targetsDefault+=("yes");

# optimized with asserts
relassIndex=${#targetsName[@]};
targetsName+=("relass");
targetsDescription+=("Optimized build with asserts. No debug symbols.");
targetsCPPFLAGS+=("-O2 -DMEMTAILOR_DEBUG -DMATHIC_DEBUG -DMATHICGB_DEBUG $mildWarn");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("");
targetsMakeArgs+=("");
targetsDefault+=("yes");

# debug with asserts
debIndex=${#targetsName[@]};
targetsName+=("deb");
targetsDescription+=("Debug build with asserts. Not optimized.");
targetsCPPFLAGS+=("-g -DMEMTAILOR_DEBUG -DMATHIC_DEBUG -DMATHICGB_DEBUG $mildWarn");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("");
targetsMakeArgs+=("");
targetsDefault+=("yes");

#debug without asserts
debnoassIndex=${#targetsName[@]};
targetsName+=("debnoass");
targetsDescription+=("Debug build without asserts. Not optimized.");
targetsCPPFLAGS+=("-g $mildWarn");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("");
targetsMakeArgs+=("");
targetsDefault+=("no");

# profile
proIndex=${#targetsName[@]};
targetsName+=("pro");
targetsDescription+=("Profile build with optimization and no asserts.");
targetsCPPFLAGS+=("-g -pg -O2 $mildWarn");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("-pg");
targetsMakeArgs+=("");
targetsDefault+=("no");

# profile with asserts
proassIndex=${#targetsName[@]};
targetsName+=("proass");
targetsDescription+=("Profile build with optimization and asserts.");
targetsCPPFLAGS+=("-g -pg -O2 -DMEMTAILOR_DEBUG -DMATHIC_DEBUG -DMATHICGB_DEBUG $mildWarn");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("-pg");
targetsMakeArgs+=("");
targetsDefault+=("no");

# analyze build: a larger amount of warnings turned on
anaIndex=${#targetsName[@]};
targetsName+=("ana");
targetsDescription+=("Build with many warnings turned on and -Werror.");
targetsCPPFLAGS+=("-Wall -Wextra -Wno-uninitialized -Wno-unused-parameter\
   -O1 -Wfloat-equal -Wundef\
  -Wno-endif-labels -Wshadow -Wlarger-than-1000 -Wpointer-arith \
  -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wsign-compare \
  -Waggregate-return -Wmissing-noreturn -Wmissing-format-attribute \
  -Wno-multichar -Wno-deprecated-declarations -Wpacked \
  -Wno-redundant-decls -Wunreachable-code -Winline \
  -Wno-invalid-offsetof -Winvalid-pch -Wlong-long \
  -Wdisabled-optimization -DMEMTAILOR_DEBUG -DMATHIC_DEBUG -DMATHICGB_DEBUG \
  -Werror");
targetsCXXFLAGS+=("");
targetsLDFLAGS+=("");
targetsMakeArgs+=("");
targetsDefault+=("no");


function makeHelpComment {
  echo "# This file was auto-generated using the script "
  echo "#   mathicgb/build/setup/make-Makefile.sh";
  echo "# As such it is recommended to change that script instead of "
  echo "# changing this file if you want to keep the changes in future."
  echo "#";
  echo "# The projects that you can build are"
  for ((i=0; i<${#projectsName[@]}; i++)); do
    echo "#";
    echo "#   Project: ${projectsName[i]}.";
    echo "#     gitUrl=${projectsGitUrl[i]}";
    echo "#     dependencies=${projectsDependencies[i]}";
  done

  echo "#";
  echo "# Within each project the targets than you can build are";
  for ((i=0; i<${#targetsName[@]}; i++)); do
    echo "#";
    echo "#   Target: ${targetsName[i]}. ${targetsDescription[i]}";
    echo "#     CPPFLAGS=${targetsCPPFLAGS[i]}";
    echo "#     CXXFLAGS=${targetsCXXFLAGS[i]}";
    echo "#     LDFLAGS=${targetsLDFLAGS[i]}";
    echo "#     makeArgs=${targetsMakeArgs[i]}";
    echo "#     buildsByDefault=${targetsDefault[i]}";
  done

  echo "#";
  echo "# To build everything that builds by default, type \"make\".";
  echo "#   (\"make\" will take quite a while to run.)";
  echo "# To build all rel targets, type \"make rel\".";
  echo "# To build target rel of mathic, type \"make mathicrel\".";
  echo "# This will automatically build memtailorrel too due to the";
  echo "# dependency.";
  echo;
}

function makeTarget {
  projectIndex="$1";
  projectName="${projectsName[projectIndex]}";
  projectDependencies="${projectsDependencies[projectIndex]}";

  targetIndex="$2";
  targetName="${targetsName[targetIndex]}";
  targetCPPFLAGS="${targetsCPPFLAGS[targetIndex]}";
  targetCXXFLAGS="${targetsCXXFLAGS[targetIndex]}";
  targetLDFLAGS="${targetsLDFLAGS[targetIndex]}";
  targetMakeArgs="${targetsMakeArgs[targetIndex]}";
  targetDefault="${targetsDefault[targetIndex]}";

  dependencies="";
  for depProject in $projectDependencies; do
    dependencies+=" $depProject$targetName";
  done

  name="$projectName$targetName"
  targetDir="$projectName/$targetName";
  prefix="\${PWD}/installed/$targetName";
  projectConfigureArgs="--prefix=\"$prefix\" \${${projectName}_conf}";

  echo "$name: ${projectName}BasicSetup $dependencies"
  echo $'\t'"rm -rf \"$targetDir\""
  echo $'\t'"mkdir -p \"$targetDir\" \"$prefix/lib/pkgconfig\";"
  echo $'\t'"( \\"
  echo $'\t'"  cd \"$targetDir\"; \\"
  echo $'\t'"  export PKG_CONFIG_PATH=\"$prefix/lib/pkgconfig\"; \\";
  echo $'\t'"  export CXXFLAGS=\"$targetCXXFLAGS\"; \\";
  echo $'\t'"  export CPPFLAGS=\"$targetCPPFLAGS\"; \\";
  echo $'\t'"  export LDFLAGS=\"$targetLDFLAGS\"; \\";
  echo $'\t'"  export GTEST_PATH=\"../..\"; \\";
  echo $'\t'"  ../configure $projectConfigureArgs; \\"
  echo $'\t'"  make $targetMakeArgs install; \\"
  echo $'\t'");"
  echo "$targetName: $name";
  if [ "$targetDefault" = "yes" ]; then
    echo "$projectName: $name";
  fi
}

function makeGTest {
  version="1.6.0";
  zipFile="gtest-$version.zip";
  url="http://googletest.googlecode.com/files/$zipFile";
  extractDir="gtest-$version/";

  echo "gtest:";
  echo $'\t'"rm -rf $zipFile";
  echo $'\t'"wget $url;";
  echo $'\t'"unzip $zipFile;";
  echo $'\t'"rm $zipFile;";
  echo $'\t'"rm -rf gtest;";
  echo $'\t'"mv $extractDir/ gtest;";
}

function makeProject {
  projectIndex="$1";
  name="${projectsName[projectIndex]}";
  gitUrl="${projectsGitUrl[projectIndex]}";

  dep="$3";
  echo "all: $name";
  echo "${name}BasicSetup: gtest"
  echo $'\t'"if [ ! -e \"$name/\" ]; then git clone $gitUrl; fi;"
  echo $'\t'"if [ ! -e \"$name/configure\" ]; then (cd $name/; ./autogen.sh;); fi;"

  for ((k=0; k<${#targetsName[@]}; k++)); do
    makeTarget "$projectIndex" "$k";
  done
}

function makeMakefile {
  echo "all:";
  makeHelpComment;
  makeGTest;

  # -j8: Causes 8 parallel tasks including within called makefiles. User
  # setting overwrites this if the user has specified -jX on the command line.
  #
  # V=0: Causes a non-verbose build.
  echo "MAKEFLAGS += -j8 V=0"

  for ((j=0; j<${#projectsName[@]}; j++)); do
    makeProject "$j";
  done
}

makeMakefile
