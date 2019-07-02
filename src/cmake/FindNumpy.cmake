#  Variables set by this script
#  NUMPY_FOUND
#  NUMPY_INCLUDE_DIRS

# Initial version of this file was from MITK
# and specifically http://docs.mitk.org/nightly/FindNumpy_8cmake_source.html
# MITK license is BSD-style (at Jan 2016). Full version is at
# http://mitk.org/git/?p=MITK.git;a=blob_plain;f=LICENSE.txt
# Copyright (c) 2003-2012 German Cancer Research Center,
# Division of Medical and Biological Informatics

# Minor modifications by Kris Thielemans
# - removed special handling of spaces in Python executable
# - renamed NUMPY_INCLUDE_DIR to NUMPY_INCLUDE_DIRS
# TODO: should find npymath core library

find_package(PackageHandleStandardArgs)
     
# search for the runtime if it's not already set
if( NOT DEFINED PYTHON_EXECUTABLE )
  FIND_PACKAGE(PythonLibs REQUIRED)
  FIND_PACKAGE(PythonInterp REQUIRED)
endif()

set(_python ${PYTHON_EXECUTABLE})
# KT commented out next 3 lines. not sure if you need them.
#if(UNIX)
#  STRING(REPLACE " " "\ " _python ${PYTHON_EXECUTABLE})
#endif()

execute_process (
  COMMAND ${_python} -c "from __future__ import print_function; import os; os.environ['DISTUTILS_USE_SDK']='1'; import numpy.distutils; print(numpy.distutils.misc_util.get_numpy_include_dirs()[0])"
  OUTPUT_VARIABLE output
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

if(DEFINED output AND EXISTS ${output} )
  set (NUMPY_INCLUDE_DIRS ${output})
endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Numpy DEFAULT_MSG NUMPY_INCLUDE_DIRS)

MARK_AS_ADVANCED (
  NUMPY_INCLUDE_DIRS
)
