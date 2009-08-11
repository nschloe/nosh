# - Try to find the NOX
# Once done this will define
#
#  CUPS_FOUND - system has Cups
#  CUPS_INCLUDE_DIR - the Cups include directory
#  CUPS_LIBRARIES - Libraries needed to use Cups
#  Set CUPS_REQUIRE_IPP_DELETE_ATTRIBUTE to TRUE if you need a version which
#  features this function (i.e. at least 1.1.19)

# Copyright (c) 2009, Nico Schlömer, <nico.schloemer@ua.ac.be>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


INCLUDE(CheckLibraryExists)

FIND_PATH(NOX_INCLUDE_DIR NOX.H )

FIND_LIBRARY(NOX_LIBRARIES NAMES nox )

IF (NOX_INCLUDE_DIR AND NOX_LIBRARIES)
   SET(NOX_FOUND TRUE)

ELSE  (NOX_INCLUDE_DIR AND NOX_LIBRARIES)
   SET(NOX_FOUND FALSE)
ENDIF (NOX_INCLUDE_DIR AND NOX_LIBRARIES)

IF (NOX_FOUND)
   IF (NOT NOX_FIND_QUIETLY)
      MESSAGE(STATUS "Found NOX: ${NOX_LIBRARIES}")
   ENDIF (NOT NOX_FIND_QUIETLY)
ELSE (NOX_FOUND)
   SET(NOX_LIBRARIES )
   IF (NOX_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could NOT find NOX")
   ENDIF (NOX_FIND_REQUIRED)
ENDIF (NOX_FOUND)


MARK_AS_ADVANCED(NOX_INCLUDE_DIR NOX_LIBRARIES)