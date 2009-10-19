# - Try to find Trilinos

# Copyright (c) 2009, Nico Schlömer, <nico.schloemer@ua.ac.be>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

INCLUDE(CheckLibraryExists)

SET(_trilinos_INCLUDE_SEARCH_DIRS
  C:/trilinos/include
  C:/trilinos
  "$ENV{ProgramFiles}/trilinos"
  /usr/local/trilinos
  /opt/trilinos
  /usr/local/trilinos/current/include
  /opt/trilinos-10.0-r1/LINUX_MPI_64BIT/include/
)

FIND_PATH( Trilinos_INCLUDE_DIR
           Trilinos_version.h
           ${_trilinos_INCLUDE_SEARCH_DIRS} )

IF (Trilinos_INCLUDE_DIR)
   SET(Trilinos_FOUND TRUE)
ELSE (Trilinos_INCLUDE_DIR)
   SET(Trilinos_FOUND FALSE)
ENDIF (Trilinos_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments
IF (Trilinos_FOUND)
  IF (NOT Trilinos_FIND_QUIETLY)
    MESSAGE(STATUS "Found Trilinos: ${Trilinos_INCLUDE_DIR}")
  ENDIF (NOT Trilinos_FIND_QUIETLY)

  SET( Trilinos_LIB_DIR ${Trilinos_INCLUDE_DIR}/../lib )

  # ----------------------------------------------------------------------------
  FOREACH(COMPONENT ${Trilinos_FIND_COMPONENTS})
    STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
    STRING(TOLOWER ${COMPONENT} LOWERCOMPONENT)

    # find library
    FIND_LIBRARY( Trilinos_${UPPERCOMPONENT}_LIBRARY
                  NAMES ${LOWERCOMPONENT}
                  HINTS ${Trilinos_LIB_DIR}
                )
    IF (NOT Trilinos_${UPPERCOMPONENT}_LIBRARY)
      MESSAGE(FATAL_ERROR "Could NOT find Trilinos component ${COMPONENT}")
    ELSE (NOT Trilinos_${UPPERCOMPONENT}_LIBRARY)
      SET( Trilinos_${UPPERCOMPONENT}_FOUND TRUE )
      SET( Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${Trilinos_${UPPERCOMPONENT}_LIBRARY})
      IF (NOT Trilinos_FIND_QUIETLY)
        MESSAGE(STATUS "Found Trilinos component: ${COMPONENT}")
      ENDIF (NOT Trilinos_FIND_QUIETLY)
    ENDIF(NOT Trilinos_${UPPERCOMPONENT}_LIBRARY)

  ENDFOREACH(COMPONENT)
  # ----------------------------------------------------------------------------

ELSE (Trilinos_FOUND)
   IF (Trilinos_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could NOT find Trilinos")
   ENDIF (Trilinos_FIND_REQUIRED)
ENDIF (Trilinos_FOUND)


MARK_AS_ADVANCED(Trilinos_INCLUDE_DIR)
