CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS}
                      ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS}
                  ${Trilinos_TPL_LIBRARY_DIRS})
# ------------------------------------------------------------------------------
SET (medley_SRCS
      dummy.f
      strong-scaling-medley.cpp
   )
SET(MY_EXECUTABLE
     "medley.exe")
ADD_EXECUTABLE(${MY_EXECUTABLE}
                ${medley_SRCS})
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                   )
# ------------------------------------------------------------------------------
SET (dot_SRCS
      dummy.f
      strong-scaling-dotproduct.cpp
   )
SET(MY_EXECUTABLE
     "dot.exe")
ADD_EXECUTABLE(${MY_EXECUTABLE}
                ${dot_SRCS})
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
                       ${Trilinos_LIBRARIES}
                       ${Trilinos_TPL_LIBRARIES}
                   )
## ------------------------------------------------------------------------------
#SET (matvec_SRCS
#      dummy.f
#      strong-scaling-matvec.cpp
#   )
#SET(MY_EXECUTABLE
#     "matvec.exe")
#ADD_EXECUTABLE(${MY_EXECUTABLE}
#                ${matvec_SRCS})
#TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
#                       ${Trilinos_LIBRARIES}
#                       ${Trilinos_TPL_LIBRARIES}
#                   )
# ------------------------------------------------------------------------------
