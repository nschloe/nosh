CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(Nosh REQUIRED)

# Set the compiler.
SET(CMAKE_CXX_COMPILER ${Nosh_CXX_COMPILER})

# Set the name of your application.
PROJECT(my-application CXX)

# Set INCLUDE directories for compilation.
INCLUDE_DIRECTORIES(${Nosh_INCLUDE_DIRS}
                    ${Nosh_Trilinos_INCLUDE_DIRS})

# Set LINK directories for linking.
LINK_DIRECTORIES(${Nosh_LIBRARY_DIRS}
                 ${Nosh_Trilinos_LIBRARY_DIRS}
                 ${Nosh_Trilinos_TPL_LIBRARY_DIRS})

# Set the name of your executable.
SET(MY_EXECUTABLE "cont")

# Set the source files.
ADD_EXECUTABLE(${MY_EXECUTABLE} "continuation.cpp" "MyScalarField.cpp")

# Specify against which libraries they will be linked.
TARGET_LINK_LIBRARIES(${MY_EXECUTABLE}
                      ${Nosh_LIBRARIES}
                      ${Nosh_Trilinos_LIBRARIES}
                      ${Nosh_Trilinos_TPL_LIBRARIES})
