# Byte-compile and install nosh_mesh_io.
FIND_PACKAGE(PythonInterp 2 REQUIRED)

SET(MODULE_NAMES "__init__;mesh_io")

FOREACH(MODULE_NAME ${MODULE_NAMES})
  GET_FILENAME_COMPONENT(ABS_MODULE_PATH ${MODULE_NAME}.py ABSOLUTE)
  ADD_CUSTOM_TARGET(${MODULE_NAME} ALL
        COMMAND ${PYTHON_EXECUTABLE} -c "import py_compile; py_compile.compile('${ABS_MODULE_PATH}')"
        COMMENT "Byte compiling ${MODULE_NAME}.py"
        VERBATIM)
  INSTALL(FILES
    ${ABS_MODULE_PATH}
    ${MODULE_NAME}.pyc
    DESTINATION ${INSTALL_LIB_DIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/nosh
    )
ENDFOREACH()
