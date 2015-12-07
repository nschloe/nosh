set(input_files
  rectanglesmall
  pacman
  screw
  shell
  sphere
  cubesmall
  brick-w-hole
  )

# Download the files
INCLUDE(ExternalData)
set(
  ExternalData_URL_TEMPLATES
  "https://downloads.sourceforge.net/project/noshdata/%(algo)/%(hash)"
  )
foreach(file ${input_files})
  ExternalData_Expand_Arguments(
    noshTestFetchData
    OUT_DATA DATA{${CMAKE_SOURCE_DIR}/test/data/${file}.e}
    )
endforeach()
ExternalData_Add_Target(noshTestFetchData)

# convert and split the downloaded files
SET(FILES2)
SET(FILES7)
foreach(file ${input_files})
  add_custom_command(
    DEPENDS ${file}.e
    COMMAND meshio-convert ${file}.e ${file}.h5m
    OUTPUT ${file}.h5m
    )

  LIST(APPEND FILES2 ${file}-2.h5m)
  add_custom_command(
    DEPENDS ${file}.h5m
    # https://bitbucket.org/fathomteam/moab/issues/27/parallel-read-of-simple-file-fails
    COMMAND mbpart 2 -m ML_RB ${file}.h5m ${file}-2.h5m
    #COMMAND mbpart 2 -m ML_KWAY ${file}.h5m ${file}-2.h5m
    OUTPUT ${file}-2.h5m
    )

  LIST(APPEND FILES7 ${file}-7.h5m)
  add_custom_command(
    DEPENDS ${file}.h5m
    # https://bitbucket.org/fathomteam/moab/issues/27/parallel-read-of-simple-file-fails
    COMMAND mbpart 7 -m ML_RB ${file}.h5m ${file}-7.h5m
    #COMMAND mbpart 7 -m ML_KWAY ${file}.h5m ${file}-7.h5m
    OUTPUT ${file}-7.h5m
    )
endforeach()

add_custom_target(split2 ALL
  SOURCES ${FILES2}
  )

add_custom_target(split7 ALL
  SOURCES ${FILES7}
  )
