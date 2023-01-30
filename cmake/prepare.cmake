## Enable Clang sanitizer for debug builds
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_DEBUG 
    "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" 
    CACHE STRING "")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG 
    "${CMAKE_EXE_LINKER_FLAGS_DEBUGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" 
    CACHE STRING "")
endif()

# macro that gets all subdirectories
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      IF(NOT ${child} STREQUAL "build")
        LIST(APPEND dirlist ${child})
      ENDIF()
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()