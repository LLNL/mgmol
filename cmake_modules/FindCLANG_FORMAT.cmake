# Find clang-format
#
# CLANG_FORMAT_EXECUTABLE   - Path to clang-format executable
# CLANG_FORMAT_FOUND        - True if clang-format executable was found
# CLANG_FORMAT_VERSION      - The version of clang-format found

find_program(CLANG_FORMAT_EXECUTABLE
  NAMES clang-format clang-format-6.0
  DOC "clang-format executable")
mark_as_advanced(CLANG_FORMAT_EXECUTABLE)

# Extract version
if (CLANG_FORMAT_EXECUTABLE)
  execute_process(COMMAND ${CLANG_FORMAT_EXECUTABLE} -version
    OUTPUT_VARIABLE clang_format_version
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (clang_format_version MATCHES "^clang-format version .*")
    # clang_format_version sample: "clang-format version 6.0.0-1ubuntu2 (tags/RELEASE_600/final)"
    string(REGEX REPLACE "clang-format version ([.0-9]+).*" "\\1" 
      CLANG_FORMAT_VERSION "${clang_format_version}")
  endif()
  if (NOT CLANG_FORMAT_VERSION STREQUAL "6.0.0" AND
      NOT CLANG_FORMAT_VERSION STREQUAL "6.0.1")
    message(SEND_ERROR "Wrong clang-format version. Please use clang-format 6.0.")
  endif()
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CLANG_FORMAT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLANG_FORMAT REQUIRED_VARS CLANG_FORMAT_EXECUTABLE VERSION_VAR CLANG_FORMAT_VERSION)
