# Additional targets to perform clang-format/clang-tidy

# Get all project files - FIXME: this should also use the list of generated targets
IF(NOT CHECK_CXX_SOURCE_FILES)
    MESSAGE(FATAL_ERROR "Variable CHECK_CXX_SOURCE_FILES not defined - set it to the list of files to auto-format")
    RETURN()
ENDIF()

# Adding clang-format check and formatter if found
FIND_PROGRAM(CLANG_FORMAT "clang-format")
IF(NOT CLANG_FORMAT)
  RETURN()
ENDIF()
EXEC_PROGRAM(${CLANG_FORMAT} ${CMAKE_CURRENT_SOURCE_DIR} ARGS --version OUTPUT_VARIABLE CLANG_VERSION)
STRING(REGEX REPLACE ".*([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" CLANG_MAJOR_VERSION ${CLANG_VERSION})

IF(${CLANG_MAJOR_VERSION} EQUAL ${CLANG_FORMAT_VERSION} OR DEFINED ${APPLE})
    MESSAGE(WARNING "Found ${CLANG_FORMAT} version ${CLANG_MAJOR_VERSION}, this might lead to incompatible formatting")
    MESSAGE(WARNING "Please use Clang ${CLANG_FORMAT_VERSION} for formatting")
    RETURN()
ELSE()
    MESSAGE(STATUS "Found ${CLANG_FORMAT} version ${CLANG_FORMAT_VERSION}, adding formatting targets")
ENDIF()


IF(CLANG_FORMAT)
    ADD_CUSTOM_TARGET(
        format
        COMMAND
        ${CLANG_FORMAT}
        -i
        -style=file
        ${CHECK_CXX_SOURCE_FILES}
        COMMENT "Auto formatting of all source files"
    )

    ADD_CUSTOM_TARGET(
        check-format
        COMMAND
        ${CLANG_FORMAT}
        -style=file
        -output-replacements-xml
        ${CHECK_CXX_SOURCE_FILES} 
        # print output
        | tee ${CMAKE_BINARY_DIR}/check_format_file.txt | grep -c "replacement " | 
                tr -d "[:cntrl:]" && echo " replacements necessary" && cat ${CMAKE_BINARY_DIR}/check_format_file.txt
        # WARNING: fix to stop with error if there are problems
        COMMAND ! grep -c "replacement " 
                  ${CMAKE_BINARY_DIR}/check_format_file.txt > /dev/null
        COMMENT "Checking format compliance"
    )
ENDIF()
