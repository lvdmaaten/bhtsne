
# Check if examples are enabled
if(NOT OPTION_BUILD_EXAMPLES)
    return()
endif()

# Example applications
add_subdirectory(bhtsne_cmd)
add_subdirectory(result2svg)
