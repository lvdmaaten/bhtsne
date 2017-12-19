
# Check if tools are enabled
if(NOT OPTION_BUILD_TOOLS)
    return()
endif()

# Example applications
add_subdirectory(bhtsne_cmd)
add_subdirectory(performance_test)
add_subdirectory(performance_visualizer)
add_subdirectory(result2svg)
