cmake_minimum_required (VERSION 2.8)
project (BunnyKiller)

# Deny in-source compilation
if ("${CMAKE_BINARY_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  message (FATAL_ERROR "Source and build directories are the same.")
endif ()

# If we do not have C++11, don't bother
include (CheckCXXCompilerFlag)

check_cxx_compiler_flag ("-std=c++11" COMPILER_SUPPORTS_CXX11)

if (NOT COMPILER_SUPPORTS_CXX11)
    message (FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support")
endif ()

set (CMAKE_CXX_STANDARD          11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS        OFF)

# tinyXML setup
add_library (tinyXML STATIC src/External/TinyXML/tinyxml2.cpp)

# Embree setup
set (EMBREE_STATIC_LIB                  ON  CACHE BOOL "Use Embree static library." FORCE)
set (EMBREE_ISPC_SUPPORT                OFF CACHE BOOL "Embree ISPC support." FORCE)
set (EMBREE_TUTORIALS                   OFF CACHE BOOL "Embree tutorials" FORCE)
set (EMBREE_STAT_COUNTERS               OFF CACHE BOOL "Statistic counters." FORCE)
set (EMBREE_RAY_MASK                    ON  CACHE BOOL "Ray mask support." FORCE)
set (EMBREE_BACKFACE_CULLING            OFF CACHE BOOL "Backface culling." FORCE)
set (EMBREE_INTERSECTION_FILTER         ON  CACHE BOOL "Intersection mask callback." FORCE)
set (EMBREE_INTERSECTION_FILTER_RESTORE OFF CACHE BOOL "Restores ray hit." FORCE)
set (EMBREE_TASKING_SYSTEM              "INTERNAL" CACHE STRING "Don't use Intel TBB" FORCE)
set (EMBREE_STATIC_RUNTIME              ON CACHE BOOL "Use static C++ library." FORCE)

add_subdirectory (src/External/Embree)

# Compiler flags
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif ()

if (WIN32)
    # MSVC flags
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /O2 /bigobj")
    
    string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    
    if (CMAKE_GENERATOR MATCHES "Visual Studio")
        set_property (GLOBAL PROPERTY USE_FOLDERS ON)
    endif ()
    
    # Stupid sanity check, pssibly CMake bugs' workaround
    if ("${MSVC_RUNTIME_LIBRARY_CONFIG}" STREQUAL "")
            set (MSVC_RUNTIME_LIBRARY_CONFIG "/MT")
    endif ()
    
    foreach (flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE)
        if ("${MSVC_RUNTIME_LIBRARY_CONFIG}" STREQUAL "/MT")
            string (REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
        elseif ("${MSVC_RUNTIME_LIBRARY_CONFIG}" STREQUAL "/MD")
            string (REPLACE "/MT" "/MD" ${flag_var} "${${flag_var}}")
        else ()
            string (REPLACE "/MD" "${MSVC_RUNTIME_LIBRARY_CONFIG}" ${flag_var} "${${flag_var}}")
            string (REPLACE "/MT" "${MSVC_RUNTIME_LIBRARY_CONFIG}" ${flag_var} "${${flag_var}}")
        endif ()
    endforeach ()
else (WIN32)
    # Rest of compilers flags
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -fvisibility=hidden -ffast-math -fno-finite-math-only")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Winvalid-pch -Wno-implicit-fallthrough")
    
        set (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
        set (CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem")
    
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
        set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=512")
    endif ()
endif (WIN32)

# Main executable
add_executable (BunnyKiller src/bunnykiller.cpp)

# Root include path for BunnyKiller's headers
target_include_directories (BunnyKiller PUBLIC src)

# Include path for External dependencies
target_include_directories (BunnyKiller SYSTEM PUBLIC src/External)

target_link_libraries (BunnyKiller tinyXML embree)

# IDE support
file (GLOB_RECURSE Headers "src/*.h")
add_custom_target (IDE SOURCES ${Headers})

