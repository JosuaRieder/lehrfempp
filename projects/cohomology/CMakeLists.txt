cmake_minimum_required(VERSION 3.12)

include(FetchContent)
if(NOT COMMAND FetchContent_MakeAvailable)
	macro(FetchContent_MakeAvailable NAME)
		FetchContent_GetProperties(${NAME})
		string(TOLOWER "${NAME}" lcName)
		if(NOT ${lcName}_POPULATED)
			FetchContent_Populate(${NAME})
			add_subdirectory("${${lcName}_SOURCE_DIR}" "${${lcName}_BINARY_DIR}" EXCLUDE_FROM_ALL)
		endif()
	endmacro()
endif()

FetchContent_Declare(
	Quiver
	GIT_REPOSITORY https://github.com/JosuaRieder/Quiver
	GIT_TAG        1b5b5b724894c46669706e752706014a387c9b5f
)
FetchContent_MakeAvailable(Quiver)

add_library(projects.cohomology INTERFACE)
target_include_directories(projects.cohomology INTERFACE "${CMAKE_CURRENT_LIST_DIR}")
target_link_libraries(projects.cohomology INTERFACE
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  Quiver
)
target_compile_features(projects.cohomology INTERFACE cxx_std_20)
set_target_properties(projects.cohomology PROPERTIES FOLDER "projects/cohomology")

option(LF_PROJECTS_COHOMOLOGY_BUILD_EXAMPLE "Whether the example for the cohomology project should be built" OFF)
if(LF_PROJECTS_COHOMOLOGY_BUILD_EXAMPLE)
  add_executable(projects.cohomology.example example.cc)
  target_link_libraries(projects.cohomology.example PUBLIC
    projects.cohomology
  )
  set_target_properties(projects.cohomology.example PROPERTIES FOLDER "projects/cohomology")
endif()

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
