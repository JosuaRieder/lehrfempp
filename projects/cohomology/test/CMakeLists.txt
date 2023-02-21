include(GoogleTest)

set(src
  test.cc
)
add_executable(projects.cohomology.test "${src}")
target_link_libraries(projects.cohomology.test PUBLIC
  projects.cohomology
  GTest::gtest_main
)
set_target_properties(projects.cohomology.test PROPERTIES FOLDER "projects/cohomology")
gtest_discover_tests(projects.cohomology.test
  SOURCES "${src}"
)
