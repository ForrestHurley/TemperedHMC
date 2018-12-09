
cmake_minimum_required(VERSION 2.8)

project(main_tests)


set(MAIN_TESTS_SOURCES exoplanet_test.cc)

add_executable(main_tests ${MAIN_TESTS_SOURCES})

#target_include_directories(main_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(main_tests gmock_main main_lib)
