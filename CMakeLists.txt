
cmake_minimum_required(VERSION 2.8)

project(main)

#temporary fix to allow gmock to build. There is
#already a pull request which should make it unnecessary (10-05-18)
add_compile_options(-std=c++11)

add_subdirectory(hamiltonian_monte_carlo)

set(MAIN_SOURCES src/main.cc src/exoplanet_model.cc src/gaussian_model.cc)

add_executable(main ${MAIN_SOURCES})

target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(main hamiltonian_monte_carlo)
