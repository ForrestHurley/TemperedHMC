
cmake_minimum_required(VERSION 2.8)
project(hamiltonian_monte_carlo)

set(HAMILTONIAN_MONTE_CARLO_SOURCES src/mcmc_diagnostics.cc)

add_library(hamiltonian_monte_carlo ${HAMILTONIAN_MONTE_CARLO_SOURCES})

target_include_directories(hamiltonian_monte_carlo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
