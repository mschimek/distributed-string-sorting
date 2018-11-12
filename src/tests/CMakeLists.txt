
add_executable(sorter_with_LcpStringPtr sorter_with_LcpStringPtr.cpp)
target_link_libraries(sorter_with_LcpStringPtr tlx)

add_executable(allgather_simple_test allgather_simple_test.cpp ../mpi/environment.cpp)
target_link_libraries(allgather_simple_test ${MPI_LIBRARIES} tlx)

add_executable(lcp_string_container lcp_string_container.cpp ../mpi/environment.cpp)
target_link_libraries(lcp_string_container ${MPI_LIBRARIES} tlx)
