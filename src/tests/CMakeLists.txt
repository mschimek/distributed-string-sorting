
add_executable(allgather_simple_test allgather_simple_test.cpp ../mpi/environment.cpp)
target_link_libraries(allgather_simple_test ${MPI_LIBRARIES} tlx)

add_executable(lcp_string_container lcp_string_container.cpp ../mpi/environment.cpp)
target_link_libraries(lcp_string_container ${MPI_LIBRARIES} tlx)

add_executable(shift_test shift_test.cpp ../mpi/environment.cpp)
target_link_libraries(shift_test ${MPI_LIBRARIES} tlx)

add_executable(distributed_sorter distributed_sorter.cpp ../mpi/environment.cpp)
target_link_libraries(distributed_sorter ${MPI_LIBRARIES} tlx)

add_executable(local_sorter local_sorter.cpp ../mpi/environment.cpp)
target_link_libraries(local_sorter ${MPI_LIBRARIES} tlx)

add_executable(is_sorted is_sorted.cpp ../mpi/environment.cpp)
target_link_libraries(is_sorted ${MPI_LIBRARIES} tlx)

add_executable(losertree losertree.cpp ../mpi/environment.cpp)
target_link_libraries(losertree ${MPI_LIBRARIES} tlx)
