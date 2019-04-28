add_executable(distributed_sorter distributed_sorter.cpp ../mpi/environment.cpp)
target_link_libraries(distributed_sorter ${MPI_LIBRARIES} tlx rbcsorting atomic)

add_executable(prefix_doubling prefix_doubling.cpp ../mpi/environment.cpp)
target_link_libraries(prefix_doubling ${MPI_LIBRARIES} tlx)
