add_executable(distributed_sorter distributed_sorter.cpp ../mpi/environment.cpp)
target_link_libraries(distributed_sorter ${MPI_LIBRARIES} tlx)

add_executable(prefix_doubling prefix_doubling.cpp ../mpi/environment.cpp)
target_link_libraries(prefix_doubling ${MPI_LIBRARIES} tlx)

add_executable(hQuicksort HQuicksort.cpp ../mpi/environment.cpp)
target_link_libraries(hQuicksort ${MPI_LIBRARIES} tlx)

add_executable(writeSuffixes writeSuffixesToFile.cpp ../mpi/environment.cpp)
target_link_libraries(writeSuffixes ${MPI_LIBRARIES} tlx)
