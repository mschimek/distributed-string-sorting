add_executable(allTests allTests.cpp ../mpi/environment.cpp)
target_link_libraries(allTests ${MPI_LIBRARIES} tlx)

add_executable(alltoallStringsTests alltoallStrings_tests.cpp ../mpi/environment.cpp)
target_link_libraries(alltoallStringsTests ${MPI_LIBRARIES} tlx)

add_executable(bloomfilterTests bloomfilter_tests.cpp ../mpi/environment.cpp)
target_link_libraries(bloomfilterTests ${MPI_LIBRARIES} tlx)
