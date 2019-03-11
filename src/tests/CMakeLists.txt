add_executable(alltoallStringsTests alltoallStrings_tests.cpp ../mpi/environment.cpp)
target_link_libraries(alltoallStringsTests ${MPI_LIBRARIES} tlx)

add_executable(bloomfilterTests bloomfilter_tests.cpp ../mpi/environment.cpp)
target_link_libraries(bloomfilterTests ${MPI_LIBRARIES} tlx)

add_executable(dnGeneratorTests DNGenerator_tests.cpp ../mpi/environment.cpp)
target_link_libraries(dnGeneratorTests ${MPI_LIBRARIES} tlx)

add_executable(measuringToolTests measuringTool_tests.cpp ../mpi/environment.cpp)
target_link_libraries(measuringToolTests ${MPI_LIBRARIES} tlx)

add_executable(pipelinedGolombEncodingTests pipelinedGolombEncoding_tests.cpp ../mpi/environment.cpp)
target_link_libraries(pipelinedGolombEncodingTests ${MPI_LIBRARIES} tlx)
