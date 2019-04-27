# CMake generated Testfile for 
# Source directory: /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC
# Build directory: /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(optimizedcolls "/usr/bin/bash" "/home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC/test/test_optimizedcolls.sh" "/home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC")
set_tests_properties(optimizedcolls PROPERTIES  DEPENDS "test_optimizedcolls" _BACKTRACE_TRIPLES "/home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC/CMakeLists.txt;108;add_test;/home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC/CMakeLists.txt;0;")
subdirs("external/tlx")
