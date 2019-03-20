#include "mpi/readInput.hpp"
#include "mpi/synchron.hpp"
#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"

int main() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dsss::mpi::environment env;
    const std::string path = "smallTest.dat";
    //auto[chars, blabla] = readFile(path);
    //std::cout << chars.size() << std::endl;
    //for (size_t i = 0; i < chars.size(); ++i) {
    //  std::cout << chars.data() + i << std::endl;
    //  while(i < chars.size() &&  chars[i] != 0)
    //    ++i;
    //}
    std::vector<unsigned char> rawStrings = readFileAndDistribute(path);
    mpi::execute_in_order([&]() {
        StringLcpContainer<StringSet> cont(std::move(rawStrings));
        StringSet ss = cont.make_string_set();
        std::cout << "rank: " << env.rank() << std::endl;
        ss.print();
        });

    env.finalize();
}
