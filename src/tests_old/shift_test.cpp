#include "mpi/shift.hpp"

int main()
{
  dsss::mpi::environment env;

  std::vector<unsigned char> send_string{65, 65, 0};
  if (env.rank() == 1)
    send_string = {66,66, 66, 0};
  if (env.rank() == 2)
    send_string = {0};
  if (env.rank() == 3)
    send_string = {67, 67, 67, 67, 0};
  if (env.rank() == 4)
    send_string = {0};

  constexpr bool shift_left = false;
  std::vector<unsigned char> recv = dss_schimek::mpi::shift_string<shift_left>(send_string.data(), env);
  std::cout << "rank: " << env.rank() << " received string: " << reinterpret_cast<char*>(recv.data()) << std::endl;
  env.finalize();
}
