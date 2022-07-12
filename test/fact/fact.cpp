#include <iostream>
#include <string>

#include <mth_util.h>

template <typename T>
void print_fact(T n, T d)
{

  using namespace mth_util;
  std::cout << '\n' << "factorial(" << n << ")    = " << factorial(n);
  std::cout << "\n---";
  std::cout << '\n' << "factorial(" << n << "/"  << d <<
                                            ")  = " << factorial(n)/
                                                         factorial(d);
  std::cout << '\n' << "factorial(" << n << ", " << d <<
                                            ") = " << factorial(n, d);
}

/*
 * Quick test program for the factorial functions
 *
 * $ fact n d
 *
 * n!/d! is computed using int, long, and double
 */
int main(int argc, char* argv[])
{

  if (argc != 3) {
    std::cerr << "\nProper use is:  " << argv[0] << " n d\n";
    return 0;
  }

  {
    int n {std::stoi(argv[1])};
    int d {std::stoi(argv[2])};
    print_fact(n, d);
  }
  std::cout << '\n';
  {
    long n {std::stol(argv[1])};
    long d {std::stol(argv[2])};
    print_fact(n, d);
  }
  std::cout << '\n';
  {
    double n {std::stod(argv[1])};
    double d {std::stod(argv[2])};
    print_fact(n, d);
  }

  std::cout << '\n';
}
