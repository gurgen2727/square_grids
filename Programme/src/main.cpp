#include "functions.hpp"
#include <iostream>

int main()
{
        int M, N;
        std::cout << "Enter count of rows    ";
        std::cin >> M;
        std::cout << "Enter count of columns ";
        std::cin >> N;
        std::cout << std::endl 
                  << " P" << N << " x P" << M << " grid graph  " 
                  << std::endl;
        generate_ham_cycle(M, N, 1);
        return 0;
}
