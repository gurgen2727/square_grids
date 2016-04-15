#ifndef _FUNCTIONS_HPP_
#define _FUNCTIONS_HPP_

#include <vector>

//print hamiltonian loop in grid
void print_vector_in_grid(std::vector<int> &, int);

//distance between two cycles 
int distance(int, int);

//print hamiltonian loop
void print_vector(std::vector<int> &);

//generate hamiltonian cycles
//int generate_ham_cycle_helper(int, int&, std::vector<int> &);
int generate_ham_cycle_helper(int, int&, std::vector<int> &, int** &, const bool* &, int&, int&);

//generate non simetric matrix for mxn grid
void generate_grid_matrix(int** &, int, int);

//visit generate_ham_cycle_helper
void generate_ham_cycle(int, int, int = 1);


#endif//_FUNCTIONS_HPP_
