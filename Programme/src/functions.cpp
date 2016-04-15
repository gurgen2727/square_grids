#include "functions.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

std::vector< std::vector<int> > cycles;
std::ofstream fout;

//#define M 50
//#define COLUMN 10

//bool visit[M] = {0,};

//bool matrix[M][M]= { 
//1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//characteristic function of hamiltonian cycle
long long unsigned int charac_function(std::vector<int> &mas, std::vector<int> &f, int column)
{
	int index = 0;
	int result = 0;
	int horilength = 0;
	int vertlength = 0;
	int sum = 0;
	int size = mas.size();
	f.clear();
	for(int i = 1; i <= size / 2; ++i){ // ciklum i erkarutyamb
		//std::cout << " f(" << setw(2) << i << ") " << "= ";
		for(int j = 0; j < size; ++j) {	// hertov cikli j-rd elementic i erkarutyamb
			vertlength = fabs((mas[(j + i)%size] / column) - (mas[j] / column));
			horilength = fabs((mas[(j + i)%size] % column) - (mas[j] % column));
			result += vertlength + horilength;
			//std::cout << setw(2) << vertlength + horilength << "+";
		}
		//std::cout << "\b = " << result << ";" << std::endl;
                //fout << ' ' << result << std::endl;
		f.push_back(result);
                sum += result;
		result = 0;
	}
	return sum;
        //std::cout << sum << std::endl;
}
//print hamiltonian loop in grid
void print_vector_in_grid(std::vector<int> &mas, int column)
{
	const char* ch1 = "___";
	char ch2 = '|';
	int size = mas.size();
	int row = size / column;
	std::cout << std::endl << "\t";
	for(int i = 1; i <= size; ++i){
		//std::cout << std::setw(3) << i;
		std::cout << ' ';
		if(i % column == 0) {
			std::cout << std::endl << "\t";
			if( i / column != row ) {
				for(int j = 0; j < column; ++j) {
					//std::cout << std::setw(3);
					for(int k = 0; k < size; ++k) {
						if(mas[k] == i - column + j) {
							if((k > 0 && mas[k-1] == i + j) ||
							   (k < size-1 && mas[k+1] == i + j) ||
							   (k == 0  && mas[size-1] == i + j)) {
								std::cout << ch2 << " ";
							}
							else {
								//std::cout << " ";
                                                                std::cout << "  ";
							}
							break;
						}	
					}
					//std::cout << "  ";
					std::cout << "  ";
				}
			}
			std::cout << std::endl << "\t";
		}
		else {
			//std::cout << ' ';
			for(int j = 0; j < size; ++j) {
				if(mas[j] == i - 1) {
					if((j > 0   && mas[j-1] == i) ||
					   (j < size-1 && mas[j+1] == i)) {
						std::cout << ch1;
					}
					else {
						//std::cout << ' ';
						std::cout << "   ";
					}
					break;
				}
			}
		}
	}
	std::cout << std::endl;
}
//distance between two cycles
int distance(int c1, int c2)
{
	int t = cycles[c1].size() - 1; 
	if(c1 == c2) {
		return t + 1;
	}
	int c = 2;
	for(int i = 1; i < t; ++i){
		for(int j = 1; j < t; ++j){
			if(cycles[c1][i] == cycles[c2][j]){
				if(cycles[c1][i+1] == cycles[c2][j+1] || cycles[c1][i+1] == cycles[c2][j-1]) {
					++c;
				}
				break;
			}
		}
	}
	return c;
}

//print hamiltonian loop
void print_vector(std::vector<int> &myvector)
{
	int t = myvector.size();
        for(int i = 0; i < t; ++i){
                std::cout << std::setw(3) << myvector[i] + 1 << ' ';
        }
        std::cout << std::endl;
}

//print characteristic function of hamiltonian cycle
void print_function(std::vector<int> &f)
{
	int t = f.size();
        for(int i = 0; i < t; ++i){
		std::cout << " f(" << std::setw(2) << i + 1 << ") " << "= " << f[i] << "   ";
        }
        std::cout << std::endl;
}

//print characteristic function of hamiltonian loop to file functions.in for the octave code
void print_function_to_file(std::vector<int> &myvector)
{
	int t = myvector.size();
        for(int i = 0; i < t; ++i){
                fout << std::setw(3) << myvector[i] << ' ';
        }
        fout << std::endl;
}

//Verify function's value for cycles
bool if_equal(std::vector<int> &f1, std::vector<int> &f2)
{
	int t = f1.size();
	bool b = true;
        for(int i = 0; i < t; ++i){
		if(f1[i] != f2[i]){
			b = false;
		}
        }
	return b;
}

//generate hamiltonian cycles
int generate_ham_cycle_helper(int v, int& w, std::vector<int> &myvector, int** &matrix, bool* &visit, int& size, int& column)
{
        //static long long unsigned int cycle_count = 0; // count of hamiltonyan cycles
        if(myvector.size() == size && matrix[myvector.back()][w] == 1){
                //visit print_vector for hamiltonyan cycle
                //print_vector(myvector);
		cycles.push_back(myvector);
                //++cycle_count;
		//std::cout << "--------->  " << cycle_count << "  <---------" << std::endl;
        }
        //visit[v] = true;
        //for(int i = 0; i < size; ++i){             // optimize, not on all size
        //        if(matrix[v][i] && !visit[i]){
        //                myvector.push_back(i);
        //                generate_ham_cycle_helper(i, w, myvector);
        //        }
        //}
        //visit[v] = false;
        //myvector.pop_back();
        visit[v] = true;
	if(v && matrix[v][v - 1] && !visit[v - 1]) {
		myvector.push_back(v - 1);
		//generate_ham_cycle_helper(v - 1, w, myvector);
                generate_ham_cycle_helper(v - 1, w, myvector, matrix, visit, size, column);
		myvector.pop_back();
	} 
	if(v >= column && matrix[v][v - column] && !visit[v - column]) {
			myvector.push_back(v - column);
			//generate_ham_cycle_helper(v - column, w, myvector);
                        generate_ham_cycle_helper(v - column, w, myvector, matrix, visit, size, column);
			myvector.pop_back();
	} 
	if(v != size - 1 && matrix[v][v + 1] && !visit[v + 1]) {
		myvector.push_back(v + 1);
		//generate_ham_cycle_helper(v + 1, w, myvector);
                generate_ham_cycle_helper(v + 1, w, myvector, matrix, visit, size, column);
		myvector.pop_back();
	} 
	if(v + column < size  && matrix[v][v + column] && !visit[v + column]) {
			myvector.push_back(v + column);
			//generate_ham_cycle_helper(v + column, w, myvector);
                        generate_ham_cycle_helper(v + column, w, myvector, matrix, visit, size, column);
			myvector.pop_back();
	} 
        visit[v] = false;
        //return cycle_count; // return count of hamiltonyan cycles
			      // return (cycle_count / 2)
                              // when the matrix is simetrical;
	return 0;
}

// generate non simetric matrix for mxn grid
void generate_grid_matrix(int** &matrix, int m, int n)
{
        int size = m * n;
        matrix = new int* [size];
        for(int i = 0; i < size; ++i){
                matrix[i] = new int [size];
        }
        for(int i = 0; i < size; ++i){
                for(int j = 0; j < size; ++j){
                        matrix[i][j] = 0;
		}
	}
        for(int s, e, j = 0; j < m; ++j){
                s = j * n;
                e = j * n + n - 1;
                for(int i = s; i < e; ++i){
                        matrix[i][i + 1] = 1;
                        matrix[i + 1][i] = 1;
                }
        }
        for(int i = 0; i < size - n; ++i){
                matrix[i][i + n] = 1;
                matrix[i + n][i] = 1;
        }
        matrix[0][n] = 0;       // for non simetrical matrix
}

// visit generate_ham_cycles_helper
void generate_ham_cycle(int m, int n, int v)
{
	int size = m * n;
	int column = n;
        int** matrix;
        generate_grid_matrix(matrix, m, n);
        ///////////////////////////////////////////////////////
	fout.open("cuda/matrix.in");
	fout << "  " << size << std::endl;
        for(int i = 0; i < size; ++i){
                for(int j = 0; j < size; ++j){
                        fout << std::setw(3) << matrix[i][j];
                }
                fout << std::endl;
        }
	fout.close();
	fout.clear();
        ///////////////////////////////////////////////////////
	//print adjacency matrix for grid graph
        //for(int i = 0; i < size; ++i){
        //        for(int j = 0; j < size; ++j){
        //                std::cout << matrix[i][j] << ' ';
        //                //std::cout << matrix[i][j] << ", ";
        //        }
        //        std::cout << std::endl;
        //}
        ///////////////////////////////////////////////////////
        bool* visit; // = {0,};
	visit = new bool[size];
        for(int i = 0; i < size; ++i){
                visit[i] = false;
        }

        int w = v - 1;		    // starting node of hamiltonian cycle
        std::vector<int> myvector;
        myvector.push_back(w);      // start of hamiltonian cycle
        std::cout << std::endl;
        generate_ham_cycle_helper(w, w, myvector, matrix, visit, size, column);
        std::cout << " Count of HAMILTONIAN CYCLES  " << cycles.size() << std::endl << std::endl;
        //std::cout << "\nCount of HAMILTONIAN CYCLES = " 
        //          << generate_ham_cycle_helper(w, w, myvector, matrix, visit, size, column) 
        //          << std::endl << std::endl;
	//print all hamiltonian cycles
	//for(int i = 0; i < cycles.size(); ++i){
	//	print_vector(cycles[i]);
	//}
	int s = cycles.size();
	if(s) {
		int t,c1,c2;
		int min = distance(0, 0);
		long unsigned int sum1 = 0;
		long unsigned int sum2 = 0;
		std::vector<int> f1;
		std::vector<int> f2;
		c1 = 0;
		c2 = 0;
		//variables for non simetrical cycles
		bool b = true;
		unsigned int count = 0;
		fout.open("octave/functions.in");
		for(int i = 0; i < s; ++i){
			//sum and f1 functions of C1 cycle
			sum1 = charac_function(cycles[i], f1, column);
			for(int j = i + 1; j < s; ++j){
				//sum and f2 functions of C2 cycle
				sum2 = charac_function(cycles[j], f2, column);
				if(if_equal(f1, f2)){
					b = false;
				// print cycles, which function are equaly
				//	print_function(f1);
				//	print_vector_in_grid(cycles[i], column);
				//	print_vector_in_grid(cycles[j], column);
				}
				t = distance(i,j);
				if(min > t){
					min = t;
					c1 = i;
					c2 = j;
				}
			}
			if(b) {
			//print non charac. functions to file functions.in
				print_function_to_file(f1);
				++count;
			}
			b = true;
		}
		fout.close();
		fout.clear();
		std::cout << " Count of non simetrical HAMILTONIAN CYCLES  " 
                          << count << std::endl << std::endl;
		std::cout << " nm   " << size << std::endl;
		std::cout << " min distance   " << min << std::endl;
		std::cout << " max distance   " << 2 * size - 2 * min << std::endl;
		//print_vector(cycles[c1]);
		//print_vector(cycles[c2]);
		print_vector_in_grid(cycles[c1], column);
		print_vector_in_grid(cycles[c2], column);
		//for(int i = 0; i < s; ++i){
		//	for(int j = i + 1; j < s; ++j){
		//		t = distance(i,j);
		//		if(min == t){
		//			std::cout << "/////////////////////////////////////////////////////////////" << std::endl;
		//			print_vector_in_grid(cycles[i], column);
		//			print_vector_in_grid(cycles[j], column);
		//		}
		//	}
		//}
	}
	///////////////////////////////////////////////////////
	//std::vector<int> v1;
	//std::vector<int> v2;
	//std::vector<int> f1;
	//std::vector<int> f2;
	//v1.push_back(0);v1.push_back(1);v1.push_back(9);v1.push_back(17);v1.push_back(18);v1.push_back(10);v1.push_back(2);v1.push_back(3);v1.push_back(11);v1.push_back(12);v1.push_back(4);v1.push_back(5);v1.push_back(6);v1.push_back(7);v1.push_back(15);v1.push_back(23);v1.push_back(31);v1.push_back(30);v1.push_back(22);v1.push_back(14);v1.push_back(13);v1.push_back(21);v1.push_back(29);v1.push_back(28);v1.push_back(20);v1.push_back(19);v1.push_back(27);v1.push_back(26);v1.push_back(25);v1.push_back(24);v1.push_back(16);v1.push_back(8);
	//v2.push_back(0);v2.push_back(1);v2.push_back(2);v2.push_back(3);v2.push_back(11);v2.push_back(12);v2.push_back(4);v2.push_back(5);v2.push_back(6);v2.push_back(7);v2.push_back(15);v2.push_back(23);v2.push_back(31);v2.push_back(30);v2.push_back(22);v2.push_back(14);v2.push_back(13);v2.push_back(21);v2.push_back(29);v2.push_back(28);v2.push_back(20);v2.push_back(19);v2.push_back(27);v2.push_back(26);v2.push_back(18);v2.push_back(10);v2.push_back(9);v2.push_back(17);v2.push_back(25);v2.push_back(24);v2.push_back(16);v2.push_back(8);
	//int sum1 = charac_function(v1, f1, 8);
	//int sum2 = charac_function(v2, f2, 8);
	//print_function(f1);
	//print_function(f2);
	//std::cout << "sum1 = " << sum1 << std::endl;
	//std::cout << "sum2 = " << sum2 << std::endl;
	//print_vector_in_grid(v1, 8);
	//print_vector_in_grid(v2, 8);
        ///////////////////////////////////////////////////////
	std::cout << " Generated" << std::endl;
	std::cout << " Adjacency matrix to file                   cuda/matrix.in" << std::endl;
        std::cout << " Characteristic functions of loops to file  octave/functions.in" << std::endl;
	//generate Kirchhorff matrix
        for(int i = 0; i < size; ++i){
                for(int j = 0; j < size; ++j){
                        matrix[i][j] = -matrix[i][j];
                }
        }
	matrix[0][n] = -1;
	matrix[0][0] = 2;
	matrix[n - 1][n - 1] = 2;
	matrix[size - n][size - n] = 2;
	matrix[size - 1][size - 1] = 2;
        for(int i = 1; i < n-1; ++i){
		matrix[i][i] = 3;
		matrix[size - i - 1][size - i - 1] = 3;
        }
        for(int i = 1; i < m - 1; ++i){
		matrix[i * n][i * n] = 3;
		matrix[n + i * n - 1][n + i * n - 1] = 3;
        }
        for(int i = 1; i < m - 1; ++i){
		for(int j = 1; j < n - 1; ++j){
			matrix[i * n + j][i * n + j] = 4;
		}
        }
	fout.open("octave/matrix.in");
        for(int i = 0; i < size; ++i){
                for(int j = 0; j < size; ++j){
                        fout << std::setw(3) << matrix[i][j];
                }
                fout << std::endl;
        }
	fout.close();
	fout.clear();
	std::cout << " Kirchhorff matrix to file                  octave/matrix.in" << std::endl << std::endl;
        ///////////////////////////////////////////////////////
        for(int i = 0; i < size; ++i){
                delete [] matrix[i];
        }
        delete [] matrix;
}
