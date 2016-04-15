//Header from standard libraries
#include <fstream>

//size of blocks in grid
#define BLOCK_SIZE 16

//kernel function  to multiply c=a*b
__global__ void Muld(float* a,float* b,float* c,int n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	float sum = 0;

	if(i>=n || j>=n){
		return;
	}

	for(int k = 0 ; k < n ; ++k) {
		sum = sum + a[i*n + k]*b[k*n + j];
	}

	c[i*n+j] = sum;
}

//host function  to multiply C=A*B
void Mul(float* A, float* B, int n,float* C)
{
        int size;
        // Load A and B to the device
        float* Ad;
        size = n * n * sizeof(float);
        cudaMalloc((void**)&Ad, size);
        cudaMemcpy(Ad, A, size, cudaMemcpyHostToDevice);
        float* Bd;
        size = n * n * sizeof(float);
        cudaMalloc((void**)&Bd, size);
        cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice);
        // Allocate C on the device
        float* Cd;
        size = n * n * sizeof(float);
        cudaMalloc((void**)&Cd, size);
        // Compute the execution configuration assuming
        // the matrix dimensions are multiples of BLOCK_SIZE
        dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
        dim3 dimGrid(n / dimBlock.x + (n%dimBlock.x!=0), n / dimBlock.y+(n%dimBlock.y!=0));
        // Launch the device computation
        Muld<<<dimGrid, dimBlock>>>(Ad, Bd, Cd , n);
        // Read C from the device
        cudaMemcpy(C, Cd, size, cudaMemcpyDeviceToHost);
        // Free device memory
        cudaFree(Ad);
        cudaFree(Bd);
        cudaFree(Cd);
}

int main(int argc,char* argv[])
{
	std::ifstream  fin(argv[1]);
	std::ofstream fout("out.txt");
	int n;
	fin >> n;
        float *hA = (float*)malloc(sizeof(float)*n*n);
        float *hB = (float*)malloc(sizeof(float)*n*n);
        for(int i = 0; i < n; ++i) {
                for(int j = 0; j < n; ++j) {
                        fin >> hA[i*n+j];
                }
        }
        for(int i = 0; i < n; ++i) {
                for(int j = 0; j < n; ++j) {
                        hB[i*n+j] = 0;
                }
		hB[i*n+i] = 1;
        }
	// in case, when count of input matrix > 1
	// hC = hA * hB;
	// Mul(hA,hB,n,hC);
	for(int i = 2; i <= n; i <<=1) {
		Mul(hA,hA,n,hA);
		if(i & n) {
			Mul(hA,hB,n,hB);
		}
	}
        for(int i = 0; i < n; ++i) {
                for(int j = 0; j < n; ++j) {
                        fout << hB[i*n+j] << ' ';
		}
		fout << std::endl;
        }
	free(hA);
	free(hB);
        return 0;
}

