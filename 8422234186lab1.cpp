#include <iostream>
#include <fstream>
#include <time.h>
 
using namespace std;

unsigned int power_of_two(unsigned int n) {
    // This function will either reutnr the same number if it's already a power of 2
    // Or return the next power of 2 otherwise

    if (n == 0) return 1; // edge case since n-- results in overflow (shouldn't be used since size can't be 0)
    
    n--; // Decrement to handle cases where n is already a power of 2
    n |= n >> 1;  // Flood right bits with 1s
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16; // Handles up to 32-bit integers
    
    return n + 1; // Adding 1 rolls it over to the next power of 2
}

void divide_and_conquer(int** A, unsigned int aRow, unsigned int aCol,
                         int** B, unsigned int bRow, unsigned int bCol,
                         int** C, unsigned int cRow, unsigned int cCol,
                         unsigned int n){
    if (n==1){ // base case
        // Once the recursive algorithm has split the matrices into 1x1, 
        // it begins the multiplication
        C[cRow][cCol] = C[cRow][cCol] + A[aRow][aCol]*B[bRow][bCol];
        return;
    }

    unsigned int half = n/2;

    // C11 = A11B11 + A12B21
    divide_and_conquer(
        A, aRow, aCol,
        B, bRow, bCol,
        C, cRow, cCol,
        half);

    divide_and_conquer(
        A, aRow, aCol + half,
        B, bRow + half, bCol,
        C, cRow, cCol,
        half);

    // C12 = A11B12 + A12B22
    divide_and_conquer(
        A, aRow, aCol,
        B, bRow, bCol + half,
        C, cRow, cCol + half,
        half);

    divide_and_conquer(
        A, aRow, aCol + half,
        B, bRow + half, bCol + half,
        C, cRow, cCol + half,
        half);

    // C21 = A21B11 + A22B21
    divide_and_conquer(
        A, aRow + half, aCol,
        B, bRow, bCol,
        C, cRow + half, cCol,
        half);

    divide_and_conquer(
        A, aRow + half, aCol + half,
        B, bRow + half, bCol,
        C, cRow + half, cCol,
        half);

    // C22 = A21B12 + A22B22
    divide_and_conquer(
        A, aRow + half, aCol,
        B, bRow, bCol + half,
        C, cRow + half, cCol + half,
        half);

    divide_and_conquer(
        A, aRow + half, aCol + half,
        B, bRow + half, bCol + half,
        C, cRow + half, cCol + half,
        half);
}

int** allocate_temp_matrix(unsigned int size){
    // This function allocates a temporary size × size matrix.
    // Filled w/ 0's
    int** matrix = new int*[size];
    for (unsigned int r = 0; r < size; r++){
        matrix[r] = new int[size]{};
    }
    return matrix;
}

void add_blocks(int** X, unsigned int xRow, unsigned int xCol,
               int** Y, unsigned int yRow, unsigned int yCol,
               int** result,
               unsigned int blockSize){
    // This function adds 2 blocks
    // result = block of X + block of Y
    for (unsigned int row = 0; row < blockSize; row++){
        for (unsigned int col = 0; col < blockSize; col++){
            result[row][col] = X[xRow + row][xCol + col] + Y[yRow + row][yCol + col];
        }
    }
}

void subtract_blocks(int** X, unsigned int xRow, unsigned int xCol,
                    int** Y, unsigned int yRow, unsigned int yCol,
                    int** result,
                    unsigned int blockSize){
    // This function subtracts 2 blocks
    // result = block of X - block of Y
    for (unsigned int row = 0; row < blockSize; row++){
        for (unsigned int col = 0; col < blockSize; col++){
            result[row][col] = X[xRow + row][xCol + col] - Y[yRow + row][yCol + col];
        }
    }
}

void deallocate_matrix(int** matrix, unsigned int size){
    // Deallocate memory (CRITICAL to avoid memory leaks)
    for (unsigned int i=0; i<size; i++) {
        delete[] matrix[i]; // Free each row array
    }
    delete[] matrix; // Free the array of pointers
}

void strassen(
    int** A, unsigned int aRow, unsigned int aCol,
    int** B, unsigned int bRow, unsigned int bCol,
    int** C, unsigned int cRow, unsigned int cCol,
    unsigned int blockSize)
{
    
    if (blockSize == 1){// Base case: multiply two 1 × 1 blocks.
        C[cRow][cCol] = A[aRow][aCol] * B[bRow][bCol];
        return;
    }

    unsigned int half = blockSize / 2;

    // Temporary matrices used for additions and subtractions.
    int** temporaryA = allocate_temp_matrix(half);
    int** temporaryB = allocate_temp_matrix(half);

    // The seven strassen products.
    int** M1 = allocate_temp_matrix(half);
    int** M2 = allocate_temp_matrix(half);
    int** M3 = allocate_temp_matrix(half);
    int** M4 = allocate_temp_matrix(half);
    int** M5 = allocate_temp_matrix(half);
    int** M6 = allocate_temp_matrix(half);
    int** M7 = allocate_temp_matrix(half);

    // M1 = (A11 + A22)(B11 + B22)
    add_blocks(
        A, aRow, aCol,
        A, aRow + half, aCol + half,
        temporaryA,
        half
    );

    add_blocks(
        B, bRow, bCol,
        B, bRow + half, bCol + half,
        temporaryB,
        half
    );

    strassen(
        temporaryA, 0, 0,
        temporaryB, 0, 0,
        M1, 0, 0,
        half
    );

    // M2 = (A21 + A22)B11
    add_blocks(
        A, aRow + half, aCol,
        A, aRow + half, aCol + half,
        temporaryA,
        half
    );

    strassen(
        temporaryA, 0, 0,
        B, bRow, bCol,
        M2, 0, 0,
        half
    );

    // M3 = A11(B12 - B22)
    subtract_blocks(
        B, bRow, bCol + half,
        B, bRow + half, bCol + half,
        temporaryB,
        half
    );

    strassen(
        A, aRow, aCol,
        temporaryB, 0, 0,
        M3, 0, 0,
        half
    );

    // M4 = A22(B21 - B11)
    subtract_blocks(
        B, bRow + half, bCol,
        B, bRow, bCol,
        temporaryB,
        half
    );

    strassen(
        A, aRow + half, aCol + half,
        temporaryB, 0, 0,
        M4, 0, 0,
        half
    );

    // M5 = (A11 + A12)B22
    add_blocks(
        A, aRow, aCol,
        A, aRow, aCol + half,
        temporaryA,
        half
    );

    strassen(
        temporaryA, 0, 0,
        B, bRow + half, bCol + half,
        M5, 0, 0,
        half
    );

    // M6 = (A21 - A11)(B11 + B12)
    subtract_blocks(
        A, aRow + half, aCol,
        A, aRow, aCol,
        temporaryA,
        half
    );

    add_blocks(
        B, bRow, bCol,
        B, bRow, bCol + half,
        temporaryB,
        half
    );

    strassen(
        temporaryA, 0, 0,
        temporaryB, 0, 0,
        M6, 0, 0,
        half
    );

    // M7 = (A12 - A22)(B21 + B22)
    subtract_blocks(
        A, aRow, aCol + half,
        A, aRow + half, aCol + half,
        temporaryA,
        half
    );

    add_blocks(
        B, bRow + half, bCol,
        B, bRow + half, bCol + half,
        temporaryB,
        half
    );

    strassen(
        temporaryA, 0, 0,
        temporaryB, 0, 0,
        M7, 0, 0,
        half
    );

    // Construct C11, C12, C21, and C22.
    for (unsigned int row = 0; row < half; row++)
    {
        for (unsigned int col = 0; col < half; col++)
        {
            // C11 = M1 + M4 - M5 + M7
            C[cRow + row][cCol + col] = M1[row][col] + M4[row][col] - M5[row][col] + M7[row][col];

            // C12 = M3 + M5
            C[cRow + row][cCol + half + col] = M3[row][col] + M5[row][col];

            // C21 = M2 + M4
            C[cRow + half + row][cCol + col] = M2[row][col] + M4[row][col];

            // C22 = M1 - M2 + M3 + M6
            C[cRow + half + row][cCol + half + col] = M1[row][col] - M2[row][col] + M3[row][col] + M6[row][col];
        }
    }

    // Delete only the matrices created by this function call.
    deallocate_matrix(temporaryA, half);
    deallocate_matrix(temporaryB, half);

    deallocate_matrix(M1, half);
    deallocate_matrix(M2, half);
    deallocate_matrix(M3, half);
    deallocate_matrix(M4, half);
    deallocate_matrix(M5, half);
    deallocate_matrix(M6, half);
    deallocate_matrix(M7, half);
}

int main(){
    // parse size from input.txt
    ifstream input("input.txt");
    unsigned int size;
    input >> size;
    // Assuming size is always a positive integer
    // Check if matrix is 2^n
    // If not power of 2, give the next integer that's a power of 2
    unsigned int padded_size = power_of_two(size);

    int** matrix_A = new int*[padded_size];
    int** matrix_B = new int*[padded_size];
    int** matrix_C1 = new int*[padded_size];
    int** matrix_C2 = new int*[padded_size];
    int next_element;

    // Load Matrix A
    for (unsigned int r=0; r<padded_size; r++){
        matrix_A[r] = new int[padded_size]; // for each row, allocate padded_size columns 
        for (unsigned int c=0; c<padded_size; c++){
            if (c<size && r<size){
                input >> next_element;
                matrix_A[r][c] = next_element;
            }
            else {
                matrix_A[r][c] = 0;
            }
        }
    }

    // Load Matrix B
    for (unsigned int r=0; r<padded_size; r++){
        matrix_B[r] = new int[padded_size]; // for each row, allocate padded_size columns 
        for (unsigned int c=0; c<padded_size; c++){
            if (c<size && r<size){
                input >> next_element;
                matrix_B[r][c] = next_element;
            }
            else {
                matrix_B[r][c] = 0;
            }
        }
    }

    // Allocate Matrices C1 & C2
    for (unsigned int r=0; r<padded_size; r++){
        matrix_C1[r] = new int[padded_size]{}; // for each row, allocate padded_size columns 
        matrix_C2[r] = new int[padded_size]{}; 
    }

    // Implement Divide and Conquer (DNQ)
    // start timer for DNQ
    time_t start_DNQ, end_DNQ;
    start_DNQ = clock();
    divide_and_conquer(matrix_A, 0, 0,
                       matrix_B, 0, 0,
                       matrix_C1, 0, 0,
                       padded_size);
    end_DNQ = clock();
    long long DNQ_runtime = static_cast<long long>(end_DNQ-start_DNQ);
    

    // Implement Strassen
    // start timer for Strassem
    time_t start_Strassen, end_Strassen;
    start_Strassen = clock();
    strassen(matrix_A, 0, 0,
             matrix_B, 0, 0,
             matrix_C2, 0, 0,
             padded_size);
    end_Strassen = clock();
    long long Strassen_runtime = static_cast<long long>(end_Strassen-start_Strassen);

    // Open the three output files.
    ofstream output_m1("output_m1.txt");
    ofstream output_m2("output_m2.txt");
    ofstream output_q3("output_q3.txt");

    // Write only the original, unpadded portion.
    // Output C1 to output_m1.txt
    for (unsigned int r=0; r<size; r++){
        for (unsigned int c=0; c<size; c++){
            output_m1 << matrix_C1[r][c] << " ";
        }
        if(r<size-1){
            output_m1 << '\n';
        }
    }


    // Output C2 to output_m2.txt
    for (unsigned int r=0; r<size; r++){
        for (unsigned int c=0; c<size; c++){
            output_m2 << matrix_C2[r][c] << " ";
        }
        if(r<size-1){
            output_m2 << '\n';
        }
    }

    // Output DNQ runtime followed by Strassen runtime
    output_q3 << DNQ_runtime << " " << Strassen_runtime;

    input.close();
    output_m1.close();
    output_m2.close();
    output_q3.close();

    // Deallocate the matrices
    deallocate_matrix(matrix_A, padded_size);
    deallocate_matrix(matrix_B, padded_size);
    deallocate_matrix(matrix_C1, padded_size);
    deallocate_matrix(matrix_C2, padded_size);
    return 0;
}