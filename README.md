# Lab1
```
Student Name: Le Duong
Student ID: 8422234186
Email: lnduong@usc.edu
```

Program Summary:
```
My code reads two square matrices from `input.txt` and dynamically allocates the arrays needed to store them. If the matrix size is not a power of two, it pads the matrices with zeros so the recursive algorithms can divide them evenly. It then multiplies the matrices using both the divide and conquer method and Strassen’s method and records the runtime of each algorithm. Finally, it writes the two resulting matrices and runtimes to their corresponding output files and deallocates all the memory that was created.
```

For the power_of_two function, I asked google AI for some methods to 
find the next integer that's a power of 2 since that's a shape
requirement for the matrices.

Prompt:
```
What are some methods to find the next integer that's a power of 2?
```

Output:
```
** It output man things but this method struck me as the most simple and
efficient way to get the next power of 2, I only changed the name
of the function: **

unsigned int nextPowerOfTwo(unsigned int n) {
    if (n == 0) return 1;
    
    n--; // Decrement to handle cases where n is already a power of 2
    n |= n >> 1;  // Flood right bits with 1s
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16; // Handles up to 32-bit integers
    
    return n + 1; // Adding 1 rolls it over to the next power of 2
}
```

For the divide_and_conquer recursive function, I asked ChatGPT to give me
a sudo code for the algorithm that would work for a dynamically 
allocated array. This ensures that I don't have to create new A and B 
sub-matrices every iteration. I gave it a screenshot of the problem (since
there are multiple recursive divide and conquer matrix mult. algorithms out 
there, strassen being one of them) to define the algorithm. Then I wrote the 
code based on the sudo code it gave me. 

This is my prompt:

```
"give me the sudo code for the divide and conquer recursive algorithm that 
would work for a dynamically allocated 2D matrix array"
```

ChatGPT's output:

```
"PROCEDURE RecursiveMultiply(
    A, aRow, aCol,
    B, bRow, bCol,
    C, cRow, cCol,
    blockSize
):

    // Base case: the blocks are 1 × 1
    IF blockSize == 1:
        C[cRow][cCol] =
            C[cRow][cCol] + A[aRow][aCol] × B[bRow][bCol]
        RETURN

    half = blockSize / 2

    // C11 = A11B11 + A12B21
    RecursiveMultiply(
        A, aRow,        aCol,
        B, bRow,        bCol,
        C, cRow,        cCol,
        half
    )

    RecursiveMultiply(
        A, aRow,        aCol + half,
        B, bRow + half, bCol,
        C, cRow,        cCol,
        half
    )

    // C12 = A11B12 + A12B22
    RecursiveMultiply(
        A, aRow,        aCol,
        B, bRow,        bCol + half,
        C, cRow,        cCol + half,
        half
    )

    RecursiveMultiply(
        A, aRow,        aCol + half,
        B, bRow + half, bCol + half,
        C, cRow,        cCol + half,
        half
    )

    // C21 = A21B11 + A22B21
    RecursiveMultiply(
        A, aRow + half, aCol,
        B, bRow,        bCol,
        C, cRow + half, cCol,
        half
    )

    RecursiveMultiply(
        A, aRow + half, aCol + half,
        B, bRow + half, bCol,
        C, cRow + half, cCol,
        half
    )

    // C22 = A21B12 + A22B22
    RecursiveMultiply(
        A, aRow + half, aCol,
        B, bRow,        bCol + half,
        C, cRow + half, cCol + half,
        half
    )

    RecursiveMultiply(
        A, aRow + half, aCol + half,
        B, bRow + half, bCol + half,
        C, cRow + half, cCol + half,
        half
    )"
```

For the function call, I asked ChatGPT where I should begin. 

Prompt:
```
"how would I call the sudo function? What rows and columns should I input?"
```
ChatGPT's output:
```
"For the initial call, all row and column positions should be 0 because you are multiplying the entire matrices starting at their top-left corners.

DivideAndConquer(
    matrix_A, 0, 0,       // A starts at row 0, column 0
    matrix_B, 0, 0,       // B starts at row 0, column 0
    matrix_C1, 0, 0,      // Store result starting at row 0, column 0
    blockSize            // Size of the matrices
);"
```

For the strassen function, I did the same thing as DNQ, where I asked ChatGPT for a sudo code, 
then write the c++ code on my own. I uploaded the strassen method as described in the lab assignment

Prompt:
```
give me the sudo code for this strassen method
```

Output:
```
PROCEDURE Strassen(
    A, aRow, aCol,
    B, bRow, bCol,
    C, cRow, cCol,
    blockSize
):

    // Base case: multiply two 1 × 1 blocks
    IF blockSize == 1:
        C[cRow][cCol] = A[aRow][aCol] × B[bRow][bCol]
        RETURN

    half = blockSize / 2

    ALLOCATE temporary matrices T1 and T2 of size half × half
    ALLOCATE matrices M1, M2, M3, M4, M5, M6, M7
        of size half × half

    // M1 = (A11 + A22)(B11 + B22)
    T1 = A11 + A22
    T2 = B11 + B22
    Strassen(T1, 0, 0, T2, 0, 0, M1, 0, 0, half)

    // M2 = (A21 + A22)B11
    T1 = A21 + A22
    Strassen(T1, 0, 0, B, bRow, bCol, M2, 0, 0, half)

    // M3 = A11(B12 - B22)
    T2 = B12 - B22
    Strassen(A, aRow, aCol, T2, 0, 0, M3, 0, 0, half)

    // M4 = A22(B21 - B11)
    T2 = B21 - B11
    Strassen(
        A, aRow + half, aCol + half,
        T2, 0, 0,
        M4, 0, 0,
        half
    )

    // M5 = (A11 + A12)B22
    T1 = A11 + A12
    Strassen(
        T1, 0, 0,
        B, bRow + half, bCol + half,
        M5, 0, 0,
        half
    )

    // M6 = (A21 - A11)(B11 + B12)
    T1 = A21 - A11
    T2 = B11 + B12
    Strassen(T1, 0, 0, T2, 0, 0, M6, 0, 0, half)

    // M7 = (A12 - A22)(B21 + B22)
    T1 = A12 - A22
    T2 = B21 + B22
    Strassen(T1, 0, 0, T2, 0, 0, M7, 0, 0, half)

    // Combine the seven products
    FOR row = 0 TO half - 1:
        FOR column = 0 TO half - 1:

            C[cRow + row][cCol + column] =
                M1[row][column]
                + M4[row][column]
                - M5[row][column]
                + M7[row][column]

            C[cRow + row][cCol + half + column] =
                M3[row][column]
                + M5[row][column]

            C[cRow + half + row][cCol + column] =
                M2[row][column]
                + M4[row][column]

            C[cRow + half + row][cCol + half + column] =
                M1[row][column]
                - M2[row][column]
                + M3[row][column]
                + M6[row][column]

    DEALLOCATE T1, T2, and M1 through M7
```

Line 377, I asked chatGPT how to convert the running time into integer. I gave it
a picture of the sample code from the lab guide. 

Prompt:
```
convert the running_time into an integer
```

Output:
```
long long running_time =
    static_cast<long long>(end - start);
```

