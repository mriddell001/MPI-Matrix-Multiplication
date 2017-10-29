/*
Filename: matrixmult.c
Author: mriddell
Name: Matthew Riddell
CSU ID: 005489481
Due: 4/5/2017
CSCI 551
Assignment 3
Contents: The contents are the requirements for the Assignment 3, write a C
program to perform the multiplication of two matrices of integers using MPI.
*/

//Notes, introduce stop points in matrix function to greatly reduce time.

#define _BSD_SOURCE
#include <unistd.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/**
 * mm - This function works with the original core to preserve matrix A
 * @param {INT} n - original n
 * @param {INT} a_stop - stop position for matrix A
 * @param {INT} form - loop order
 * @param {INT} a_a - portion of matrix A
 * @param {INT} a_b - portion of matrix B
 * @param {INT} a_c - matrix C
 */
void mm(int n, int a_stop, int form, int *a_a, int *a_b, int *a_c) {
  int x = 0, y = 0, z = 0;
  /*Note, the order of the indexes in array B will be reversed because the array
   has been transposed and the rows and columns are switched. */
  switch (form) {//C[ij] = A[ik]*B[kj]
   case 0: {//ijk
     for (int i = 0; i < a_stop; i++) {
       for (int j = 0; j < n; j++) {
         for (int k = 0; k < n; k++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 1: {//ikj
     for (int i = 0; i < a_stop; i++) {
       for (int k = 0; k < n; k++) {
         for (int j = 0; j < n; j++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 2: {//jik
     for (int j = 0; j < n; j++) {
       for (int i = 0; i < a_stop; i++) {
         for (int k = 0; k < n; k++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 3: {//jki
     for (int j = 0; j < n; j++) {
       for (int k = 0; k < n; k++) {
         for (int i = 0; i < a_stop; i++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 4: {//kij
     for (int k = 0; k < n; k++) {
       for (int i = 0; i < a_stop; i++) {
         for (int j = 0; j < n; j++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 5: {//kji
     for (int k = 0; k < n; k++) {
       for (int j = 0; j < n; j++) {
         for (int i = 0; i < a_stop; i++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   default: {break;}
  }
}

/**
 * mma - This function works with a partial portion of Matrix A
 * @param {INT} n - original n
 * @param {INT} a_stop - stop position for matrix A
 * @param {INT} i_set - start position for matrix A
 * @param {INT} form - loop order
 * @param {INT} a_a - portion of matrix A
 * @param {INT} a_b - portion of matrix B
 * @param {INT} a_c - matrix C
 */
void mma(int n, int a_stop, int i_set, int form, int *a_a, int *a_b, int *a_c) {
  int x = 0, y = 0, z = 0;
  /*Note, the order of the indexes in array B will be reversed because the array
   has been transposed and the rows and columns are switched. */
  switch (form) {//C[ij] = A[ik]*B[kj]
   case 0: {//ijk
     for (int i = i_set; i < a_stop; i++) {
       for (int j = 0; j < n; j++) {
         for (int k = 0; k < n; k++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 1: {//ikj
     for (int i = i_set; i < a_stop; i++) {
       for (int k = 0; k < n; k++) {
         for (int j = 0; j < n; j++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 2: {//jik
     for (int j = 0; j < n; j++) {
       for (int i = i_set; i < a_stop; i++) {
         for (int k = 0; k < n; k++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 3: {//jki
     for (int j = 0; j < n; j++) {
       for (int k = 0; k < n; k++) {
         for (int i = i_set; i < a_stop; i++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 4: {//kij
     for (int k = 0; k < n; k++) {
       for (int i = i_set; i < a_stop; i++) {
         for (int j = 0; j < n; j++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   case 5: {//kji
     for (int k = 0; k < n; k++) {
       for (int j = 0; j < n; j++) {
         for (int i = i_set; i < a_stop; i++) {
           x= *(a_a + i*n + k); y= *(a_b + j*n + k); z= *(a_c + i*n + j);
           *(a_c + i*n + j) = x * y + z;}}}
     break;}
   default: {break;}
  }
}

/**
 * split - This function determines if the split of the data is even
 *         based on the n-value and total number of processes.
 * @param {INT} n - matrix column and rows dimensions
 * @param {INT} c_sz - Number of processes
 * @returns {INT} - returns the split type.
 */
int split(int n, int c_sz) {
  int mod;
  mod = n % c_sz;
  if (mod == 0) {return 0;}
  else {return 1;}
}

/**
 * order - This function determines what case the matrix form type is.
 * @param {CHAR} a - first input to form type.
 * @param {CHAR} b - second input to form type.
 * @param {CHAR} c - third input to form type.
 * @returns {INT} - returns the for loop order type.
 */
int order(char a, char b, char c) {
  int x,y,z,sum=0,form=0;
  x = a - '0'; y = b - '0'; z = c - '0';
  x = x - 57;  y = y - 57;  z = z - 57;
  sum = x*100+y*10+z;
  switch (sum) {
    case 12:{form=0;break;}
    case 21:{form=1;break;}
    case 102:{form=2;break;}
    case 120:{form=3;break;}
    case 201:{form=4;break;}
    case 210:{form=5;break;}
  }
  return form;
}

/**
 * main - This function accepts input from the user using stdin and initializes
 *        variables as needed. It also splits the data into chunks based on the
 *        number of processors and divides the work.
 *        once work is divided, each process will call the function to compute
 *        the math for the matrixes.
 * @returns {INT} - returns value of zero on sucessful exit.
 */
int main() {
  char ch, a, b, c, type;
  int m_rank, c_sz, n, x, tmp=0, form=0, rows, p_rows;
  int *arr_a=NULL, *arr_b=NULL, *arr_c=NULL, *a_start=NULL, *arows=NULL;
  int init_array[4];
  double start, end;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &c_sz);

  //The following processes input to record the order in which the loops will
  //incriment as well as the input mode (Input/Random)
  if (m_rank == 0) {
    for (int i = 0; i < 6; i++) {
      ch = getchar();
      switch (i) {
        case 0: {a = ch; break;}
        case 1: {b = ch; break;}
        case 2: {c = ch; break;}
        case 4: {type = ch; break;}
        default:{break;}
      }
    }
    scanf("%i", &n);
    ch = getchar();
    char r = 'R';                                 //Comparison character
    init_array[0] = n;                            //Side length
    init_array[1] = split(n, c_sz);               //Even split
    init_array[2] = n % c_sz;                     //Ranks < mod receive larger
    init_array[3] = order(a,b,c);                 //Order the loops are executed
    form = init_array[3];                         //Set the form for rank 0
    while ((!arr_a)||(!arr_b)||(!arr_c)) {
      if (!arr_a) {
        arr_a = malloc(sizeof *arr_a * n * n);//The initial matrix to hold A
      }
      if (!arr_b) {
        arr_b = malloc(sizeof *arr_b * n * n);//The initial matrix to hold B
      }
      if (!arr_c) {
        arr_c = malloc(sizeof *arr_c * n * n);//The initial matrix to hold C
      }
    }
    if (type == r) {                              //The input method is random
      time_t t;
      srand((unsigned) time(&t));
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          tmp = rand() % 101;
          *(arr_a + i*n + j) = tmp;
          tmp = rand() % 101;
          *(arr_b + j*n + i) = tmp;
          *(arr_c + i*n + j) = 0;
        }
      }
    }
    else {                                        //The input method is inputed
      for (int i = 0; i < n*n; i++) {             //Record the first matrix: A
        scanf("%i", &tmp);                        //Scan input for integer.
        *(arr_a + i) = tmp;
        *(arr_c + i) = 0;
        ch = getchar();                           //Ignore whitespace: '\n', ' '
      }
      for (int i = 0; i < n; i++) {               //Record the second matrix: B
        for (int j = 0; j < n; j++) {
          scanf("%i", &tmp);                      //Scan input for integer.
          *(arr_b + j*n + i) = tmp;
          ch = getchar();                         //Ignore whitespace: '\n', ' '
        }
      }
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);                   //Start Barrier for MPI timing
  if (m_rank == 0) {start = MPI_Wtime();}

  //If there's only a single processor run multiplication, collect timing, print
  //information, delete arrays and call MPI finalize.
  if (c_sz == 1) {
    mm(n, n, form, arr_a, arr_b, arr_c);
    MPI_Barrier (MPI_COMM_WORLD);
    end = MPI_Wtime();
    double total_time = end - start;
    printf("\nrunning on %i processors.\n", c_sz);
    printf("elapsed time = %.6e seconds\n", total_time);
    char r = 'R';
    if (type != r) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          printf("%d", *(arr_c + i*n + j));
          if (j < n-1) {printf(" ");}
          else {printf("\n");}
        }
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    if (arr_a) {free(arr_a);}
    if (arr_b) {free(arr_b);}
    if (arr_c) {free(arr_c);}
    MPI_Finalize();
    return 0;
  }
  else {
    //This is the case where there's more than one processor working on the
    //matrix multiplication. Matrix B has been transposed, so the rows need to
    //be divided between the processors. If rows cannot be evenly split, lower
    //order processors take extra first until the remainder is all taken care
    //of. Send each processor to the matrix_mult function and add together all
    //the resulting matrixes.
    MPI_Bcast(init_array, 4, MPI_INT, 0, MPI_COMM_WORLD);
    //For processors > 0, set up receiving matrices for incoming matrixes.
    if (m_rank) {
      n = init_array[0];
      form = init_array[3];
      //If n can be split evenly then the work can be distributed evenly.
      if (init_array[1] == 0) {rows = n/c_sz;}
      else {
        // If there is an uneven split, either some processors get zero work, in
        // the case where n<c_sz, or some receive extra, which means that those
        // with rank < n%c_sz get an extra row, and this creates the best split
        // possible given the number of processors and the size of the matrix.
        if (n > c_sz) {
          // If n is larger than the comm size, ex. n=3 and c_sz=2
          rows = n/c_sz;
          if (m_rank < (n%c_sz)) {rows++;}
        }
        else {
          // IF n < c_sz, determine if this process receives work.
          // EX: n = 2 | c_sz = 3 | m_rank = 1
          // Rank 1 receives work.
          // EX: n = 2 | c_sz = 3 | m_rank = 2
          // Rank 2 receives no work.
          if (m_rank < n) {
            rows = 1;
          }
          else {rows = 0;} // This case needs a flag to prevent faulty MPI error
        }
      }
      while ((!arr_a)||(!arr_b)||(!arr_c)) {
        if (!arr_a) {
          arr_a = malloc(sizeof *arr_a * n * n);//The initial matrix to hold A
        }
        if (!arr_b) {
          arr_b = malloc(sizeof *arr_b * n * n);//The initial matrix to hold B
        }
        if (!arr_c) {
          arr_c = malloc(sizeof *arr_c * n * n);//The initial matrix to hold C
        }
      }
      int array_size = n * n;
      for (int i = 0; i < array_size; i++) {
        *(arr_a + i) = 0;
        *(arr_b + i) = 0;
        *(arr_c + i) = 0;
      }
    }
    //Wait until all processors have completed the previous section or it cannot
    //be guarenteed that each processor will have matrices ready to receive.
    MPI_Barrier (MPI_COMM_WORLD);
    //Consider MPI_Bcast for arr_b
    MPI_Bcast(arr_b, n*n, MPI_INT, 0, MPI_COMM_WORLD);
    if (m_rank) {
      p_rows = 0;
      //First determine just how many rows will be received.
      if (m_rank > n) {rows = 0;} //If there's less rows than there are processors.
      else {
        //Spit rows evenly by processors.
        rows = n/c_sz;
        p_rows = m_rank * rows;
        //If the split isn't even.
        if (init_array[1] != 0) {
          //Then the lower rank processors get work first.
          if (m_rank < init_array[2]) {rows++; p_rows = p_rows + m_rank;}
          else {p_rows = p_rows + init_array[2];}
        }
      }
      x = rows;
      MPI_Recv(arr_a + p_rows*n, x*n, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (rows > 0) {
        int received = 1;
        MPI_Send(&received, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        mma(n, rows, p_rows, form, arr_a, arr_b, arr_c);
      }
    }
    else {
      int index = (n/c_sz) + (init_array[1] != 0);
      int *arr_start;
      a_start = malloc (sizeof *a_start * c_sz);
      arows = malloc (sizeof *arows * c_sz);
      *(a_start) = 0;
      int wait_count = 0;
      for (int i = 1; i < c_sz; i++) {
        if (init_array[1] == 0) {rows = n/c_sz;}
        else {
          //This accounts for cases where c_sz > n.
          rows = (n/c_sz) + (i<init_array[2]);
        }
        arr_start = arr_a + index*n;
        *(a_start + i) = *arr_start;
        x = rows;
        *(arows + i) = x;
        index += x;
        MPI_Send(arr_start, x*n, MPI_INT, i, 0, MPI_COMM_WORLD);
        int received = 0;
        MPI_Recv(&received, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        while (!received) {
          wait_count = wait_count + 1;
          if (wait_count > 2147483646) {
            printf("Waiting on process %d to confirm received\n", i);
            wait_count = 0;
          }
        }
        wait_count = 0;
      }
      //Now set the remaining that rank 0 will cover.
      if (init_array[1] == 0) {rows = n/c_sz;}
      else {rows = (n/c_sz) + 1;} //This accounts for cases where c_sz > n.

      if (rows > 0) {mm(n, rows, form, arr_a, arr_b, arr_c);}
    }
    MPI_Barrier (MPI_COMM_WORLD);
    if (m_rank != 0) {
      MPI_Send(arr_c + p_rows*n, x*n, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else {
      for (int p = 1; p < c_sz; p++) {
        MPI_Recv(arr_c + *(a_start + p)*n, *(arows + p)*n, MPI_INT, p, 0,
          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    //At this point each process that isn't zero has completed its work.
    MPI_Barrier (MPI_COMM_WORLD);
    if (m_rank == 0) {
      //If the rank is zero, check if results need to be printed, collect the
      //time and print needed information.
      end = MPI_Wtime();
      double total_time = end - start;
      printf("\nrunning on %i processors.\n", c_sz);
      printf("elapsed time = %.6e seconds\n", total_time);
      char r = 'R';
      if (type != r) {
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            printf("%d", *(arr_c + i*n + j));
            if (j < n-1) {printf(" ");}
            else {printf("\n");}
          }
        }
      }
    }
    //At this point, m_rank 0 should have received all matrix data and it
    //should be ready to finish the closing pieces and finish the MPI finalization.
    MPI_Barrier (MPI_COMM_WORLD);
    if (arr_a) {free(arr_a);}
    if (arr_b) {free(arr_b);}
    if (arr_c) {free(arr_c);}
    if (a_start) {free(a_start);}
    if (arows) {free(arows);}
    MPI_Barrier (MPI_COMM_WORLD);
    if (!m_rank) {MPI_Finalize();}
    return 0;
  }
}
