#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <random>
#include <set>
#include <mpi.h>
#include <random>

int ProcNum = 0; 
int ProcRank = 0; 
const int Start_Array_Size = 1e7;

std::mutex Mutex;

void Print_Array(int* Array, int Size) {
    for (int i = 0; i < Size; i++) {
        std::cout << Array[i] << ' ';
    }
    std::cout << '\n';
}

bool Start_Array_Initialization(int* Start_Array) {

    for (int i = 0; i < Start_Array_Size; i++) {
        Start_Array[i] = std::rand();
    }
    //Print_Array(Start_Array, 40);

    return true;
}

void Heapify(int* Array, int n, int i)
{
    int Largest = i;
    int L = 2 * i + 1;
    int R = 2 * i + 2;

    if (L < n && Array[L] > Array[Largest])
        Largest = L;

    if (R < n && Array[R] > Array[Largest])
        Largest = R;

    if (Largest != i) {
        std::swap(Array[i], Array[Largest]);
        Heapify(Array, n, Largest);
    }
}

void HeapSort(int *Array, int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
        Heapify(Array, n, i);

    for (int i = n - 1; i > 0; i--) {
        std::swap(Array[0], Array[i]);
        Heapify(Array, i, 0);
    }
}

bool Start_Parallel_Sort(int* Start_Array) {
    int Work_Array_Size = Start_Array_Size / ProcNum;

    //std::cout << "Work Array Size = " << Work_Array_Size << '\n';

    int* Work_Array = Start_Array + ProcRank * Work_Array_Size;

    if (ProcRank == ProcNum - 1) Work_Array_Size =
        Start_Array_Size - Work_Array_Size * (ProcNum - 1);
    HeapSort(Work_Array, Work_Array_Size);

    if (ProcRank == 0 && ProcNum > 1) {
        MPI_Status Status;
        int* Received_Work_Array = Start_Array + Work_Array_Size;
        
        int Total_Received_Values = Work_Array_Size;

        for (int i = 1; i < ProcNum - 1; i ++) {
            MPI_Recv(Received_Work_Array, Work_Array_Size, MPI_INT, i, 0, MPI_COMM_WORLD, &Status);
            Received_Work_Array += Work_Array_Size;
            Total_Received_Values += Work_Array_Size;
        }

        Work_Array_Size = Start_Array_Size - Total_Received_Values;
        MPI_Recv(Received_Work_Array, Work_Array_Size, MPI_INT, ProcNum - 1, 0, MPI_COMM_WORLD, &Status);
    }
    else if (ProcRank > 0 && ProcNum > 1) {
        MPI_Send(Work_Array, Work_Array_Size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    return true;
}

bool Merge_Work_Arrays(int* Start_Array) {
    if (ProcRank) return false;

    int* New_Array = new int[Start_Array_Size];
    int* Work_Arrays_Starts = new int[ProcNum];
    int* Work_Arrays_Ends = new int[ProcNum];

    int Work_Array_Size = Start_Array_Size / ProcNum;

    for (int i = 0; i < ProcNum - 1; i++) {
        Work_Arrays_Starts[i] = Work_Array_Size * i;
    }
    Work_Arrays_Starts[ProcNum - 1] = Start_Array_Size - (Start_Array_Size - Work_Array_Size * (ProcNum - 1));

    for (int i = 0; i < ProcNum - 1; i++) {
        Work_Arrays_Ends[i] = Work_Array_Size * (i + 1);
    }
    Work_Arrays_Ends[ProcNum - 1] = Start_Array_Size;

    Work_Arrays_Starts[0] = 0;

    for (int i = 0; i < Start_Array_Size; i++) {
        int Min_Value = INT_MAX;
        int Min_Value_Array_Index = 0;

        for (int j = 0; j < ProcNum; j++) {
            if (Work_Arrays_Starts[j] < Work_Arrays_Ends[j]
                && Min_Value > Start_Array[Work_Arrays_Starts[j]]) {
                Min_Value = Start_Array[Work_Arrays_Starts[j]];
                Min_Value_Array_Index = j;
            }
        }

        New_Array[i] = Min_Value;
        Work_Arrays_Starts[Min_Value_Array_Index] ++;
    }
    
    for (int i = 0; i < Start_Array_Size; i++) {
        Start_Array[i] = New_Array[i];
    }
  
    delete[] New_Array;
}

int main(int argc, char* argv[]) {
    
    int* Start_Array = new int[Start_Array_Size];;
    int* Array_for_Test = nullptr;

    double Start, Finish, Duration;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) Start_Array_Initialization(Start_Array);
    if (ProcRank == 0) {
        Array_for_Test = new int[Start_Array_Size];
        for (int i = 0; i < Start_Array_Size; i++) {
            Array_for_Test[i] = Start_Array[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
   
    Start = MPI_Wtime();
    MPI_Bcast(Start_Array, Start_Array_Size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    Start_Parallel_Sort(Start_Array);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (ProcRank == 0 && ProcNum > 1) {
        Merge_Work_Arrays(Start_Array);
    }

    Finish = MPI_Wtime();
    Duration = Finish - Start;

    if (ProcRank == 0) {
        printf("Time of execution = %f\n", Duration);
    }

    //if (ProcRank == 0) Print_Array(Start_Array, 40);

    MPI_Finalize();

    HeapSort(Array_for_Test, Start_Array_Size);

}
