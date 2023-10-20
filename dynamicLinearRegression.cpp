#include <iostream>
#include <lapacke.h>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;
#define nrhs 1

int main() {
    int MAX_COLS;
    cout << "Enter the dimension of A: ";
    cin >> MAX_COLS;
    int NN=MAX_COLS-3;

    int MAX_ROWS;
    cout << "Enter the number of data's rows and data's columns in order: ";
    cin >> MAX_ROWS;

    int info;
    int LDA = NN; // LDA>=max(1,NN)
    float test = 0.2;
    float train = 1.0 - test;
    int LDB = NN;  // LDB>=max(1,NN)
    int N = static_cast<int>(MAX_ROWS * train);
    int NdarNplus1Dovom = (NN*(NN+1))/2;
    int ipiv[NN];
    srand(time(NULL));

    // Declare a memory block of size MAX_COLS for A[i] B[i] ...
    double** dataSet = new double*[N];
for (int i = 0; i < N; i++) {
    dataSet[i] = new double[MAX_COLS];
    for (int j = 0; j < MAX_COLS; j++) {
        dataSet[i][j] = 0.0;
    }
}
    double** data = new double*[MAX_ROWS];
    for (int i = 0; i < MAX_ROWS; i++) {
        data[i] = new double[MAX_COLS];
        for (int j = 0; j < MAX_COLS; j++) {
            data[i][j] = 0.0;
        }
    }

    ifstream readdata("Book2.txt");
    for (int i = 0; i < MAX_ROWS; i++) {
        for (int j = 0; j < MAX_COLS; j++) {
            readdata >> data[i][j];
        }
    }
    double dummy;
    int random1, random2;
    for (int i = 0; i <  MAX_ROWS; i++) { //there was a 50*MAX_ROWS, I erased the 50
        random1 = rand() % MAX_ROWS;
        random2 = rand() % MAX_ROWS;
        for (int j = 0  ; j < MAX_COLS; j++) {
            dummy = data[random1][j];
            data[random1][j] = data[random2][j];
            data[random2][j] = dummy;
        }
    }

    ofstream out("output2.txt");

    double* avDataSet = new double[MAX_COLS];
        for (int i = 0; i < MAX_COLS; i++) {
            avDataSet[i] = 0.0;
        }
        for (int j = 0; j < MAX_COLS; j++) {
            for (int i = 0; i < N; i++) {
                dataSet[i][j] = data[i][j];
                avDataSet[j] += dataSet[i][j];
            }
            avDataSet[j]/=N;
        }
//for (int j = 0; j < MAX_COLS; j++)cout<<avDataSet[j]<<'\n';
    for(int i=0;i<N;i++){
        for(int j=0;j<(MAX_COLS);j++){
            dataSet[i][j]-=avDataSet[j];
            }
        }


        // Dynamic memory allocation of Left hand(A) and right hand(B) matrices

    //defining Left hands: avX2, avXY,avY2 and so on



double* MLeftHand = new double[NN*NN];
for (int j = 0; j < NN*NN; j++) {
    MLeftHand[j]=0.0;
}
for (int i = 0; i < NN; i++) {
    for(int j = i; j < NN; j++){
        for (int k = 0; k < N; k++) {
            MLeftHand[NN*i+j] += (dataSet[k][j] * dataSet[k][i]);//adds 1 if j>0

        }
    MLeftHand[NN*i+j] /= (N*1e7);
     if(j!=i) MLeftHand[NN*j+i]= MLeftHand[NN*i+j];
     cout<<MLeftHand[NN*i+j]<<'\n';
    }

}
//MLeftHand[0]=0;
//MLeftHand[1]=1;
//MLeftHand[2]=1;
//MLeftHand[3]=0;




//for (int i = 0; i < NN; i++) {
//    for(int j = i; j < NN; j++){
//        if(j==i)MLeftHand[NN*i+j]=1;else MLeftHand[NN*i+j]=0;
//    }
//}



// Define Right hands: avXK, avYK, avZK

double* MRightHand = new double[NN];


for (int j = 0; j < NN; j++) {
    MRightHand[j] = 0.0;
    for (int i = 0; i < N; i++) {
        MRightHand[j] += (dataSet[i][j] * dataSet[i][MAX_COLS-1]);
    }
    MRightHand[j] /= (N*1e7);

    cout<<MRightHand[j]<<'\n';
}

    info=LAPACKE_dgesv(LAPACK_COL_MAJOR,NN, nrhs, MLeftHand, LDA, ipiv, MRightHand , LDB);


    if (info > 0)
    {
        cout << "The diagonal element of the triangular factor of A, ";
        cout << "U(" << info << "," << info << ") is zero, so the ";
        cout << "solution could not be computed." << endl;
        return 1;
    }
    else if (info < 0)
    {
        cout << "Argument " << -info << " had an illegal value." << endl;
        return 1;
    }
    else
    {
    cout << "The solution is:" << endl;
    for (int i = 0; i < NN; i++)
    {
        cout << "x[" << i << "] = " << MRightHand[i] << endl;
    }

    }


            for (int i = N; i < MAX_ROWS; i++) {
                cout<<data[i][MAX_COLS-1]-avDataSet[MAX_COLS-1]<<'\t'<<  MRightHand[0]*(data[i][0]-avDataSet[0])+     MRightHand[1] *(data[i][1]-avDataSet[1])<<'\n';
            }


    for(int i=0;i<N;i++)
        out<<dataSet[i][0]<<'\t'<<dataSet[i][1]<<'\t'<<dataSet[i][2]<<'\n'<<MLeftHand[0] <<'\t'<< MRightHand[i]<< '\t' << MLeftHand[i] << endl;
    out.close();
    double ave=0.0;
    cout<<"actual price" << '\t' << "our pre" << '\t'<< '\t'<< "fabs" << '\n';
    for(int j=0; j< MAX_COLS;j++){
        for(int i=N+1;i<MAX_ROWS;i++){
            cout << data[i][4]<<'\t'<<'\t' <<MRightHand[0]*data[i][0]+MRightHand[1]*data[i][3]<<'\t'<< '\t'<<fabs(-data[i][4]+MRightHand[0]*data[i][0]+MRightHand[1]*data[i][3])<<'\n';
            ave+=pow((-data[i][4]+MRightHand[0]*data[i][0]+MRightHand[1]*data[i][3]),2.0);
    }
}
    ave/=(MAX_ROWS-N-2);
    ave=sqrt(ave);
    cout<<"sigma= "<<ave<<'\n';

    // Clean up dynamically allocated memory
   for (int i = 0; i < N; i++) {
        delete[] dataSet[i];
    }
    delete[] dataSet;

    delete[] MLeftHand;
    delete[] MRightHand;

    for (int i = 0; i < MAX_ROWS; i++) {
        delete[] data[i];
    }
    delete[] data;

    delete[] avDataSet;

    return 0;
}
