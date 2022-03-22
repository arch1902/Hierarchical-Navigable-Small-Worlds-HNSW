#include <string>
#include <mpi.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
#include <string>
#include <omp.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]){
    ifstream fs;

    int max_level;
    fs.open("bin_inputs/max_level", ios::binary | ios::in);
    fs.read((char *)&max_level, sizeof(max_level));
    fs.close();
    cout << max_level <<'\n';

    int ep;
    fs.open("bin_inputs/ep", ios::binary | ios::in);
    fs.read((char *)&ep, sizeof(ep));
    fs.close();
    cout << ep <<'\n';    

    vector<int> level;
    fs.open("bin_inputs/level",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    int file_size = fs.tellg()/4; 
    cout << file_size <<'\n';
    fs.close();
    fs.open("bin_inputs/level",ios::binary | ios::in); 
    if (fs.is_open()){
        int tp;
        for (int j=0;j<file_size;j++){
            fs.read((char *)&tp, sizeof(tp));
            level.push_back(tp);
        }
        fs.close(); 
    }
    cout << level.size()<< "  "<<level[file_size-1] <<'\n';

    vector<int> index;
    fs.open("bin_inputs/index",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    file_size = fs.tellg()/4; 
    cout << file_size <<'\n';
    fs.close();
    fs.open("bin_inputs/index",ios::binary | ios::in);  
    if (fs.is_open()){
        int tp;
        for (int j=0;j<file_size;j++){
            fs.read((char *)&tp, sizeof(tp));
            index.push_back(tp);
        }
        fs.close(); 
    }
    cout << index.size()<< "  "<<index[file_size-1] <<'\n';


    vector<int> indptr;
    fs.open("bin_inputs/indptr",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    file_size = fs.tellg()/4; 
    cout << file_size <<'\n';
    fs.close();
    fs.open("bin_inputs/indptr",ios::binary | ios::in);
    if (fs.is_open()){
        int tp;
        for (int j=0;j<file_size;j++){
            fs.read((char *)&tp, sizeof(tp));
            indptr.push_back(tp);
        }
        fs.close(); 
    }  
    cout << indptr.size()<< "  "<<indptr[file_size-1] <<'\n';

    vector<int> level_offset;
    fs.open("bin_inputs/level_offset",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    file_size = fs.tellg()/4; 
    cout << file_size <<'\n';
    fs.close();
    fs.open("bin_inputs/level_offset",ios::binary | ios::in);
    if (fs.is_open()){
        int tp;
        for (int j=0;j<file_size;j++){
            fs.read((char *)&tp, sizeof(tp));
            level_offset.push_back(tp);
        }
        fs.close(); 
    }  
    cout << level_offset.size()<< "  "<<level_offset[file_size-1] <<'\n';

    fs.open("bin_inputs/arguments",ios::binary | ios::in); 
    int L, D, U;
    fo.read((char*)&U, sizeof U);
    fo.read((char*)&L, sizeof L);
    fo.read((char*)&D, sizeof D);
    fo.close();

    int rank, size;
    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



    rbuf = (int *)malloc(gsize*100*sizeof(int)); 
    MPI_Gather( sendarray, 100, MPI_INT, rbuf, 100, MPI_INT, root, comm);

    MPI_Finalize();


}
