#include <string>
#include <mpi.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
#include <string>
#include <omp.h>
#include <iostream>

using namespace std;

/*

Path of the to_students folder from personal scratch
../../../../../scratch/cse/phd/anz198717/TA/COL380/A3/to_students

*/

int L, D, U;

double cosine_dist(vector<double> &U, double* V)
{
    double dotp = 0.0, norm_u = 0.0, norm_v = 0.0 ;
    for(int i = 0; i < D; i++) {
        dotp += U[i] * V[i] ;
        norm_u += U[i] * U[i] ;
        norm_v += V[i] * V[i] ;
    }
    return 1 - dotp/(sqrt(norm_u*norm_v)) ;
}

void SearchLayer(vector<double> &q,priority_queue<pair<double,int>> &topk, vector<int> &indptr, vector<int> &index, 
                    vector<int> &level_offset, int lc, vector<int> &visited, double* vect, int k){
    priority_queue<pair<double,int>> candidates;
    vector<pair<double,int>> temp;
    while(!topk.empty()){
        temp.push_back(topk.top());
        candidates.push(topk.top());
        topk.pop();
    }
    for(auto x:temp){
        topk.push(x);
    }

    int ep, start, end;
    while(!candidates.empty()){
        ep = candidates.top().second;
        candidates.pop();
        start = indptr[ep] + level_offset[lc];
        end = indptr[ep] + level_offset[lc + 1];
        for(int i=start;i<end;i++){
            int node = index[i];
            if( node==-1 or visited[node]==1){
                continue;
            }
            visited[node]=1;
            double dist = cosine_dist(q,vect+ node*D);
            if(dist > topk.top().first and topk.size()>=k){
                continue;
            }
            topk.push({dist,node});
            while(topk.size()>k) topk.pop();
            candidates.push({dist,node});
        }
    }
}

void QueryHNSW(vector<double> &q,priority_queue<pair<double,int>> &topk, int ep, vector<int> &indptr, vector<int> &index, 
                    vector<int> &level_offset, int max_level, double* vect, int k){
    topk.push({cosine_dist(q,vect+ep*D),ep});
    vector<int> visited(L,0);
    visited[ep] = 1;
    for(int level = max_level;level>=0;level--){
        SearchLayer(q, topk, indptr, index, level_offset, level, visited, vect, k);
    }
}

void helper(int rank, int size, vector<vector<double>> &user, int ep, vector<int> &indptr, vector<int> &index, 
            vector<int> &level_offset, int max_level, double* vect, int k, vector<vector<int>> &results){

    int step = U/size;
    int start = rank*step;
    int end = (rank + 1)*step;
    if(rank == size-1){
        end = U;
    }

    //cout<<"Starting for rank "<<rank<<"/"<<size-1<<endl;
    int num_threads = omp_get_num_threads();
    //cout<<omp_get_num_threads()<<endl;


    int step__ = (end-start)/num_threads;
    for(int j=0;j<num_threads;j++){
        //cout<<"Start Thread "<<j+1<<"/"<<num_threads<<endl;
        int start__ = start + j*step__;
        int end__ = start + (j+1)*step__;
        if(j==num_threads-1) end__ = end;
    
        #pragma omp task shared(results)
        for(int i=start__;i<end__;i++)
        {
            //cout<<"Start for user "<<i<<endl;
            priority_queue<pair<double,int>> topk;
            QueryHNSW(user[i-start], topk, ep, indptr, index, level_offset, max_level, vect, k);
            vector<int> temp;
            while(!topk.empty()){
                temp.push_back(topk.top().second);
                topk.pop();
            }
            results[i] = temp;

            //cout<<"User "<<i<<" "<<temp[1]<<" "<<temp[0]<<endl;
        }
        //cout<<"Done Thread "<<j+1<<"/"<<num_threads<<endl;
    }
    #pragma omp taskwait
}

int main(int argc, char* argv[]){

    ifstream fs;

    int k = atoi(argv[1]);
    int C = atoi(argv[2]);

    int max_level;
    fs.open("bin_inputs/max_level", ios::binary | ios::in);
    fs.read((char *)&max_level, sizeof(max_level));
    fs.close();
    //cout << max_level <<'\n';

    int ep;
    fs.open("bin_inputs/ep", ios::binary | ios::in);
    fs.read((char *)&ep, sizeof(ep));
    fs.close();
    //cout << ep <<'\n';    

    vector<int> level;
    fs.open("bin_inputs/level",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    int file_size = fs.tellg()/4; 
    //cout << file_size <<'\n';
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
    //cout << level.size()<< "  "<<level[file_size-1] <<'\n';

    vector<int> index;
    fs.open("bin_inputs/index",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    file_size = fs.tellg()/4; 
    //cout << file_size <<'\n';
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
    //cout << index.size()<< "  "<<index[file_size-1] <<'\n';


    vector<int> indptr;
    fs.open("bin_inputs/indptr",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    file_size = fs.tellg()/4; 
    //cout << file_size <<'\n';
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
    //cout << indptr.size()<< "  "<<indptr[file_size-1] <<'\n';

    vector<int> level_offset;
    fs.open("bin_inputs/level_offset",ios::binary | ios::in); 
    fs.seekg(0, ios::end);
    file_size = fs.tellg()/4; 
    //cout << file_size <<'\n';
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
    //cout << level_offset.size()<< "  "<<level_offset[file_size-1] <<'\n';

    int rank, size;
    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    fs.open("bin_inputs/arguments",ios::binary | ios::in); 
    
    fs.read((char*)&U, sizeof U);
    fs.read((char*)&L, sizeof L);
    fs.read((char*)&D, sizeof D);
    fs.close();

    //cout<<L<<" "<<D<<" "<<U<<endl;


    int step = U/size;
    int start = rank*step;
    int end = (rank + 1)*step;
    if(rank == size-1){
        end = U;
    }
    vector<vector<double>> user(end-start,vector<double>());
    fs.open("bin_inputs/user",ios::binary | ios::in);
    fs.seekg(step*rank*D*sizeof(double), std::ios::beg);
    if (fs.is_open()){
        double tp;
        for(int j=start;j<end;j++){
            for(int k=0;k<D;k++){
                fs.read((char *)&tp, sizeof(tp));
                user[j-start].push_back(tp);
            }
        }

        fs.close(); 
    }
    //cout<<"At rank "<<rank<<": "<<user.size()<<" "<<user[4999].size()<<endl;

    // step = L/size;
    // start = rank*step;
    // end = (rank + 1)*step;
    // if(rank == size-1){
    //     end = L;
    // }


    // double* vect = new double[L*D];
    // double* temp = new double[(end-start)*D];
    // int* recvcount = new int[size];
    // int* disp = new int[size];
    // fs.open("bin_inputs/vect",ios::binary | ios::in);
    // fs.seekg(step*rank*D*sizeof(double), std::ios::beg);
    // if (fs.is_open()){
    //     double tp;
    //     for(int j=start;j<end;j++){
    //         for(int k=0;k<D;k++){
    //             fs.read((char *)&tp, sizeof(tp));
    //             temp[(j-start)*D + k] = tp;
    //         }
    //     }
    //     fs.close(); 
    // }

    // for(int i=0;i<size;i++){
    //     if(i==size-1){
    //         recvcount[i] = L*D - i*step*D;
    //     }else {
    //         recvcount[i] = step*D;
    //         disp[i] = i*step*D;
    //     }
    // }


    // MPI_Allgatherv(temp, (end-start)*D, MPI_DOUBLE , vect, recvcount, disp, MPI_DOUBLE, MPI_COMM_WORLD);


    double* vect = new double[L*D];
    fs.open("bin_inputs/vect",ios::binary | ios::in);
    if (fs.is_open()){
        double tp;
        for(int j=0;j<L;j++){
            for(int k=0;k<D;k++){
                fs.read((char *)&tp, sizeof(tp));
                vect[j*D + k] = tp;
            }
        }
        fs.close(); 
    }

    //cout<<"Reading done !"<<endl;


    step = U/size;
    start = rank*step;
    end = (rank + 1)*step;
    if(rank == size-1){
        end = U;
    }
    
    // #pragma omp parallel num_threads(C)
    // {
    //     #pragma omp single

    
    vector<vector<int>> results(U);

    

    #pragma omp parallel num_threads(C)
    {
        #pragma omp single
        {
            helper(rank, size, user, ep, indptr, index, level_offset, max_level, vect, k, results);
            //cout<<"DONE-------------------------------------------------------------"<<endl;
        }
    }
    //cout<<"Hello I am rank "<<rank<<endl;

    //cout<<"Rank "<<rank<<": "<<results.size()<<endl;

    //cout<<"I WAS HERE"<<endl;

    //cout<<results.size()<<" "<<results[0].size()<<endl;

    

    fstream wf("user_prediction", ios::out | ios::binary);
    wf.seekp(rank*step*k*sizeof(int), std::ios::beg);
    for(int i=start;i<end;i++){
        for(int j=k-1;j>=0;j--){
            //cout<<rank<<" "<<i<<" "<<j<<endl;
            wf.write((char*)&results[i][j], sizeof(int));
        }
        //if(rank==1) cout<<"User "<<i<<": "<<results[i-start][1]<<" "<<results[i-start][0]<<endl;
    }
    wf.close();
    

    //cout<<"I CAME HERE"<<endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0){

        fstream fsout;
        fsout.open("user_prediction.txt",ios::out);
        std::fstream ifs("user_prediction", std::ios::in | std::ios::binary);
        int tp;
        for(int i=0;i<U;i++){
            for(int j=0;j<k;j++){
                ifs.read((char *)&tp, sizeof(tp));
                fsout<<tp<<" ";
            }
            fsout<<endl;
        }
        fsout.close();
        ifs.close();
    }


    MPI_Finalize();
}
