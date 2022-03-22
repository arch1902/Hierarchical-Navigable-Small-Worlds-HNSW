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

double cosine_dist(vector<double> &U, vector<double> &V)
{
    double dotp = 0.0, norm_u = 0.0, norm_v = 0.0 ;
    for(int i = 0; i < U.size(); i++) {
        dotp += U[i] * V[i] ;
        norm_u += U[i] * V[i] ;
        norm_v += U[i] * V[i] ;
    }
    return 1 - dotp/(sqrt(norm_u*norm_v)) ;
}

priority_queue<pair<double,int>> SearchLayer(vector<double> &q,priority_queue<pair<double,int>> candidates, vector<int> &indptr, vector<int> &index, 
                                                                vector<int> &level_offset, int lc, vector<int> &visited, vector<vector<double>> &vect, int k){
    priority_queue<pair<double,int>> topk = candidates;
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
            double dist = cosine_dist(q,vect[node]);
            if(dist > topk.top().first and topk.size()==k){
                continue;
            }
            topk.push({dist,node});
            if(topk.size()>k) topk.pop();
            candidates.push({dist,node});
        }
    }
    return topk;
}

void QueryHNSW(vector<double> &q,priority_queue<pair<double,int>> &topk, int ep, vector<int> &indptr, vector<int> &index, 
                    vector<int> &level_offset, int max_level, vector<vector<double>> &vect, int k){
    topk.push({cosine_dist(q,vect[ep]),ep});
    vector<int> visited(vect.size(),0);
    visited[ep] = 1;
    for(int level = max_level;level>=0;level--){
        topk = SearchLayer(q, topk, indptr, index, level_offset, level, visited, vect, k);
    }
}

void helper(int rank, int size, vector<vector<double>> &user, int ep, vector<int> &indptr, vector<int> &index, 
            vector<int> &level_offset, int max_level, vector<vector<double>> &vect, int k, vector<vector<int>> &results){

    int step = user.size()/size;
    int start = rank*step;
    int end = (rank + 1)*step;
    if(rank == size-1){
        end = user.size();
    }

    cout<<"Starting for rank "<<rank<<"/"<<size-1<<endl;
    for(int i=start;i<end;i++){

        #pragma omp task
        {
            cout<<"Start for user "<<i<<endl;
            priority_queue<pair<double,int>> topk;
            QueryHNSW(user[i], topk, ep, indptr, index, level_offset, max_level, vect, k);
            vector<int> temp;
            while(!topk.empty()){
                temp.push_back(topk.top().second);
                topk.pop();
            }
            results[i] = temp;
            cout<<"Done for user "<<i<<endl;
        }
    }
    #pragma omp taskwait
}

int main(int argc, char* argv[]){

    fstream fs;

    int k = atoi(argv[1]);
    int C = atoi(argv[2]);

    cout<<k<<" recommendations using "<<C<<" threads !"<<endl;

    int max_level;
    fs.open("to_students/max_level.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            max_level = stoi(tp);
            break;
        }
        fs.close(); 
    }
    cout<<"1"<<endl;

    int ep;
    fs.open("to_students/ep.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            ep = stoi(tp);
            break;
        }
        fs.close(); 
    }
    cout<<"2"<<endl;

    vector<int> level;
    fs.open("to_students/level.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            level.push_back(stoi(tp));
        }
        fs.close(); 
    }
    cout<<"3"<<endl;

    vector<int> index;
    fs.open("to_students/index.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            index.push_back(stoi(tp));
        }
        fs.close(); 
    }
    cout<<"4"<<endl;

    vector<int> indptr;
    fs.open("to_students/indptr.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            indptr.push_back(stoi(tp));
        }
        fs.close(); 
    }
    cout<<"5"<<endl;

    vector<int> level_offset;
    fs.open("to_students/level_offset.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            level_offset.push_back(stoi(tp));
        }
        fs.close(); 
    }
    cout<<"6"<<endl;

    vector<vector<double>> vect;
    fs.open("to_students/vect.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;

        while(getline(fs, tp)){ 
            if(tp=="") break;
            istringstream ss(tp);
    		string word;
            vector<double> temp;
    		while (ss >> word) 
    		{
        		temp.push_back(stod(word));
    		}
            vect.push_back(temp);
        }
        fs.close(); 
    }
    cout<<"7"<<endl;

    vector<vector<double>> user;
    fs.open("to_students/user.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            istringstream ss(tp);
    		string word;
            vector<double> temp;
    		while (ss >> word) 
    		{
        		temp.push_back(stod(word));
    		}
            user.push_back(temp);
        }
        fs.close(); 
    }
    cout<<"8"<<endl;

    cout<<"---------------------------- Data read"<<endl;


    // cout<<level.size()<<endl;
    // cout<<index.size()<<endl;
    // cout<<indptr.size()<<endl;
    // cout<<level_offset.size()<<endl;
    // cout<<vect.size()<<endl;
    // cout<<user.size()<<endl;


    int rank, size;
    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<vector<int>> results(user.size());

    #pragma omp parallel num_threads(C)
    {
        #pragma omp single
        helper(rank, size, user, ep, indptr, index, level_offset, max_level, vect, k, results);
    }

    MPI_Finalize();

    fstream fsout;
	fsout.open("out.txt",ios::out);

    for(int i=0;i<user.size();i++){
        for(int j=0;j<k;j++){
            fsout<<results[i][j]<<" ";
        }
        fsout<<endl;
    }
    fsout.close();
}