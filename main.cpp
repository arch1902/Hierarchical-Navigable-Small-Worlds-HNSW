#include <string>
#include <mpi.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

/*

Path of the to_students folder from personal scratch
../../../../../scratch/cse/phd/anz198717/TA/COL380/A3/to_students

*/

// void f(int tid, ){

// }

int main(int argc, char* argv[]){

    fstream fs;

    int max_level;
    fs.open("to_students/max_level.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            max_level = stoi(tp);
            break;
        }
        fs.close(); //close the file object.
    }
    //cout<<"1"<<endl;

    int ep;
    fs.open("to_students/ep.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            ep = stoi(tp);
            break;
        }
        fs.close(); //close the file object.
    }
    //cout<<"2"<<endl;

    vector<int> level;
    fs.open("to_students/level.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            level.push_back(stoi(tp));
        }
        fs.close(); //close the file object.
    }
    //cout<<"3"<<endl;

    vector<int> index;
    fs.open("to_students/index.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            index.push_back(stoi(tp));
        }
        fs.close(); //close the file object.
    }
    //cout<<"4"<<endl;

    vector<int> indptr;
    fs.open("to_students/indptr.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            indptr.push_back(stoi(tp));
        }
        fs.close(); //close the file object.
    }
    //cout<<"5"<<endl;

    vector<int> level_offset;
    fs.open("to_students/level_offset.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            level_offset.push_back(stoi(tp));
        }
        fs.close(); //close the file object.
    }
    //cout<<"6"<<endl;

    vector<vector<int>> vect;
    fs.open("to_students/vect.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;

        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            istringstream ss(tp);
    		string word;
            vector<int> temp;
    		while (ss >> word) 
    		{
        		temp.push_back(stoi(word));
    		}
            vect.push_back(temp);
        }
        fs.close(); //close the file object.
    }
    //cout<<"7"<<endl;

    vector<vector<int>> user;
    fs.open("to_students/user.txt",ios::in); //open a file to perform read operation using file object
    if (fs.is_open()){ //checking whether the file is open
        string tp;
        while(getline(fs, tp)){ //read data from file object and put it into string.
            if(tp=="") break;
            istringstream ss(tp);
    		string word;
            vector<int> temp;
    		while (ss >> word) 
    		{
        		temp.push_back(stoi(word));
    		}
            user.push_back(temp);
        }
        fs.close(); //close the file object.
    }
    //cout<<"8"<<endl;


    int rank, size;
    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // f();

    MPI_Finalize();
 
}
