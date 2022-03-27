#include <string>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]){

    string in_path = argv[1];
    string out_path = argv[2];
    fstream fs;
    ofstream fo;

    // cout<<k<<" recommendations using "<<C<<" threads !"<<endl;

    int max_level;
    fs.open(in_path+"/max_level.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            max_level = stoi(tp);
            break;
        }
        fs.close(); 
    }
    fo.open(out_path+"/max_level",std::ios::binary);
    fo.write((char*)&max_level, sizeof max_level);
    fo.close();
    cout<<"1"<<endl;

    int ep;
    fs.open(in_path+"/ep.txt",ios::in); 
    if (fs.is_open()){ 
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            ep = stoi(tp);
            break;
        }
        fs.close(); 
    }
    fo.open(out_path+"/ep",std::ios::binary);
    fo.write((char*)&ep, sizeof ep);
    fo.close();
    cout<<"2"<<endl;

    fs.open(in_path+"/level.txt",ios::in);
    fo.open(out_path+"/level",std::ios::binary); 
    if (fs.is_open()){ 
        int ti;
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            ti = stoi(tp);
            fo.write((char*)&ti, sizeof ti);
        }
        fo.close();
        fs.close(); 
    }
    cout<<"3"<<endl;

    fs.open(in_path+"/index.txt",ios::in); 
    fo.open(out_path+"/index",std::ios::binary); 
    int ind_size = 0;
    if (fs.is_open()){ 
        int ti;
        string tp;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            ind_size ++;
            ti = stoi(tp);
            fo.write((char*)&ti, sizeof ti);
        }
        fo.close();
        fs.close(); 
    }
    cout<<"4"<<endl;

    fs.open(in_path+"/indptr.txt",ios::in); 
    fo.open(out_path+"/indptr",std::ios::binary); 
    if (fs.is_open()){ 
        string tp;
        int ti;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            ti = stoi(tp);
            fo.write((char*)&ti, sizeof ti);
        }
        fo.close();
        fs.close(); 
    }
    cout<<"5"<<endl;

    fs.open(in_path+"/level_offset.txt",ios::in); 
    fo.open(out_path+"/level_offset",std::ios::binary); 
    if (fs.is_open()){ 
        string tp;
        int ti;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            ti = stoi(tp);
            fo.write((char*)&ti, sizeof ti);
        }
        fs.close(); 
        fo.close();
    }
    cout<<"6"<<endl;

    int L=0, D=0, U=0;
    bool flag = true;
    fs.open(in_path+"/vect.txt",ios::in); 
    fo.open(out_path+"/vect",std::ios::binary);
    if (fs.is_open()){ 
        string tp;

        while(getline(fs, tp)){ 
            if(tp=="") break;
            istringstream ss(tp);
    		string word;
            double td;
    		while (ss >> word) 
    		{
                if(flag){D += 1;}
                td = stod(word);
                fo.write((char*)&td, sizeof td);
    		}
            L +=1;
            flag = false;
        }
        fo.close();
        fs.close(); 
    }
    cout<<"7"<<endl;

    // fs.open(in_path+"/user.txt",ios::in); 
    // fo.open(out_path+"/user",std::ios::binary);
    // if (fs.is_open()){ 
    //     string tp;
    //     double td;
    //     while(getline(fs, tp)){ 
    //         if(tp=="") break;
    //         istringstream ss(tp);
    // 		string word;
    // 		while (ss >> word) 
    // 		{
    //             td = stod(word);
    //             fo.write((char*)&td, sizeof td);
    // 		}
    //         U +=1;
    //     }
    //     fs.close(); 
    //     fo.close();
    // }
    // cout<<"8"<<endl;
    fo.open(out_path+"/arguments",std::ios::binary);
    fo.write((char*)&U, sizeof U);
    fo.write((char*)&L, sizeof L);
    fo.write((char*)&D, sizeof D);
    fo.write((char*)&ind_size, sizeof ind_size);
    fo.close();
}