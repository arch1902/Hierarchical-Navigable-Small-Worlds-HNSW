#include <string>
#include <bits/stdc++.h>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]){

    fstream fs;
    ofstream fo;

    // cout<<k<<" recommendations using "<<C<<" threads !"<<endl;

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
    fo.open("bin_inputs/max_level",std::ios::binary);
    fo.write((char*)&max_level, sizeof max_level);
    fo.close();
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
    fo.open("bin_inputs/ep",std::ios::binary);
    fo.write((char*)&ep, sizeof ep);
    fo.close();
    cout<<"2"<<endl;

    fs.open("to_students/level.txt",ios::in);
    fo.open("bin_inputs/level",std::ios::binary); 
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

    fs.open("to_students/index.txt",ios::in); 
    fo.open("bin_inputs/index",std::ios::binary); 
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
    cout<<"4"<<endl;

    fs.open("to_students/indptr.txt",ios::in); 
    fo.open("bin_inputs/indptr",std::ios::binary); 
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

    fs.open("to_students/level_offset.txt",ios::in); 
    fo.open("bin_inputs/level_offset",std::ios::binary); 
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

    int L, D, U;
    bool flag = false;
    fs.open("to_students/vect.txt",ios::in); 
    fo.open("bin_inputs/vect",std::ios::binary);
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

    fs.open("to_students/user.txt",ios::in); 
    fo.open("bin_inputs/user",std::ios::binary);
    if (fs.is_open()){ 
        string tp;
        double td;
        while(getline(fs, tp)){ 
            if(tp=="") break;
            istringstream ss(tp);
    		string word;
    		while (ss >> word) 
    		{
                td = stod(word);
                fo.write((char*)&td, sizeof td);
    		}
            U +=1;
        }
        fs.close(); 
        fo.close();
    }
    cout<<"8"<<endl;
    fo.open("bin_inputs/arguments",std::ios::binary);
    fo.write((char*)&U, sizeof U);
    fo.write((char*)&L, sizeof L);
    fo.write((char*)&D, sizeof D);
    fo.close();
}