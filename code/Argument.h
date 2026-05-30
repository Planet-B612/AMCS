#pragma once
#include <string>
#include <iostream>
#include "CommonStruc.h"
#include "graph.h"
using namespace std;

// global variables
Graph R_graph, O_graph;  
vector<double> Inv_inDeg;
int __numV_left;
double __eta_left;
int root_num;
double residual=0.0, decimal=1.0;
vector<bool> __Activated;
vector<int> seed_set;
vector<float> cost;
int round_num=0;
bool use_UB = true;
int num_gen = 0;
int num_addback = 0;


class Argument{
public:
    float eta_0 = 10000;
    string model = "IC";
    bool Rnd_cost = false;
    float eps = 0.7;
    double delta=0.01;
    vector<string> dataset = {"dblp"};
    // vector<int> data={2};
    int dataset_No = 0;
    string graph_path="../dataset/";
    string pw_path="../dataset/";
    int run_times=1;
    int times=0;
    int batch=2;
    bool seed_out=true;
    bool gene_ini_pw=false;
    bool real_time_pw=false;
    int numV=1;
    float left_num = 100.0;
    float over_pnodes = 10.0;
    
    Argument()
    {
    }
    void arg_update(int argn, char** argv)
    {
        for (int i = 0; i < argn; i++)
        {
            if (argv[i] == string("-dataset_No"))
                dataset_No = stoi(argv[i + 1]);
            if (argv[i] == string("-eta_0"))
                eta_0 = stof(argv[i + 1]);
            if (argv[i] == string("-eps"))
                eps = stof(argv[i + 1]);
            if (argv[i] == string("-batch"))
                batch = stof(argv[i + 1]);
            if (argv[i] == string("-run_times"))
                run_times = stoi(argv[i + 1]);
            if (argv[i] == string("-times")) 
                times = stoi(argv[i + 1]);
            if (argv[i] == string("-gene_ini_pw"))
            {
                if(argv[i + 1] == string("true"))
                    gene_ini_pw = true;
            }
            if (argv[i] == string("-real_time_pw"))
            {
                if(argv[i + 1] == string("true"))
                    real_time_pw = true;
            }
        }      
    }   
    void Initialization()
    {
        __Activated.clear();
        seed_set.clear();
        cost.clear();
        if(batch<1)
        {
            cout << "The batch size should be larger than 0, please check the input." << endl;
            exit(1);
        }
        string parameter_str = "# The parameters are: Rnd_cost=" + to_string(Rnd_cost) + ", model=" + model+", eps_Inf="+to_string(eps);
        cout << parameter_str << endl;
    }
    
    void load_cost_graph(int k)
    {
        string dataset_dir = graph_path + dataset[k];
        GraphBase::load_graph_directly_nbr_sorted(dataset_dir, O_graph, R_graph);
        numV=static_cast<int>(R_graph.size());
        __numV_left=numV;
        round_num=0;

        Inv_inDeg.resize(numV);
        for(int i=0;i<numV;i++)
        {
            Inv_inDeg[i]=static_cast<float>(1.0/R_graph[i].size());
        }

        cost.resize(numV);
        fill(cost.begin(), cost.end(), 1.0);
        __Activated.resize(numV);
        fill(__Activated.begin(), __Activated.end(), false);
        string cost_file;
        if (Rnd_cost)
        {
            cost_file = graph_path + dataset[k] + "_cost_Rand.txt";
            cout << "Using random cost file:  " << cost_file << endl;
        }
        else
        {
            cost_file = graph_path + dataset[k] + "_cost_001DEG.txt";
            cout << "Using degree based cost file: " << cost_file << endl;
        }
        std::ifstream inFile;
        inFile.open(cost_file);
        if (!inFile)
        {
            cout << "cannot open the cost file at " << cost_file << endl;
            exit(1);
        }
        inFile.seekg(0, std::ios_base::beg);
        for (int i = 0; i < numV; i++)
        {
            inFile >> cost[i];
        }
        inFile.close();
        if(gene_ini_pw)
        {
            for (int i = 0; i < run_times; i++){
                vector<vector<int>> PO(numV);
                ofstream out_pw;
                if(model=="IC")
                {   
                    out_pw.open(pw_path + dataset[dataset_No] + "_pw_ic" + to_string(i) + ".txt");
                    assert((!out_pw.fail()));
                    for(int i=0;i<(numV);i++)
                    {
                        auto nbrs=(O_graph)[i];
                        for(auto edge:nbrs)
                        {
                            if((dsfmt_gv_genrand_open_close())<Inv_inDeg[edge])
                            {
                                PO[i].push_back(edge);
                            }
                        }
                    }
                }
                else
                {
                    out_pw.open(pw_path + dataset[dataset_No] + "_pw_lt" + to_string(i) + ".txt");
                    assert((!out_pw.fail()));
                    for(int v=0;v<(numV);v++)
                    {
                        auto nbrs=(R_graph)[v];
                        auto nbrs_size=nbrs.size();
                        if(nbrs_size==0)
                        {
                            continue;
                        }
                        int index = dsfmt_gv_genrand_uint32_range(nbrs_size);
                        int u = nbrs[index];
                        PO[u].push_back(v);
                    }
                }
                for(int i=0;i<numV;i++)
                {
                    for(long unsigned int j=0;j<PO[i].size();j++)
                    {
                        out_pw<<i<<" "<<PO[i][j]<<endl;
                    }
                }
                out_pw.close();
            }
            cout << "generate PO Done" <<endl;
            exit(0);
        }
        return;
    }
    
    void seed_record(vector<int> seeds,int k,int i)
    {
        ofstream out_seeds("../results/seed/MINE_" + dataset[k] + "_" + to_string(eta_0*numV) + "_" + to_string(batch) + "_" + to_string(i) + ".txt", ios::out);
        assert((!out_seeds.fail()));
        for (auto node : seeds)
        {
            out_seeds << node << endl;
        }
        out_seeds.close();
    }

    bool graph_sort_check(const Graph &graph)
    {
        for (ulint i = 0; i < graph.size(); i++)
        {
            auto node = graph[i];
            if(node.size() <2)
            {
                continue;
            }
            for (ulint j = 0; j < node.size() - 1; j++)
            {
                if (node[j] > node[j + 1])
                {
                    cout << "Error: the graph is not sorted!" << endl;
                    return false;
                }
            }
        }
        return true;
    }
};