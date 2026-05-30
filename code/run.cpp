#include "../dSFMT/dSFMT.h"
#include "graph.h"
#include <iostream>
#include <vector>
#include "CommonStruc.h"
#include <cstring>
#include "Timer.h"
#include "Memory.h"
#include "MemoryUsage.h"
#include "Algorithm.h"
// #include "Argument.h"
#include <queue>
using namespace std;


int main(int argn, char **argv)
{
    
    Argument arg;  // claimed in Argument.h
    R_graph.clear(), O_graph.clear();  // global variables
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));  // the type must be uint32_t, to be accord with the function definition
    arg.arg_update(argn, argv);
    double avg_cost = 0.0;
    double avg_time = 0.0;
    while(arg.times < arg.run_times)
    {
        arg.Initialization();
        auto k = arg.dataset_No;
        arg.load_cost_graph(k);
        double original_eta=0;
        if (arg.eta_0 > 1) 
        {
            original_eta = arg.eta_0;
            arg.eta_0 = arg.eta_0/(1.0*arg.numV);
        }
        else
        {
            original_eta=arg.eta_0*arg.numV;
        }
        double total_cost=0.0;
        __eta_left = (arg.eta_0 ) * arg.numV;
        root_num=ceil(1.0/arg.eta_0);  // initial root number
        cout << " Running EMASS at eta = " << __eta_left << ", dataset = " << arg.dataset[k] << ", # node = " << arg.numV << ", eta = " << __eta_left <<", batch = "<<arg.batch<<", eps = "<<arg.eps<<", model = "<<arg.model <<", Rnd_cost = "<<arg.Rnd_cost<< ", real_time_pw = "<<arg.real_time_pw << endl;
        TAlg Alg(arg);

        vector<int> seeds;
        auto RR_info = Alg.AdaptiveSelect();
        seeds = seed_set;
        for (auto node : seeds)
        {
            total_cost += cost[node];
        }
        avg_cost += total_cost;
        avg_time += get<4>(RR_info);
        string results;
        results = "(" + arg.dataset[k] + ", eta = " + to_string(original_eta) + ", cost = " + to_string(total_cost) + ", % = " + to_string(1.0) + ", time = " + to_string(get<4>(RR_info)) +")" ;
        // result_bk << results << endl;
        cout << results << endl;
        results = arg.dataset[k] + " " + to_string(original_eta) + " EMASS" + " " + to_string(total_cost) + " " + to_string(1.0) + " " + to_string(get<4>(RR_info));
        cout << results << endl;
        Alg.release_memory();
        arg.times++;
    }
    cout << "Average cost: " << avg_cost / arg.run_times << endl;
    cout << "Average time: " << avg_time / arg.run_times << endl;

    return 0;
}