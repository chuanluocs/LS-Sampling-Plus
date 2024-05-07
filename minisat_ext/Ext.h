#ifndef ExtMinisat_h
#define ExtMinisat_h

#include <vector>
#include <numeric>
#include <random>
#include <utility>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>

using std::vector;
using std::pair;
using std::string;
using std::to_string;
using std::ifstream;
using std::ofstream;

namespace ExtMinisat {

class SamplingSolver{
public:
    SamplingSolver(int nvar, const vector<vector<int> >& clauses, int seed, bool do_shuffle, int enable_VSIDS);
    virtual ~SamplingSolver();
    void set_prob(const vector<pair<int, int> >& sample_prob);
    void set_prob2(string cnf_path, const vector<pair<int, int> >& sample_prob, int seed, int size);
    void get_solution(vector<int>& tc);
    void get_solution2(string cnf_path, vector<vector<int> >& tc, int seed, int size);
    vector<int> get_solution();
protected:
    void* internal_solver;
    std::mt19937_64 gen;
    bool enable_shuffling;
    vector<int> vecc;
};

}

#endif