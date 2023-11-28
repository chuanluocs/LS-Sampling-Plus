#include "slstestcasesampler.h"
#include <unistd.h>
#include <bits/stdc++.h>
#include <pthread.h>

using namespace std;

const int kMaxPbOCCSATSeed = 10000000;

SLSTestcaseSampler::SLSTestcaseSampler(string cnf_file_path, int seed): rng_2(seed)
{
    cnf_file_path_ = cnf_file_path;
    seed_ = seed;
    rng_.seed(seed_);
    SetDefaultPara();
}

SLSTestcaseSampler::~SLSTestcaseSampler()
{
    delete pbo_solver_;
}

void SLSTestcaseSampler::SetDefaultPara()
{
    testcase_set_size_ = 100;
    candidate_set_size_ = 100;

    testcase_set_.resize(testcase_set_size_);
    candidate_sample_init_solution_set_.resize(candidate_set_size_);
    candidate_testcase_set_.resize(candidate_set_size_);

    p_init_tuple_info_ = &SLSTestcaseSampler::Init2TupleInfo;
    p_update_tuple_info_ = &SLSTestcaseSampler::Update2TupleInfo;
    p_calculate_tuple_coverage_ = &SLSTestcaseSampler::Cal2TupleCoverage;

    flag_use_weighted_sampling_ = true;
    p_init_sample_weight_ = &SLSTestcaseSampler::InitSampleWeightByAppearance;
    p_update_sample_weight_ = &SLSTestcaseSampler::UpdateSampleWeightByAppearance;

    flag_use_context_aware_ = true;
    p_init_context_aware_flip_priority_ = &SLSTestcaseSampler::InitContextAwareFlipPriority;
    p_update_context_aware_flip_priority_ = &SLSTestcaseSampler::UpdateContextAwareFlipPriorityBySampleWeight;
    flag_use_cnf_reduction_ = true;
    p_reduce_cnf_ = &SLSTestcaseSampler::ReduceCNF;
    int pos = cnf_file_path_.find_last_of('/');
    cnf_instance_name_ = cnf_file_path_.substr(pos + 1);
    cnf_instance_name_.replace(cnf_instance_name_.find(".cnf"), 4, "");
    flag_reduced_cnf_as_temp_ = true;
    reduced_cnf_file_path_ = "/tmp/" + get_random_prefix() + "_reduced.cnf";
    testcase_set_save_path_ = "./" + cnf_instance_name_ + "_testcase_set.txt";

    p_get_gain = &SLSTestcaseSampler::get_gain_2_wise;

    flag_fix_t_wise_optimize_ = false;
    sample_cnt_ = 1000000;

    parallel_num = 1;
    use_cdcl = true;
    stage_cha = false;
}

void SLSTestcaseSampler::SetTestcaseSetSize(int testcase_set_size)
{
    testcase_set_size_ = testcase_set_size;
    testcase_set_.resize(testcase_set_size_);
}

void SLSTestcaseSampler::SetCandidateSetSize(int candidate_set_size)
{
    candidate_set_size_ = candidate_set_size;
    candidate_sample_init_solution_set_.resize(candidate_set_size_);
    candidate_testcase_set_.resize(candidate_set_size_);
}

void SLSTestcaseSampler::SetWeightedSamplingMethod(bool use_weighted_sampling)
{
    flag_use_weighted_sampling_ = use_weighted_sampling;
    if (flag_use_weighted_sampling_)
    {
        p_init_sample_weight_ = &SLSTestcaseSampler::InitSampleWeightByAppearance;
        p_update_sample_weight_ = &SLSTestcaseSampler::UpdateSampleWeightByAppearance;
    }
    else
    {
        p_init_sample_weight_ = &SLSTestcaseSampler::InitSampleWeightUniformly;
        p_update_sample_weight_ = &SLSTestcaseSampler::EmptyFunRetVoid;
    }
}

void SLSTestcaseSampler::SetContextAwareMethod(bool use_context_aware)
{
    flag_use_context_aware_ = use_context_aware;
    if (flag_use_context_aware_)
    {
        p_init_context_aware_flip_priority_ = &SLSTestcaseSampler::InitContextAwareFlipPriority;
        p_update_context_aware_flip_priority_ = &SLSTestcaseSampler::UpdateContextAwareFlipPriorityBySampleWeight;
    }
    else
    {
        p_init_context_aware_flip_priority_ = &SLSTestcaseSampler::EmptyFunRetVoid;
        p_update_context_aware_flip_priority_ = &SLSTestcaseSampler::EmptyFunRetVoid;
    }
}

void SLSTestcaseSampler::SetCNFReductionMethod(bool use_cnf_reduction)
{
    flag_use_cnf_reduction_ = use_cnf_reduction;
    if (flag_use_cnf_reduction_)
    {
        p_reduce_cnf_ = &SLSTestcaseSampler::ReduceCNF;
    }
    else
    {
        p_reduce_cnf_ = &SLSTestcaseSampler::EmptyFunRetVoid;
    }
}

void SLSTestcaseSampler::GenerateInitTestcase()
{
    bool is_sat = pbo_solver_->solve();
    num_var_ = pbo_solver_->get_var_num();
    if (is_sat)
    {
        testcase_set_[0] = pbo_solver_->get_sat_solution();
        num_generated_testcase_ = 1;
    }
    else
    {
        cout << "c PbOCCSAT Failed to Find Initial Sat Solution!" << endl;
    }
}

long long SLSTestcaseSampler::Get1TupleMapIndex(long i, long v_i)
{
    long long index_comb = v_i;
    long long index = index_comb * num_combination_all_possible_ + i;
    return index;
}

long long SLSTestcaseSampler::Get2TupleMapIndex(long i, long v_i, long j, long v_j)
{
    if (i >= j)
    {
        cout << "c Wrong index order!" << endl;
        return -1;
    }
    else
    {
        long long index_comb = v_j + v_i * 2;
        long long index = index_comb * num_combination_all_possible_ + (2ll * num_var_ - i - 1) * i / 2 + j - i - 1;
        return index;
    }
}

long long SLSTestcaseSampler::Get3TupleMapIndex(long i, long v_i, long j, long v_j, long k, long v_k) //3wise_no_ASF
{
    if (i >= j || j >= k)
    {
        cout << "c Wrong index order!";
        return -1;
    }
    else
    {
        long long num_var = num_var_;
        long long index_comb = v_k + v_j * 2 + v_i * 4;
        long long index = index_comb * num_combination_all_possible_
        + (num_var * (i * num_var - (i + 2) * i) + i * (i + 1) * (i + 2) / 3) / 2 + (2 * num_var - i - j - 2) * (j - i - 1) / 2 + k - j - 1;
        return index;
    }
}

int SLSTestcaseSampler::GetHammingDistance(const vector<int> &vec_1, const vector<int> &vec_2)
{
    int hamming_dis = 0;
    for (int v = 0; v < num_var_; v++)
    {
        hamming_dis += abs(vec_1[v] - vec_2[v]);
    }
    return hamming_dis;
}

void SLSTestcaseSampler::Init1TupleInfo()
{
    num_combination_all_possible_ = (long long)num_var_ ;
    num_tuple_all_possible_ = num_combination_all_possible_ * 2;
    already_t_wise = SimpleBitSet(num_tuple_all_possible_);
    count_each_var_positive_uncovered_.resize(num_var_,  2);
    count_each_var_negative_uncovered_.resize(num_var_,  2);
    num_tuple_ = 0;
}

void SLSTestcaseSampler::Update1TupleInfo()
{
    int index_testcase = num_generated_testcase_ - 1;
    const vector<int>& testcase = testcase_set_[index_testcase];
    for (int i = 0; i < num_var_ ; i++)
    {
        long long index_tuple = Get1TupleMapIndex(i, testcase[i]);
        int res = already_t_wise.set(index_tuple);
        if (res)
        {
            num_tuple_++;
            if (testcase[i] == 1)
            {
                count_each_var_positive_uncovered_[i]--;
            }
            else
            {
                count_each_var_negative_uncovered_[i]--;
            }
        }
    }
}

void SLSTestcaseSampler::Init2TupleInfo()
{
    num_combination_all_possible_ = (long long)num_var_ * (num_var_ - 1) / 2;
    num_tuple_all_possible_ = num_combination_all_possible_ * 4;
    already_t_wise = SimpleBitSet(num_tuple_all_possible_);
    count_each_var_positive_uncovered_.resize(num_var_, (num_var_ - 1) * 2);
    count_each_var_negative_uncovered_.resize(num_var_, (num_var_ - 1) * 2);
    num_tuple_ = 0;
}

void SLSTestcaseSampler::Update2TupleInfo()
{
    int index_testcase = num_generated_testcase_ - 1;
    const vector<int>& testcase = testcase_set_[index_testcase];
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, testcase[i], j, testcase[j]);
            int res = already_t_wise.set(index_tuple);
            if (res)
            {
                num_tuple_++;
                if (testcase[i] == 1)
                {
                    count_each_var_positive_uncovered_[i]--;
                }
                else
                {
                    count_each_var_negative_uncovered_[i]--;
                }
                if (testcase[j] == 1)
                {
                    count_each_var_positive_uncovered_[j]--;
                }
                else
                {
                    count_each_var_negative_uncovered_[j]--;
                }
            }
        }
    }
}

void SLSTestcaseSampler::Init3TupleInfo()  //3wise_no_ASF
{
    long long num_var = num_var_;
    num_combination_all_possible_ = num_var * (num_var - 1) * (num_var - 2) / 6;
    num_tuple_all_possible_ = num_combination_all_possible_ * 8;
    already_t_wise = SimpleBitSet(num_tuple_all_possible_);
    count_each_var_positive_uncovered_.resize(num_var, (num_var - 1) * (num_var - 2) * 2);
    count_each_var_negative_uncovered_.resize(num_var, (num_var - 1) * (num_var - 2) * 2);
    num_tuple_ = 0;
}

void SLSTestcaseSampler::Update3TupleInfo()  //3wise_no_ASF
{
    int index_testcase = num_generated_testcase_ - 1;
    const vector<int>& testcase = testcase_set_[index_testcase];
    for (int i = 0; i < num_var_ - 2; i++)
    {
        for (int j = i + 1; j < num_var_ - 1; j++)
        {
            for (int k = j + 1; k < num_var_; k++)
            {
                long long index_tuple = Get3TupleMapIndex(i, testcase[i], j, testcase[j], k, testcase[k]);
                int res = already_t_wise.set(index_tuple);
                if (res)
                {
                    num_tuple_++;
                }
            }
        }
    }
}

void SLSTestcaseSampler::UpdateXTupleInfo()
{
    
    if (selected_candidate_bitset_.init_iter()){
        do {
            size_t tp_id = selected_candidate_bitset_.iter_get();
            already_t_wise.set(tp_id);
            ++num_tuple_;
        } while (selected_candidate_bitset_.iter_next());
    }

    size_t uncovered_cnt = sample_uncovered.size();
    for (size_t i = 0; i + 1 < uncovered_cnt; ){
        if (selected_candidate_bitset_.get(sample_uncovered[i])){
            swap(sample_uncovered[i], sample_uncovered[uncovered_cnt - 1]);
            sample_uncovered.pop_back();
            --uncovered_cnt;
        } else ++i;
    }
    if (uncovered_cnt > 0 &&
        selected_candidate_bitset_.get(sample_uncovered.back())){
        sample_uncovered.pop_back();
    }
}

void SLSTestcaseSampler::clear_pq()
{
    for (auto& bs: pq){
        bs.second.difference_(selected_candidate_bitset_);
    }

    sort(pq.begin(), pq.end(), [](const pair<vector<int>, SimpleBitSet>& si, const pair<vector<int>, SimpleBitSet>& sj){
        return si.second.count() > sj.second.count();
    });

    int cur_pqsize = (int) pq.size();
    while (cur_pqsize > candidate_set_size_ || (cur_pqsize > 0 && pq.back().second.count() == 0)){
        pq.pop_back();
        pq_idx.pop_back();
        --cur_pqsize;
    }
}

void SLSTestcaseSampler::InitSampleWeightByAppearance()
{
    var_positive_appearance_count_.resize(num_var_);
    var_positive_sample_weight_.resize(num_var_);
}

void SLSTestcaseSampler::UpdateSampleWeightByAppearance()
{
    int new_testcase_index = num_generated_testcase_ - 1;
    vector<int> new_testcase = testcase_set_[new_testcase_index];
    for (int v = 0; v < num_var_; v++)
    {
        var_positive_appearance_count_[v] += new_testcase[v];
        var_positive_sample_weight_[v] = 1. - double(var_positive_appearance_count_[v]) / num_generated_testcase_;
    }
}

void SLSTestcaseSampler::InitSampleWeightByUncoveredCount()
{
    var_positive_sample_weight_.resize(num_var_);
}

void SLSTestcaseSampler::UpdateSampleWeightByUncoveredCount()
{
    for (int v = 0; v < num_var_; v++)
    {
        var_positive_sample_weight_[v] = double(count_each_var_positive_uncovered_[v]) /
                                         (count_each_var_positive_uncovered_[v] + count_each_var_negative_uncovered_[v]);
    }
}

void SLSTestcaseSampler::InitSampleWeightUniformly()
{
    var_positive_sample_weight_.resize(num_var_, 0.5);
}

void SLSTestcaseSampler::InitContextAwareFlipPriority()
{
    context_aware_flip_priority_.resize(num_var_, 0.);
}

void SLSTestcaseSampler::UpdateContextAwareFlipPriority()
{
    vector<int> init_solution = candidate_sample_init_solution_set_[selected_candidate_index_];
    vector<int> sat_testcase = candidate_testcase_set_[selected_candidate_index_];
    vector<int> var_flip_count(num_var_);
    for (int v = 0; v < num_var_; v++)
    {
        var_flip_count[v] = abs(init_solution[v] - sat_testcase[v]);
    }
    for (int v = 0; v < num_var_; v++)
    {
        context_aware_flip_priority_[v] += var_flip_count[v];
    }
    pbo_solver_->set_var_flip_priority_ass_unaware(context_aware_flip_priority_);
}

void SLSTestcaseSampler::UpdateContextAwareFlipPriorityBySampleWeight(const vector<int> &init_solution)
{
    for (int v = 0; v < num_var_; v++)
    {
        if (init_solution[v])
            context_aware_flip_priority_[v] = var_positive_sample_weight_[v];
        else
            context_aware_flip_priority_[v] = 1. - var_positive_sample_weight_[v];
    }

    pbo_solver_->set_var_flip_priority_ass_aware(context_aware_flip_priority_);
}

void SLSTestcaseSampler::ReduceCNF()
{
    string cmd = "./bin/coprocessor -enabled_cp3 -up -subsimp -no-bve -no-bce"
                 " -no-dense -dimacs=" +
                 reduced_cnf_file_path_ + " " + cnf_file_path_;

    int return_val = system(cmd.c_str());

    cnf_file_path_ = reduced_cnf_file_path_;
}

void SLSTestcaseSampler::InitPbOCCSATSolver()
{
    pbo_solver_ = new PbOCCSATSolver(cnf_file_path_, rng_.next(kMaxPbOCCSATSeed));
}

void SLSTestcaseSampler::Init()
{
    (this->*p_reduce_cnf_)();
    InitPbOCCSATSolver();
    GenerateInitTestcase();
    (this->*p_init_tuple_info_)();
    (this->*p_init_sample_weight_)();
    (this->*p_init_context_aware_flip_priority_)();
}

vector<int> SLSTestcaseSampler::GetSatTestcaseWithGivenInitSolution(const vector<int> &init_solution)
{
    int pbo_seed = rng_.next(kMaxPbOCCSATSeed);
    pbo_solver_->set_seed(pbo_seed);

    pbo_solver_->set_init_solution(init_solution);

    (this->*p_update_context_aware_flip_priority_)(init_solution);

    bool is_sat = pbo_solver_->solve();

    if (is_sat)
    {
        return pbo_solver_->get_sat_solution();
    }
    else
    {
        cout << "c PbOCCSAT Failed to Find Sat Solution!" << endl;
    }
}

vector<int> SLSTestcaseSampler::GetWeightedSampleInitSolution()
{
    vector<int> weighted_sample_init_solution(num_var_, 0);
    for (int v = 0; v < num_var_; v++)
    {
        if (rng_.nextClosed() < var_positive_sample_weight_[v])
        {
            weighted_sample_init_solution[v] = 1;
        }
    }
    return weighted_sample_init_solution;
}

struct para {
    ExtMinisat::SamplingSolver *cdcl_sampler;
    vector<pair<int, int> > *sample_prob;
    std::vector<std::vector<int>> *candidate_testcase_set_;
    int start,end;
};

void *SLSTestcaseSampler::pthread_func3(void *arg){
    struct para *p = (struct para*)arg;
    int i = p->start;
    int j = p->end;

    for(; i<j; i++){
        p->cdcl_sampler->set_prob(*(p->sample_prob));
        p->cdcl_sampler->get_solution(p->candidate_testcase_set_->at(i));
    }
    free(p);
    pthread_exit(NULL);
}

void SLSTestcaseSampler::GenerateCandidateTestcaseSet()
{
    //
    vector<pair<int, int> > res;
    res.reserve(num_var_);
    for (int i = 0; i < num_var_; ++i){
        if (!flag_use_context_aware_){
            int v1 = 1;
            int v2 = 1;
            res.emplace_back(v1, v2);
            // cout << "not context aware." << endl;
        }
        else {
            int v1 = (num_generated_testcase_ - var_positive_appearance_count_[i]);
            int v2 = var_positive_appearance_count_[i];
            res.emplace_back(v1, v2);
        }
    }
    //

    for (int i = 0; i < candidate_set_size_; i++)
    {
        if(use_cdcl){
            cdcl_sampler->set_prob(res);
            cdcl_sampler->get_solution(candidate_testcase_set_[i]);
        }
        else{
            candidate_sample_init_solution_set_[i] = GetWeightedSampleInitSolution();
            candidate_testcase_set_[i] = GetSatTestcaseWithGivenInitSolution(candidate_sample_init_solution_set_[i]);
        }
    }
}

int SLSTestcaseSampler::SelectTestcaseFromCandidateSet()
{
    int max_hamming_distance = 0;
    int best_testcase_index = 0;
    for (int i = 0; i < candidate_set_size_; i++)
    {
        int hamming_distance = 0;
        for (int j = 0; j < num_generated_testcase_; j++)
        {
            hamming_distance += GetHammingDistance(candidate_testcase_set_[i], testcase_set_[j]);
        }
        if (hamming_distance > max_hamming_distance)
        {
            best_testcase_index = i;
            max_hamming_distance = hamming_distance;
        }
    }
    return best_testcase_index;
}

int SLSTestcaseSampler::SelectTestcaseFromCandidateSetByTupleNum()
{
    for (int i = 0; i < candidate_set_size_; ++i){
        pq.emplace_back(candidate_testcase_set_[i], (this->*p_get_gain)(candidate_testcase_set_[i]));
        pq_idx.push_back(0);
    }
    
    iota(pq_idx.begin(), pq_idx.end(), 0);
    sort(pq_idx.begin(), pq_idx.end(), [&](int i, int j){
        if (pq[i].second.count() != pq[j].second.count()) return pq[i].second.count() > pq[j].second.count();
        return i < j;
    });

    int res = pq_idx[0];
    return res;
}

struct para3 {
    std::vector<std::vector<int>> *candidate_testcase_set_;
    SLSTestcaseSampler *th;
    int paral;
    int start,end;
    bool sta;
    vector< vector<SimpleBitSet> > *ret;
};

void *SLSTestcaseSampler::pthread_func4(void *arg){
    struct para3 *p = (struct para3*)arg;
    int i = p->start;
    int j = p->end;
    bool s = p->sta;
    std::vector<std::vector<int>> *candidate_testcase_set_ = p->candidate_testcase_set_;
    SLSTestcaseSampler *th=p->th;
    vector<SimpleBitSet> ret1;
    SimpleBitSet res;
    for(; i<j; i++){
        if(s)
            res = th->get_gain_by_samples(candidate_testcase_set_->at(i));
        else
            res = th->get_gain_2_wise(candidate_testcase_set_->at(i));
        ret1.emplace_back(res);
    }
    p->ret->at(p->paral) = ret1;
    free(p);
    pthread_exit(NULL);
}

int SLSTestcaseSampler::SelectTestcaseFromCandidateSetByTupleNum2(){
    int div = candidate_set_size_ / parallel_num;
    pthread_t tid[parallel_num] = {0};
    vector< vector<SimpleBitSet> > ret_table(parallel_num);
    for (int i = 0; i < parallel_num; i++){
        int *p = (int*)malloc(sizeof(*p));
        *p = i;
        struct para3 *myp = (struct para3*)malloc(sizeof(struct para3));
        myp->candidate_testcase_set_ = &candidate_testcase_set_;
        myp->start = *p * div;
        myp->th = this;
        myp->sta = stage_cha;
        myp->paral = *p;
        myp->ret = &(ret_table);
        if( *p == parallel_num - 1 )myp->end = candidate_set_size_;
        else myp->end = (*p + 1) * div;
        pthread_create(&tid[*p],NULL,pthread_func4,(void*)myp);
        free(p);
    }
    for (int i = 0; i < parallel_num; i++){
        pthread_join(tid[i],NULL);
        int start = i*div;
        int k = start;
        int l;
        if( i == parallel_num - 1 )
            l = candidate_set_size_;
        else 
            l = ( i + 1 ) * div;
        for(; k<l; k++){
            pq.emplace_back(candidate_testcase_set_.at(k), ret_table.at(i).at(k-start));
            pq_idx.push_back(0);
        }
    }
    
    // sort the cache by the scores in descending order
    iota(pq_idx.begin(), pq_idx.end(), 0);
    sort(pq_idx.begin(), pq_idx.end(), [&](int i, int j){
        if (pq[i].second.count() != pq[j].second.count()) return pq[i].second.count() > pq[j].second.count();
        return i < j;
    });
    return pq_idx[0];
}

void SLSTestcaseSampler::GenerateTestcase()
{
    GenerateCandidateTestcaseSet();
    selected_candidate_index_ = SelectTestcaseFromCandidateSetByTupleNum();
    int testcase_index = num_generated_testcase_;
    testcase_set_[testcase_index] = pq[selected_candidate_index_].first;
    selected_candidate_bitset_ = pq[selected_candidate_index_].second;
}

void SLSTestcaseSampler::stage_change(int cur_phase){
    stage_cha = true;
    sample_dim_ = cur_phase;
    tuple_sampler(sample_dim_, sample_cnt_, tp_sample);

    already_t_wise = SimpleBitSet(sample_cnt_);
    num_tuple_ = 0;
    num_tuple_all_possible_ = sample_cnt_;

    p_update_tuple_info_ = &SLSTestcaseSampler::UpdateXTupleInfo;
    p_get_gain = &SLSTestcaseSampler::get_gain_by_samples;
    
    pq.clear();
    pq_idx.clear();

    sample_uncovered.clear();
    for (int i = 0; i < sample_cnt_; ++i){
        bool ok = false;
        for (int j = 0; j < num_generated_testcase_; ++j){
            if (check_tuple_in_testcase(tp_sample[i], testcase_set_[j])){
                ok = true;
                break;
            }
        }
        if (!ok){
            sample_uncovered.emplace_back(i);
        } else {
            ++num_tuple_;
        }
    }
}

void SLSTestcaseSampler::read_cnf_header(ifstream& ifs, int& nvar, int& nclauses){
    string line;
    istringstream iss;
    int this_num_var;
    
    while (getline(ifs, line)){
        if (line.substr(0, 1) == "c")
            continue;
        else if (line.substr(0, 1) == "p"){
            string tempstr1, tempstr2;
            iss.clear();
            iss.str(line);
            iss.seekg(0, ios::beg);
            iss >> tempstr1 >> tempstr2 >> this_num_var >> nclauses;
            if (nvar < this_num_var)
                nvar = this_num_var;
            break;
        }
    }
}

bool SLSTestcaseSampler::read_cnf(){
    ifstream fin(cnf_file_path_.c_str());
    if (!fin.is_open()) return false;

    read_cnf_header(fin, num_var_, num_clauses_);

    if (num_var_ <= 0 || num_clauses_ < 0){
        fin.close();
        return false;
    }

    clauses.resize(num_clauses_);

    for (int c = 0; c < num_clauses_; ++c)
    {
        int cur_lit;
        fin >> cur_lit;
        while (cur_lit != 0)
        {
            int v = abs(cur_lit);
            clauses[c].emplace_back(cur_lit);
            fin >> cur_lit;
        }
    }

    fin.close();

    return true;
}

void SLSTestcaseSampler::GenerateTestCaseSet()
{
    clock_t start_time = clock();
    
    int cur_phase;
    if (t_wise_optimize_ == 1)
        flag_fix_t_wise_optimize_ = false;
    if (!flag_fix_t_wise_optimize_){      //3wise no asf
        if (t_wise_optimize_ == 3)
        {
            cout << "Not use ASF 3wise." << endl;
            cur_phase = 3;
            p_init_tuple_info_ = &SLSTestcaseSampler::Init3TupleInfo;
            p_update_tuple_info_ = &SLSTestcaseSampler::Update3TupleInfo;
            p_get_gain = &SLSTestcaseSampler::get_gain_3_wise;
        }
        else if (cur_phase == 1){
            cur_phase = 1;
            cout << "1wise" << endl;
            p_init_tuple_info_ = &SLSTestcaseSampler::Init1TupleInfo;
            p_update_tuple_info_ = &SLSTestcaseSampler::Update1TupleInfo;
            p_get_gain = &SLSTestcaseSampler::get_gain_1_wise;
        }
    }
    Init();
    if(num_var_ >= 10000){
        switch(t_wise_optimize_){
            case 2:
                parallel_num = 8;
                break;
            case 3:
                parallel_num = 6;
                break;
            case 4:
                parallel_num = 4;
                break;
            case 5:
                parallel_num = 3;
                break;
            case 6:
                parallel_num = 2;
                break;
            default:
                break;
        }
        cout << "parallel:" << parallel_num << endl;
    }

    tuple_solver = new CaDiCaL::Solver;
    tuple_solver->read_dimacs(cnf_file_path_.c_str(), num_var_);

    if ((t_wise_optimize_ == 2 && num_combination_all_possible_ <= sample_cnt_)||(t_wise_optimize_ == 2 && !flag_fix_t_wise_optimize_))
        parallel_num = 2;

    //update begin

    read_cnf();
    cdcl_sampler = new ExtMinisat::SamplingSolver(num_var_, clauses, seed_, true, 0); 
    cout << "use_cdcl:" << use_cdcl << endl;
    cdcl_sampler_list.emplace_back(cdcl_sampler);
    for (int i = 1; i < parallel_num; i++){
        cdcl_sampler_list.emplace_back(new ExtMinisat::SamplingSolver(num_var_, clauses, seed_*i, true, 0));
    }

    //update end

    int stage_c = 0;


    if (flag_fix_t_wise_optimize_){
        cur_phase = t_wise_optimize_;
        if (cur_phase > 2 || num_combination_all_possible_ > sample_cnt_){
        // if (cur_phase > 2 ){
            stage_change(cur_phase);
            stage_c = 1;
            clock_t stage_change_time = clock();
            double cpu_time_pre = GetCPUTime(start_time, stage_change_time);
            cout << "c tuple sampling finished ..." << endl;
            cout << "c current time: " << cpu_time_pre << endl;
        }
    } else {
        if (t_wise_optimize_ == 2){
            cur_phase = 2;
        }
        else if (t_wise_optimize_ == 3 || t_wise_optimize_ == 1);
        else {
            cout << "Not use ASF twise error!" << endl;
            exit(0);
        }
    }

    for (num_generated_testcase_ = 1; num_generated_testcase_ < testcase_set_size_; num_generated_testcase_++)
    {
        if (num_generated_testcase_ > 1 || (cur_phase == 1 && stage_c != 1) || (cur_phase == 2 && stage_c != 1) || (cur_phase == 3 && stage_c != 1)){
            (this->*p_update_tuple_info_)();
        }
        if (num_generated_testcase_ > 1) clear_pq();
        cout << num_generated_testcase_ << ": " << num_tuple_ << endl;
        (this->*p_update_sample_weight_)();
        if (parallel_num == 1)
            GenerateTestcase();
        else {
            vector<pair<int, int> > res;
            res.reserve(num_var_);
            for (int i = 0; i < num_var_; ++i){
                int v1 = (num_generated_testcase_ - var_positive_appearance_count_[i]);
                int v2 = var_positive_appearance_count_[i];
                res.emplace_back(v1, v2);
            }
            int div = (candidate_set_size_) / parallel_num;
            pthread_t tid[parallel_num] = {0};
            for (int i = 0; i < parallel_num; i++){
                int *p = (int*)malloc(sizeof(*p));
                *p = i;
                struct para *myp = (struct para*)malloc(sizeof(struct para));
                myp->cdcl_sampler = cdcl_sampler_list[*p];
                myp->sample_prob = &res;
                myp->candidate_testcase_set_ = &candidate_testcase_set_;
                myp->start = *p * div;
                if( *p == parallel_num - 1 )myp->end = candidate_set_size_;
                else myp->end = (*p + 1) * div;
                pthread_create(&tid[*p],NULL,pthread_func3,(void*)myp);
                free(p);
            }
            for (int i = 0; i < parallel_num; i++)pthread_join(tid[i],NULL);

            selected_candidate_index_ = SelectTestcaseFromCandidateSetByTupleNum2();
            int testcase_index = num_generated_testcase_;
            testcase_set_[testcase_index] = pq[selected_candidate_index_].first;
            selected_candidate_bitset_ = pq[selected_candidate_index_].second;
        }
    }
    clock_t end_time = clock();
    (this->*p_update_tuple_info_)();
    cpu_time_ = GetCPUTime(start_time, end_time);
    cout << "c Generate testcase set finished!" << endl;
    cout << "c CPU time cost by generating testcase set: " << cpu_time_ << " seconds" << endl;
    SaveTestcaseSet(testcase_set_save_path_);

    cout << "c " << cur_phase << "-tuple number of generated testcase set: " << num_tuple_ << endl;

    remove_temp_files();

    delete tuple_solver;
    delete cdcl_sampler;
    for (int i = 1; i < parallel_num; i++)delete(cdcl_sampler_list[i]);

}

void SLSTestcaseSampler::GenerateTestCaseSetFirstReachGiven2TupleNum(long long num_2_tuple_given)
{
    cpu_time_ = 0.;
    clock_t start_time = clock();
    Init();
    clock_t end_time = clock();
    cpu_time_ += GetCPUTime(start_time, end_time);

    for (num_generated_testcase_ = 1; num_generated_testcase_ < testcase_set_size_; num_generated_testcase_++)
    {
        start_time = clock();
        Update2TupleInfo();
        if (num_tuple_ >= num_2_tuple_given)
        {
            break;
        }
        (this->*p_update_sample_weight_)();
        GenerateTestcase();
        end_time = clock();
        cpu_time_ += GetCPUTime(start_time, end_time);
    }

    Update2TupleInfo();

    if (num_tuple_ >= num_2_tuple_given)
    {
        cout << "c Generated testcase set has reached given 2-tuple number!" << endl;
        cout << "c Generated testcase number: " << num_generated_testcase_ << endl;
        cout << "c CPU time cost by generating testcase set: " << cpu_time_ << " seconds" << endl;
        cout << "c Given 2-tuple number: " << num_2_tuple_given << endl;
    }
    else
    {
        cout << "c Generated testcase set failed to reached given 2-tuple number!" << endl;
        cout << "c Generated testcase number: " << num_generated_testcase_ << endl;
        cout << "c CPU time cost by generating testcase set: " << cpu_time_ << " seconds" << endl;
        cout << "c Given 2-tuple number: " << num_2_tuple_given << endl;
    }

    SaveTestcaseSet(testcase_set_save_path_);
}

void SLSTestcaseSampler::SaveTestcaseSet(string result_path)
{
    ofstream res_file(result_path);
    for (int i = 0; i < testcase_set_size_; i++)
    {
        for (int v = 0; v < num_var_; v++)
        {
            res_file << testcase_set_[i][v] << " ";
        }
        res_file << endl;
    }
    res_file.close();
    cout << "c Testcase set saved in " << result_path << endl;
}

void SLSTestcaseSampler::Cal2TupleCoverage()
{
    cout << endl
         << "c Calculate all possible 2-tuple number & coverage ..." << endl;

    CaDiCaL::Solver *cadical_solver = new CaDiCaL::Solver;

    cadical_solver->read_dimacs(cnf_file_path_.c_str(), num_var_);

    num_tuple_all_exact_ = num_tuple_all_possible_;

    long long index_tuple;

    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            for (int v_i = 0; v_i < 2; v_i++)
            {
                for (int v_j = 0; v_j < 2; v_j++)
                {
                    index_tuple = Get2TupleMapIndex(i, v_i, j, v_j);
                    if (!already_t_wise.get(index_tuple))
                    {
                        if (v_i)
                        {
                            cadical_solver->assume(i + 1);
                        }
                        else
                        {
                            cadical_solver->assume(-i - 1);
                        }
                        if (v_j)
                        {
                            cadical_solver->assume(j + 1);
                        }
                        else
                        {
                            cadical_solver->assume(-j - 1);
                        }
                        int res = cadical_solver->solve();
                        if (res == 20)
                        {
                            num_tuple_all_exact_--;
                        }
                        cadical_solver->reset_assumptions();
                    }
                }
            }
        }
    }

    coverage_tuple_ = double(num_tuple_) / num_tuple_all_exact_;

    cout << "c All possible 2-tuple number: " << num_tuple_all_exact_ << endl;
    cout << "c 2-tuple coverage: " << coverage_tuple_ << endl;

    delete cadical_solver;
}

SimpleBitSet SLSTestcaseSampler::get_gain_1_wise(const vector<int>& asgn){
    SimpleBitSet res(num_tuple_all_possible_);
    for (int i = 0; i < num_var_ ; i++)
    {
        long long index_tuple = Get1TupleMapIndex(i, asgn[i]);
        if (!already_t_wise.get(index_tuple)){
            res.set(index_tuple);
        }
    }
    return res;
}

SimpleBitSet SLSTestcaseSampler::get_gain_2_wise(const vector<int>& asgn){
    SimpleBitSet res(num_tuple_all_possible_);
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, asgn[i], j, asgn[j]);
            if (!already_t_wise.get(index_tuple)){
                res.set(index_tuple);
            }
        }
    }
    return res;
}

SimpleBitSet SLSTestcaseSampler::get_gain_3_wise(const vector<int>& asgn){ //3wise no asf
    SimpleBitSet res(num_tuple_all_possible_);
    for (int i = 0; i < num_var_; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            for (int k = j + 1; k < num_var_; k++)
            {
                long long index_tuple = Get3TupleMapIndex(i, asgn[i], j, asgn[j], k, asgn[k]);
                if (!already_t_wise.get(index_tuple)){
                    res.set(index_tuple);
                }
            }
        }
    }
    return res;
}

bool SLSTestcaseSampler::check_tuple_in_testcase(const vector<int>& utp, const vector<int>& asgn){
    for (int i = 0; i < sample_dim_; ++i){
        if (asgn[utp[i << 1]] != utp[i << 1 | 1]){
            return false;
        }
    }
    return true;
}

SimpleBitSet SLSTestcaseSampler::get_gain_by_samples(const vector<int>& asgn){
    SimpleBitSet res(num_tuple_all_possible_);
    for (int ii: sample_uncovered){
        if (check_tuple_in_testcase(tp_sample[ii], asgn)) res.set(ii);
    }
    return res;
}

void SLSTestcaseSampler::remove_temp_files(){
    string cmd;
    int ret;
    if (flag_use_cnf_reduction_ && flag_reduced_cnf_as_temp_){
        cmd = "rm " + reduced_cnf_file_path_;
        ret = system(cmd.c_str());
    }
}

string SLSTestcaseSampler::get_random_prefix()
{
    return cnf_instance_name_ + to_string(getpid()) + to_string(rng_2());
}

void SLSTestcaseSampler::tuple_sampler(int dim, int tuplecnt, vector<vector<int> >& tp_sample)
{
    tp_sample.clear();

    vector<int> tp(dim * 2, 0);
    vector<int> vec(num_var_, 0);
    std::iota(vec.begin(), vec.end(), 1);

    if( num_var_ < 10000 ){
        cout << "tuple sampler num_var < 10000 " << endl;
        for (int i = 0; i < tuplecnt; ++i){
            int tries = 0;
            for (; tries < 100000; ++tries){
                std::shuffle(vec.begin(), vec.end(), rng_2);
                for (int j = 0; j < dim; ++j){
                    tp[j << 1] = vec[j] - 1;
                    if (rng_2() % 2) tp[j << 1 | 1] = 1, tuple_solver->assume(vec[j]);
                    else tp[j << 1 | 1] = 0, tuple_solver->assume(-vec[j]);
                }
                
                int res = tuple_solver->solve();
                tuple_solver->reset_assumptions();
                if (res == 10){
                    tp_sample.emplace_back(tp);
                    break ;
                }
            }
            if (tries == 100000){
                printf("Cannot sample for %d tuples!\n", tuplecnt);
                return ;
            }
        }
        cout << "sample over " << endl;
    }
    else{
        cout << "tuple sampler num_var >= 10000 " << endl;
        for (int i = 0; i < tuplecnt; ++i){
            std::shuffle(vec.begin(), vec.end(), rng_2);
            for (int j = 0; j < dim; ++j){
                tp[j << 1] = vec[j] - 1;
                if (rng_2() % 2) tp[j << 1 | 1] = 1;
                else tp[j << 1 | 1] = 0;
            }            
            tp_sample.emplace_back(tp);
        }
        cout << "sample over " << endl;
        return ;
    }
}