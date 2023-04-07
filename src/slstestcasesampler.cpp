#include "slstestcasesampler.h"
#include <unistd.h>
#include <bits/stdc++.h>

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
        long long index = index_comb * num_combination_all_possible_ + (2 * num_var_ - i - 1) * i / 2 + j - i - 1;
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

void SLSTestcaseSampler::Init2TupleInfo()
{
    num_combination_all_possible_ = num_var_ * (num_var_ - 1) / 2;
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
            int res = already_t_wise.set(static_cast<int>(index_tuple));
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

void SLSTestcaseSampler::GenerateCandidateTestcaseSet()
{
    for (int i = 0; i < candidate_set_size_; i++)
    {
        candidate_sample_init_solution_set_[i] = GetWeightedSampleInitSolution();
        candidate_testcase_set_[i] = GetSatTestcaseWithGivenInitSolution(candidate_sample_init_solution_set_[i]);
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

void SLSTestcaseSampler::GenerateTestcase()
{
    GenerateCandidateTestcaseSet();
    selected_candidate_index_ = SelectTestcaseFromCandidateSetByTupleNum();
    int testcase_index = num_generated_testcase_;
    testcase_set_[testcase_index] = pq[selected_candidate_index_].first;
    selected_candidate_bitset_ = pq[selected_candidate_index_].second;
}

void SLSTestcaseSampler::stage_change(int cur_phase){
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

void SLSTestcaseSampler::GenerateTestCaseSet()
{
    clock_t start_time = clock();
    Init();

    tuple_solver = new CaDiCaL::Solver;
    tuple_solver->read_dimacs(cnf_file_path_.c_str(), num_var_);

    int cur_phase;

    if (flag_fix_t_wise_optimize_){
        cur_phase = t_wise_optimize_;
        if (cur_phase > 2){
            stage_change(cur_phase);
            clock_t stage_change_time = clock();
            double cpu_time_pre = GetCPUTime(start_time, stage_change_time);
            cout << "c tuple sampling finished ..." << endl;
            cout << "c current time: " << cpu_time_pre << endl;
        }
    } else {
        cur_phase = 2;
    }

    for (num_generated_testcase_ = 1; num_generated_testcase_ < testcase_set_size_; num_generated_testcase_++)
    {
        if (num_generated_testcase_ > 1 || cur_phase == 2){
            (this->*p_update_tuple_info_)();
        }
        if (num_generated_testcase_ > 1) clear_pq();
        cout << num_generated_testcase_ << ": " << num_tuple_ << endl;
        (this->*p_update_sample_weight_)();
        GenerateTestcase();
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

SimpleBitSet SLSTestcaseSampler::get_gain_2_wise(const vector<int>& asgn){
    SimpleBitSet res(num_tuple_all_possible_);
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            int index_tuple = static_cast<int>(Get2TupleMapIndex(i, asgn[i], j, asgn[j]));
            if (!already_t_wise.get(index_tuple)){
                res.set(index_tuple);
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
}