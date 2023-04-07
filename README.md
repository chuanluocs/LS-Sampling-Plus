# *LS-Sampling*: An Effective Local Search Based Sampling Approach for Achieving High t-wise Coverage

*LS-Sampling* is an effective local search based sampling approach for solving the t-wise coverage maximum (t-wise CovMax) problem. This repository includes the implementation of *LS-Sampling* and testing benchmarks. We note that *LS-Sampling* should be run on an operating system of `GNU/Linux (64-bit)`.

## Developers
- Chuan Luo (<chuanluophd@outlook.com>)
- Qiyuan Zhao (<zqy1018@hotmail.com>)
- Binqi Sun (<binqi.sun@tum.de>)

## How to Build *LS-Sampling*

```
make
```

## How to Run *LS-Sampling*

```
./LS-Sampling -input_cnf_path [BENCHMARK_PATH] <optional_parameters>
```

`-input_cnf_path` is the only required parameter. The input file for *LS-Sampling* must be in [Dimacs format](http://www.satcompetition.org/2011/format-benchmarks2011.html). The directory named `cnf_benchmarks/` contains all 123 testing benchmarks, which can be taken as input. 

For the optional parameters, we list them as follows. 

| Parameter | Allowed Values | Default Value | Description | 
| - | - | - | - |
| `-output_testcase_path` | string | `./[BENCHMARK_NAME]_testcase_set.txt` | path to which the generated test suite is saved |
| `-seed` | integer | 1 | random seed | 
| `-k` | positive integer | 100 | the size of the test suite | 
| `-lambda` | positive integer | 100 | the number of candidates per iteration | 
| `-use_formula_simplification` | 0 or 1 | 1 | 1 if the input will be simplified with `bin/coprocessor`, 0 otherwise |
| `-use_dynamic_updating_sampling_prob` | 0 or 1 | 1 | 1 if the dynamic mechanism for updating sampling probabilities is enabled, 0 otherwise |
| `-use_diversity_aware_heuristic_search` | 0 or 1 | 1 | 1 if the diversity-aware heuristic search is enabled, 0 otherwise |
| `-t_wise` | integer larger than 1 | 2 | testing strength |
| `-delta` | positive integer | 1000000 | the initial cardinality of measuring set (will be ignored when using the exact scoring function) |

## Example Command of Calling *LS-Sampling*

### Example 1

```
./LS-Sampling -input_cnf_path cnf_benchmarks/freebsd-icse11.cnf -seed 1 -k 100 -lambda 100 -use_formula_simplification 1 -use_dynamic_updating_sampling_prob 1 -use_diversity_aware_heuristic_search 1
```

The above example command represents calling *LS-Sampling* to solve benchmark `cnf_benchmarks/freebsd-icse11.cnf` by setting the random seed to 1, the allowed number of test cases to 100 and the number of candidates per iteration to 100. The output will be like: 

```
[... some omitted output]
1: 973710
2: 1927320
3: 2358585
4: 2685594
[... some omitted output]
96: 3733704
97: 3733919
98: 3734093
99: 3735170
c Generate testcase set finished!
c CPU time cost by generating testcase set: 63.2704 seconds
c Testcase set saved in ./freebsd-icse11_testcase_set.txt
c 2-tuple number of generated testcase set: 3735341
```

Here, the lines in the form of `number1: number2` represent the progress of constructing the test suite: for the first `number1` test cases, they have covered `number2` 2-wise tuples. The last line reports that the whole test suite covers 3,735,341 2-wise tuples. 

### Example 2

```
./LS-Sampling -input_cnf_path cnf_benchmarks/axtls.cnf -seed 2 -k 100 -lambda 100 -use_formula_simplification 1 -use_dynamic_updating_sampling_prob 1 -use_diversity_aware_heuristic_search 1 -t_wise 6 -delta 10000
```

The above example command represents calling *LS-Sampling* to solve benchmark `cnf_benchmarks/axtls.cnf` by setting the random seed to 2, the allowed number of test cases to 100 and the number of candidates per iteration to 100. Also, the target testing strength is set to 6, and the initial size of measuring set is set to 10000. The output will be like: 

```
[... some omitted output]
c tuple sampling finished ...
c current time: 0.087385
1: 294
2: 596
3: 882
4: 1161
[... some omitted output]
96: 8328
97: 8350
98: 8371
99: 8391
c Generate testcase set finished!
c CPU time cost by generating testcase set: 1.02321 seconds
c Testcase set saved in ./axtls_testcase_set.txt
c 6-tuple number of generated testcase set: 8411
```

Here, again, the lines in the form of `number1: number2` represent the progress of constructing the test suite: for the first `number1` test cases, they have covered `number2` 6-wise tuples in the initial measuring set. The last line reports that the whole test suite covers 8,411 6-wise tuples in the initial measuring set. 

## Implementation of *LS-Sampling*

- The implementation of *LS-Sampling* can be found in the directory entitled `src`.


## Testing Benchmarks for Evaluating *LS-Sampling*

- The directory entitled `cnf_benchmarks/` includes all testing benchmarks.
