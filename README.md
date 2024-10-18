# *LS-Sampling-Plus*: An Effective Local Search Based Sampling Approach for Achieving High t-wise Coverage

*LS-Sampling-Plus* is an effective Local search based sampling approach for solving the t-wise coverage maximum (t-wise CovMax) problem. This repository includes the implementation of *LS-Sampling-Plus* and testing benchmarks. We note that *LS-Sampling-Plus* should be run on an operating system of `GNU/Linux (64-bit)`.

## Developers
- Chuan Luo (<chuanluophd@outlook.com>)
- Qiyuan Zhao (<zqy1018@hotmail.com>)
- Binqi Sun (<binqi.sun@tum.de>)

## How to Clone this Repository

```
git clone --recursive https://github.com/chuanluocs/LS-Sampling-Plus.git
```

## How to Build *LS-Sampling-Plus*

```
make
```

## How to Run *LS-Sampling-Plus*

```
./LS-Sampling-Plus -input_cnf_path [BENCHMARK_PATH] <optional_parameters>
```

`-input_cnf_path` is the only required parameter. The input file for *LS-Sampling-Plus* must be in [Dimacs format](http://www.satcompetition.org/2011/format-benchmarks2011.html). The directory named `cnf_benchmarks/` contains all 123 testing benchmarks, which can be taken as input. 

For the optional parameters, we list them as follows. 

| Parameter | Allowed Values | Default Value | Description | 
| - | - | - | - |
| `-output_testcase_path` | string | `./[BENCHMARK_NAME]_testcase_set.txt` | path to which the generated test suite is saved |
| `-seed` | integer | 1 | random seed | 
| `-k` | positive integer | 100 | the size of the test suite | 
| `-lambda` | positive integer | 100 | the number of candidates per iteration | 
| `-use_dynamic_updating_sampling_prob` | 0 or 1 | 1 | 1 if the dynamic mechanism for updating sampling probabilities is enabled, 0 otherwise |\
| `-t_wise` | integer larger than 1 | 2 | testing strength |
| `-delta` | positive integer | 1000000 | the initial cardinality of measuring set (will be ignored when using the exact scoring function) |

## Example Command of Calling *LS-Sampling-Plus*

### Example 1

```
./LS-Sampling-Plus -input_cnf_path cnf_benchmarks/adderII.cnf -seed 1 -k 100 -lambda 100
```

The above example command represents calling *LS-Sampling-Plus* to solve benchmark `cnf_benchmarks/adderII.cnf` by setting the random seed to 1, the allowed number of test cases to 100 and the number of candidates per iteration to 100. The output will be like: 

```
[... some omitted output]
1: 813450
2: 1623245
3: 2007756
4: 2291953
[... some omitted output]
96: 3059575
97: 3059579
98: 3059584
99: 3059589
c Generate testcase set finished!
c CPU time cost by generating testcase set: 19.3237 seconds
c Testcase set saved in ./adderII_testcase_set.txt
c 2-tuple number of generated testcase set: 3059593
```

Here, the lines in the form of `number1: number2` represent the progress of constructing the test suite: for the first `number1` test cases, they have covered `number2` 2-wise tuples. The last line reports that the whole test suite covers 3,059,593 2-wise tuples. 

### Example 2

```
./LS-Sampling-Plus -input_cnf_path cnf_benchmarks/axtls.cnf -seed 2 -k 100 -lambda 100 -t_wise 6 -delta 10000
```

The above example command represents calling *LS-Sampling-Plus* to solve benchmark `cnf_benchmarks/axtls.cnf` by setting the random seed to 2, the allowed number of test cases to 100 and the number of candidates per iteration to 100. Also, the target testing strength is set to 6, and the initial size of measuring set is set to 10000. The output will be like: 

```
[... some omitted output]
c tuple sampling finished ...
c current time: 0.069298
1: 294
2: 598
3: 889
4: 1155
[... some omitted output]
96: 8411
97: 8434
98: 8456
99: 8478
c Generate testcase set finished!
c CPU time cost by generating testcase set: 0.584444 seconds
c Testcase set saved in ./axtls_testcase_set.txt
c 6-tuple number of generated testcase set: 8499
```

Here, again, the lines in the form of `number1: number2` represent the progress of constructing the test suite: for the first `number1` test cases, they have covered `number2` 6-wise tuples in the initial measuring set. The last line reports that the whole test suite covers 8,499 6-wise tuples in the initial measuring set. 

## Implementation of *LS-Sampling-Plus*

- The implementation of *LS-Sampling-Plus* can be found in the directory entitled `src`.


## Testing Benchmarks for Evaluating *LS-Sampling-Plus*

- The directories entitled `cnf_benchmarks/`, `cnfForRQ7`, `cnfForRQ8` and `cnfForRQ10` include all testing benchmarks.
- The file entitled 'benchmark.csv' shows the numbers of variables and clauses for our primary benchmark set of 123 binary benchmarks.


## Experimental Results

- The file entitled 'Appendix.pdf' shows the results of the significance test we conducted, as well as the box plots and scatter plots drawn based on the experimental results.

