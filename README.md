# AMCS

In this repository, we present the implementation of the EMASS framework, proposed in the paper titled "Efficient Approximation Algorithms for Adaptive Minimum Cost Seed Selection via mRR-set Updates".

## Required Environment
Please run the code on Linux based platforms, where g++ is installed and the version of C++ should be at least C++ 11.

## Input
### 1. The Dataset
A network example DBLP is given in the path `../dataset`, and all network-related information will be stored in this directory.

### 2. Generate User Costs
The code of generating user costs is in the file `Gene_userInfo.cpp`. Please compile the file by `make Gene_userInfo`.

Then, the generation of user costs based on the degree (i.e., $c(u)=0.01(d_{in}(u)+1), \forall u\in V$) can be done by

`./Gene_userInfo`.

### 3. Generate Realizations
To generate $k$ realizations, please make the file using `make run` under the directory of `run.cpp`. Then, the realizations can be derived by

`./run -gene_ini_pw true -run_times k`

## Execute the Code
With the above preparatory work done, the EMASS algorithm can be run by making the file with `make run` and executing the code by `./run`. Meanwhile, there is more optional arguments for selection:

$\bullet$ eta_0. $\eta$ in the paper, which defines the influence threshold to be reached.

$\bullet$ eps. $\epsilon$ in the paper, which is the estimation error of mRR-sets.

$\bullet$ run_times. The number of realizations that EMASS needs to be run on.

$\bullet$ batch. The number of seeds to be selected in each round of seed selection.

A typical command to execute EMASS is 

`./run -eta_0 10000 -run_times 10 -batch 2 -eps 0.7`

## Output Format

The format of the output is like

(dblp, eta=10000, cost=102.13, %=1.0, time=1589.36)

which means the dataset running on is dblp; the total cost spent is 102.13; the percentage of diffusions reaching $\eta$ is 100%; and the running time is 1589.36 seconds.