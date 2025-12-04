# AMCS

In this repository, we present the implementation of the EMASS framework, proposed in the paper titled "Efficient Approximation Algorithms for Adaptive Minimum Cost Seed Selection with Improved Guarantees".

## Required Environment
Please run the code on Linux based platforms, where g++ also is installed and the version of C++ should be at least C++ 11.

## Input
### 1. The Dataset
A network example DBLP is given in the path `../dataset`, and all network-related information will be stored in this directory.

### 2. Generate User Costs
The code of generating user costs is in the file `Gene_userInfo.cpp`. Please compile the file by `make Gene_userInfo`.

Then, to generate the costs based on user degree, i.e., $c(u)=0.01(d_{in}(u)+1), \forall u\in V$, we execute

`./Gene_userInfo -Rnd_cost false`.

For the random costs ($\forall u\in V, c(u)\sim U(0,1)$), we execute

`./Gene_userInfo -Rnd_cost true`.

### 3. Generate Realizations
To generate $k$ realizations, please make the file using `make run` under the directory of `run.cpp` Then, the realizations can be derived by simply executing 

`./run -gene_ini_pw true -run_times k`

## Execute the Code
With the above preparatory work done, we are ready to run the EMASS algorithm by making the file with `make run` and executing the code by `./run`. Meanwhile, there is more optional arguments for selection:

$\bullet$ eta_0. $\eta$ in the paper, which defines the influence threshold to be reached.

$\bullet$ eps. $\epsilon$ in the paper, which is the estimation error of mRR-sets.

$\bullet$ run_times. The number of realizations that EMASS needs to be run on.

$\bullet$ Rnd_cost. Specifies which user costs to use, based on degree (`-Rnd_cost false`), or random (`-Rnd_cost true`).

$\bullet$ batch. The number of seeds to be selected in each round of seed selection.

$\bullet$ do_verify. Determines whether to execute the EMASS algorithm, or verify the correctness of mRR-sets revision. 

$\bullet$ MC_round. Determines the number of Monte Carlo simulations in doing verification.

A typical command to execute EMASS is 

`./run -eta_0 10000 -run_times 10 -Rnd_cost false -batch 2 -do_verify false`

The command for conducting correctness verification is

`./run -do_verify true -MC_round 100000`

## Output Format

The format of the output is like

(dblp,eta=10000, cost=1.13, per=1.0, time=589.36)

which means the dataset running on is dblp; the total cost spent is 1.13; the percentage of diffusions reaching $\eta$ is 100%; and the running time is 589.36 seconds.

 <!-- This can also be done by passing the arguments to.
please set the arguments `gene_ini_pw=true` and `run_times=`$k$ in `Argument.h`. Next,  -->