# The-Clonalescent
## user manual and guide
--------
## Introduction
These are a collection of Python scripts to perform simulations and estimate population-level parameters under a coalescent model of sexual and clonal reproduction. To visualize the results, R scripts are also included. Broadly, these tools are designed to detect and simulate sex events in the genealogy of clonal populations using the frequency spectrum of genotypes. Specific uses include:
* Inferring the effective rate of sex in populations
* Estimating the expected number of unique individuals
* Calculating the statistics *psi* and *D psi*
* Comparing the log-likelihoods of estimates of *beta* (from the pareto distribution) and *psi* for a given genotype frequency spectrum
* Saving the posterior probability ditribution of *psi* using MCMC sampling
* Performing coalescent simulations to generate genotype frequency spectrums in a population given a set of parameters

<img src="http://static1.squarespace.com/static/54ad6922e4b0ab38fefa18b1/t/5613cd02e4b0dc9c6cde8cee/1444138244362/?format=750w">

## Dependencies
* Python 2.7.6
    * scipy
    * numpy
    * numdifftools
    * other required libraries should be included in most standard Python distributions. All dependencies can be installed using [pip][1]

## Quick Start
* run `python clonalescent.py -h` to view help and test that all dependencies are installed correctly 
* run `python clonalescent.py` from the directory it is located and use options `-i` to specify the name & location of the input file and `-o` to specify the name & location of the program output. A flag, either `-D`, `-P`, `-S` or `-L`, must be used to select the specific analysis for the program to use. These are explained in depth below.

## Detailed Usage
### Coalescent simulations
To perform coalescent simulations, use the `-S` option when running the program.

`python clonalescent.py -i infile.txt -o outfile.sim -S`

An input file containing three parameters is required: *Ne* - The effective population size, *Ge* - The number of unique genotypes, *n* - the number of individuals sampled. The input file contains exactly three lines in the order *Ne*, *Ge*, *n*. For example the following file:

```
   5000
   50
   20
```

would contain the parameters to run a coalescent simulation of 20 individuals in a genealogy with an effective population size of 5000, 50 unique genotypes. Often, it is of interest to compare the results of simulations using a variety of different parameters. For instance, comparing the genotype frequency sepctrum produced from simulating a genealogy with 50 unique genotypes to simulating a genealogy with 500 unique genotypes. It is possible to specify multiple values for each input parameter and simulations will be performed for each combination. Each additional value is placed on the same line and separated with a tab. Again, for example the following file:

```
   5000
   50 250   1500
   20 100   500
```

would specify coalescent simulations using the following combinations of (*Ne*,*Ge*,*n*): (5000,50,20), (5000,50,100), (5000, 50, 500), (5000,250,20), (5000,250,100) ... (5000,1500,500). By default, the program will perform 10 coalescent simulations and use one CPU process. To change either default value, use the following options:

```
   -n, --n_sims      Number of coalescent simulations
   -p, --processes   Number of CPU processes
```

The `-o` specifies the prefix to use for the file that the raw results of the simulations will be written to and will end in the extension `.sim`. Each line of the output file will start with an ID comprised of the input parameters followed by the simulated genotype frequency spectrum (GFS). The total number of lines will be equal to the total number of simulations performed. The `.sim` output file can be used in downstream analyses and visualization in the included R scripts. 

###Estimating the distribution of *Psi* using MCMC sampling
To generate the posterior probability distribution of the *Psi* statistic, use the `-P` option when running the program.

`python clonalescent.py -i infile.txt -o outfile.mcmc -P`

For each GFS included in the input file, this program will write the results of `-n` MCMC runs to `outfile.mcmc` and display the mean, variance and acceptance rate to standard output. The input file has the same format as the coalescent simulation output, where exactly one GFS and ID is contained per line in the following format:

```
ID GFS
```

Each character must be separated by a tab character and the ID cannot contain any white space. An example input file may appear as:

```
30.1.114 4  3  1  1  1  1  1  1  1  
50.5.127 12 5  3  1
```
This file would include two ID's with their corresponding GFS. The GFS contains the counts of each unique genotype and is of length *S*, where *S* is the total number of unique genotypes observed. In this example, the sample ID 30.1.114 would have a total of 9 unique genotypes, with 4 individuals having genotype 1, 3 individuals having genotype 2 and the rest are singletons. The frequencies do not need to be in order. By default, the program will write every sampled value of *Psi* to the output file so further analyses can be performed on this posterior probability distribution. The output will contain two columns, with one value of *Psi* per line in the following format:

```
Psi   ID
```

The output will contain the results for all ID's included as input. Optionally, an estimate of the effective population size, *Ne*, can be included using the `-N` option, and an additional column will be included in the output that contains an estimate of the expected number of unique genotypes. This optional output will have the following format:

```
Psi E_Ng ID
```
If a single value of *Ne* is used as input, the same *Ne* will be used for all ID's in the input file. Alternativley, if a different *Ne* is required for each ID, the values can be input as a comma seperated list, with an *Ne* value for each ID listed in the same order as included in the input GFS file. Using the example above, if the following command line input

```
python clonalescent.py -i infile.txt -o outfile.mcmc -P -N 500,1000
```
would use 500 as the estimate of *Ne* for ID 30.1.114 and 1000 for ID 50.5.127. The outfile can further be used for visualization using the included R scripts. To control parameters for the MCMC runs, the following options may be used:

```
-n, --n_sims      Number of MCMC iterations to run
-p, --processes   Number of CPU processes
-b, --burnin      Burn in for MCMC
-t, --thin        Thinning parameter for MCMC
```

### Estimating the distribution of *D psi* using coalescent simulations
To generate the distribution of the *D psi* statistic, use the `-D` option when running the program.

`python clonalescent.py -i infile.txt -o outfile.mcmc -D`

The input file has the same format as used to estimate the distribution of *Psi*. The output file has a similar format as the distribution of *Psi* output, where each line contains the result of one simulation and all ID's included in the input are included in the same output. Similarly, a value of *Ne* may be included to provide an estimate of the sex rate of the simulation in the second column of the output. The same options as the coalescent simulations may be used:

```
   -n, --n_sims      Number of coalescent simulations
   -p, --processes   Number of CPU processes
```

The output can be used for further analyses and visualiztion using the included R scripts.The output will contain two columns, with one value of *D Psi* per line in the following format:

```
DPsi   ID
```

###Comparing likelihoods of *beta* and *Psi*
It may be of interest to compare the likelihoods of the coalescent model included here (*Psi*) with the pareto distribution based statistic, *beta*, for a GFS. To compare these, use the `-L` option:

`python clonalescent.py -i infile.txt -L`

Here, only an input file is needed containing a GFS for each ID, using the same format as with previous options. The results will be written to standard output and contain one line for each ID. Four columns are included in the following order:

```
ID Beta_slope Lkhd_beta Lkhd_psi
```

## Authors
* Based on Wham, Fuller and LaJeunesse. (2015) The Coalescent Theory of Clonal Organisms. (In prep. but available shortly)
* Python scripts: Zach Fuller (zlf105@psu.edu)
* R scripts: Drew Wham (fcwham@gmail.com)
* For reports of any bugs, questions or comments please email zlf105@psu.edu (Python script) fcwham@gmail.com (Method and Theory)
[1]:https://pypi.python.org/pypi/pip
      
