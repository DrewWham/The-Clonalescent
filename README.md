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

## Authors
* Python scripts: Zach Fuller (zlf105@psu.edu)
* R scripts: Drew Wham (fcwham@gmail.com)
[1]:https://pypi.python.org/pypi/pip
      
