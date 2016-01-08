#Required libraries
from optparse import OptionParser
import math
from scipy import optimize
import numpy as np
from random import random
import random as rm
import collections
import Queue as queue
import multiprocessing as mp
from numdifftools import Hessian
import itertools
import operator as op
from scipy import stats

#Random number seed
rm.seed(100)

#Display usage, defaults, description and eplig
USAGE= """Usage: %prog [options] -i infile.txt"""
OPT_DEFAULTS={'infile':'-', 'n_sims':10, 'threads':1, 'Dpsi':False,'Ne':None, 'Psi':False, 'thin':10}
DESCRIPTION="""Program description: Clonalescent.py
    This python program will estimate"""
EPILOG="""Requirements:
    An input file in the correct GFS format. Multiple genotype frequency spectrums can be used in a single file, with one per line"""

#Get command line options
def get_options(defaults, usage, description='',epilog=''):
    """Get options, print usage text."""
    parser=OptionParser(usage=usage,description=description,epilog=epilog)
    parser.add_option("-i","--infile",action="store",dest="infile",type="string",
                      default=defaults.get('infile'),
                      help='Name of input file')
    parser.add_option("-o","--outfile",action="store",dest="outfile",type="string",
                      default=defaults.get('outfile'),
                      help='Name of output file')
    parser.add_option("-n","--n_sims",action="store",dest="n_sims",type="int",
                      default=defaults.get('n_sims'),
                      help='Number of coalescent simulations or mcmc iterations to run')
    parser.add_option("-p","--processes",action="store",dest="processes",type="int",
                      default=defaults.get('threads'),
                      help='Number of CPU processes')
    parser.add_option("-N","--Ne",action="store",dest="Ne",type="string",
                      default=defaults.get('Ne'),
                      help='Effective population size. If included the sex rate will be calculated for each simulation and output as the second column in the output file from the D-psi estimation. If multiple values, comma seperate')
    parser.add_option("-b","--burnin",action="store",dest="burnin",type="int",
                      default=defaults.get('burnin'),
                      help='Burn in for MCMC. Default is 10% of total simulations')
    parser.add_option("-t","--thin",action="store",dest="thin",type="int",
                      default=defaults.get('thin'),
                      help='Thinning parameter for MCMC. Default is 10')
    parser.add_option("-D","--Dpsi",action="store_true",dest="Dpsi",
                      default=defaults.get('Dpsi'),
                      help='Flag to output the results of "n" coalescent simulations containing the output of the D-psi statistic. Default false')
    parser.add_option("-P","--Psi",action="store_true",dest="Psi",
                      default=defaults.get('Psi'),
                      help='Flag to output the results of "n" MCMC runs containing the posterior probability distribution of Psi. Default false')
    parser.add_option("-S","--Sim",action="store_true",dest="Sim",
                      default=defaults.get('Sim'),
                      help='Flag to perform simulations of all combinations contained in the input file. For each parameter combination, the simulations will be performed "n" times.Output file will contain the results of all simulations.')
    parser.add_option("-L","--Lklh",action="store_true",dest="Lklh",
                      default=defaults.get('Lklh'),
                      help='Flag to generate only the log likelihoods of beta and psi for each GFS in input file.')

    (options,args)=parser.parse_args()
    
    return (options, args)

#Read input file containing data for parameter estimation
def read_input(infile):
    names, gfs = [], []
    with open(infile, mode="rU") as fp:
        for line in fp:
            line_gfs = []
            line = line.strip("\r\n").split("\t")
            names.append(line[0])
            for x in line[1:]:
                if x:
                    line_gfs.append(int(x))
            gfs.append(line_gfs)
    return names, gfs

#Read input containing parameters for simulations
def read_inparams(infile):
    params = []
    with open(infile, mode="rU") as fp:
        for line in fp:
            params.append(line.strip("\r\n").split("\t"))
    return params

#Get n choose r
def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

#Get theta likelihood
def theta_likelihood(theta, S, J):
    #If any of the values are 0 or negative return likelihood that will get rejected
    if theta <= 0 or S <= 0 or J <= 0:
        return 10000000
    else:
        return -(S * math.log(theta) + math.lgamma(theta) - math.lgamma(theta + J))

#Perform a one-dimensional optimizaiton of a defined function
def optimize_theta(func, S, J):
    return optimize.fminbound(func,1,10000,args=(S, J))

#Return a list containing the descendandts for a node
def find_descendants(node, nodes, descendant_list):
    if nodes[node, 2] == 0:
        descendant_list.append(node)
    else:
        descendant_list = find_descendants(nodes[node,2], nodes, descendant_list)
        descendant_list = find_descendants(nodes[node,3], nodes, descendant_list)
    return descendant_list

#Perform coalescent simulations
def coalescent(sample, theta, sites, rho, rmap, nrun, nTree, sim=False):
    nodes = np.zeros((2*sample-1, 5))
    nodes[:,0] = [x for x in xrange(2*sample-1)]
    k = sample
    klist = [x for x in xrange(1,sample+1)]
    time = 0
    current_node = k + 1
    
    while k > 1:
        rate = float(k*(k-1))/2
        dt = -float(math.log(random()))/rate
        time += dt
        l1 = int(math.ceil(random() * k))
        tmp = klist[l1 -1]
        klist[l1 -1] = klist[k-1]
        klist[k-1] = tmp
        l2 = int(math.ceil(random() * (k-1)))
        nodes[current_node - 1,1] = time
        nodes[current_node - 1,2] = klist[k-1];
        nodes[current_node - 1,3] = klist[l2-1];
        
        nodes[nodes[current_node - 1 ,2],4] = np.random.poisson(theta*(time - float(nodes[nodes[current_node -1 ,2],1]))/2)
        nodes[nodes[current_node - 1 ,3],4] = np.random.poisson(theta*(time - float(nodes[nodes[current_node -1 ,3],1]))/2)
        klist[l2 -1] = current_node

        current_node += 1
        k -= 1
        klist = klist[0:k]
    
    n_mtns=int(sum(nodes[:,4]))
    mut_L = [random() for x in xrange(n_mtns)]
    seqs = np.zeros((sample, n_mtns))
    
    cum_mut = 0
    for node in xrange(2*sample-2):
        if nodes[node,4] > 0:
            who_list = find_descendants(node, nodes, [])
            for mut in xrange(int(nodes[node,4])):
                seqs[who_list,cum_mut] = 1
                cum_mut += 1
        else:
            pass

    x_uniq_cnts = (collections.Counter(map(tuple, seqs))).most_common()
    counts = [x[1] for x in x_uniq_cnts]
    if sim:
        return seqs
    else:
        return counts

#Return a list containing only unique values and removing repeats
def unique(seq, idfun=None): 
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

#Calculate the expected pairwise differences
def ex_pwd(gfs):
    n, S = sum(gfs), len(gfs)
    tpw = float(n*(n-1))/2
    return sum([float(i*(i-1))/2 for i in gfs])/tpw

#Perform a simulation of pairwise differences and compare to expected
def dpi(gfs, n_sims, t_c, q, Ne):
    print "Process number", t_c, "starting simulation"
    pwd_L = []
    n, S = sum(gfs), len(gfs)
    theta = optimize_theta(theta_likelihood, S, n)
    for i in xrange(n_sims):
        s1 = coalescent(n, theta, 10, 0, [], 1, 5)
        pwd_L.append(ex_pwd(s1))
    opi = ex_pwd(gfs)
    ex_dist = float(sum(pwd_L))/len(pwd_L)
    sim_results = [(opi-x)/ex_dist for x in pwd_L]
    if Ne != None:
        sex_results = [nsex(Ne, x) for x in pwd_L]
        sim_results = zip(sim_results,sex_results)
    print "Process number", t_c, "ending simulation"
    q.put(sim_results)
   
def exp_ng(w, Ne):
    #print w, Ne
    return sum([float(w)/(w+x) for x in xrange(Ne-1)])

#Calculate the effective rate of sex    
def nsex(Ne, pi):
    return 1 - (float(Ne * pi)/((Ne*pi*((Ne-1)*pi + 1))**(.5)))

#Base function to perform MCMC estimation of a defined function
def MCMC(initial, iter_num, burn_in, thin, S, J,Ne):
    guess = [initial]
    A = [guess]

    stepsizes = [0.2]
    accepted = 0.0
    fun = lambda x: theta_likelihood(x, S, J)
    hess = -1 * np.array(Hessian(fun)(initial))
    mat = [float(1)/hess]
    for n in xrange(iter_num):
        old_alpha = A[len(A)-1]
        old_loglik = -1 * theta_likelihood(float(old_alpha[0]), S, J)
        new_alpha = np.zeros([len(old_alpha)])
        for i in xrange(len(old_alpha)):
            new_alpha[i] = rm.gauss(old_alpha[i], math.sqrt(-1*mat[0]))
        new_loglik = -1 * theta_likelihood(new_alpha, S, J)

        if (new_loglik > old_loglik):
            A.append(new_alpha)
            accepted = accepted + 1.0
        else:
            u = rm.uniform(0.0,1.0)
            if (u < (math.exp((new_loglik - old_loglik)))):               
                A.append(new_alpha)
                accepted = accepted + 1.0
            else:
                A.append(old_alpha)
    
    print "Acceptance rate = "+str(accepted/iter_num)          
    Clean = []
    for n in xrange(burn_in,iter_num):
        if (n % 10 == 0):
            Clean.append(A[n][0])
    print "Mean:  "+str(np.mean(Clean))
    print "Sigma: "+str(np.std(Clean))
    if Ne != None:
        ng_results = [exp_ng(x, Ne) for x in Clean]
        Clean = zip(Clean,ng_results)
    return Clean

#Calculate the propotion of pairwise differences
def calc_pairwise(site_freq):
    total = []
    for i in itertools.combinations(site_freq,2):
        comp = float(i[0] * i[1])
        total.append(comp/ncr(sum(site_freq),2))
    return float(sum(total))

#Run simulations for all combinations of input parameters
def run_sims(Ne_L,n_L,Ng_L,n_sims,q,t_c):
    sim = []
    print "Starting simulations for process", t_c
    for j in xrange(n_sims):
        for N in Ne_L:
            J = int(N)
            for n_v in n_L:
                n_v = int(n_v)
                for g in Ng_L:
                    g = int(g)
                    S = g
                    psi = optimize_theta(theta_likelihood, S, J)
                    
                    site_freq = coalescent(n_v, psi, 10, 0, [], 1, 5)
                    sites = len(site_freq)
                    output_line = ([".".join([str(n_v), str(g), str(J)]), "\t".join([str(x) for x in site_freq])])
                    #pi = calc_pairwise(site_freq)
                    #print n_v, g, J, pi
                    #print (("%s\n")%("\t".join(output_line)))
                    
                    sim.append(output_line)
    q.put(sim)
    print "Ending simulations for process", t_c
    
def ne_parse(Ne, names):
    Ne_L = []
    if "," in Ne:
        Ne_L = Ne.split(",")
        Ne_L = [int(x) for x in Ne_L]
    else:
        for x in xrange(len(names)):
            Ne_L.append(int(Ne))
    print "Ne values are:", Ne_L
    return Ne_L

def get_beta(gfs):
    n = len(gfs)
    gfs = sorted(gfs)
    uniq_vals = [i for i, e in enumerate(gfs) if gfs.index(e) == i]
    distincts = [math.log(gfs[x]) for x in uniq_vals]
    u_p = [math.log(float(n+1-(i+1))/n) for i in uniq_vals]
    if len(uniq_vals) <= 2:
        slope, intercept, r_value, p_value, std_err, log_lkhd = "NA","NA","NA","NA","NA","NA"
    else:
        slope, intercept, r_value, p_value= stats.linregress(distincts,u_p)[:-1]
        slope = (-1*slope)+1
        log_lkhd = beta_likelihood(gfs,slope)
        slope = slope - 1
    return (slope, intercept, r_value, p_value, log_lkhd)
    
def beta_likelihood(gfs,slope):
    prefactor = math.log(slope-1) - math.log(1)
    return sum([prefactor - slope*(math.log(x)-math.log(1)) for x in gfs])

def main():
    (options,args)=get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    outfile = options.outfile
    n_sims = options.n_sims
    processes = options.processes
    Dpsi = options.Dpsi
    Psi = options.Psi
    Ne = options.Ne
    thin = options.thin
    burnin = options.burnin
    Sim = options.Sim
    Lklh = options.Lklh

    inputs = read_input(infile)
    names, gfs_L = inputs[0], inputs[1]
    
    Ne_L = [None] * len(names)
    if Ne != None:
        Ne_L = ne_parse(Ne,names)
    
    if Lklh:
        for name, gfs in zip(names,gfs_L):
            b_slope, b_intercept, b_r_value, b_p_value, b_log_lkhd = get_beta(gfs)
            n, S = float(sum(gfs)), float(len(gfs))
            theta = optimize_theta(theta_likelihood, S, n)
            psi_lklh = theta_likelihood(theta, S, n)
            print name, b_slope, b_log_lkhd, -psi_lklh
            
    if Dpsi:
        outfile_name = (("%s.dpsi")%(outfile))
        output = open(outfile_name, mode = "w")
        sim_count = 0
        print n_sims, "simulations will be performed for", len(names), "IDs found in input GFS file..."
        for j, gfs in enumerate(gfs_L):
            print "Begining simulations to estimate Dpsi for ID", names[sim_count], "..."
            print "Will use", processes, "processors to perform simulations..."
            sim_per_proc, remain = int(n_sims)/processes, n_sims%processes
            proc_L = [sim_per_proc] * (processes - 1)
            proc_L.append(sim_per_proc + remain)
            
            t_c=0
            processes_L=[]
                
            queues=[mp.Queue() for i in proc_L]
            
            for i in proc_L:
                processes_L.append(mp.Process(target=dpi, args=(gfs, i, t_c,queues[t_c],Ne_L[j])))
                t_c+=1
    
            process_results=[]
            
            #Start the processes in the job queue
            for p in processes_L:
                p.start()
            
            #Return the output from each process in the queue and place in list
            for q in queues:
                process_results.append(q.get())
                
            #Coordinate and join each process
            for p in processes_L:
                p.join()
            print "Joining together final results for", names[sim_count],"..."
            final_results=[item for sublist in process_results for item in sublist]
            
           
            print "Writing the results of", n_sims, "simulations for ID",names[sim_count], "to output file", outfile_name
            if Ne == None:
                for i in final_results:
                    output.write(("%s\t%s\n")%(i,names[sim_count]))
                
            else:
                for i in final_results:
                    output.write(("%s\t%s\t%s\n")%(i[0],i[1],names[sim_count]))
            sim_count += 1
        print "All simulations are completed and written to output"
        
    if Psi:
        outfile_name = (("%s.mcmc")%(outfile))
        output = open(outfile_name, mode = "w")
        sim_count = 0
        if burnin == None:
            burnin = int(n_sims * 0.1)

        for j, gfs in enumerate(gfs_L):
            print "Beginning MCMC sampling of posterior distribution for ID", names[sim_count]
            n, S = float(sum(gfs)), float(len(gfs))
            initial = optimize_theta(theta_likelihood, S, n)
            print "MCMC summary for", names[sim_count], ": "
            print "Ne used for sampling is: ", Ne_L[j]
            post = MCMC(initial, n_sims, burnin, thin, S, n, Ne_L[j])
            print "Writing results to output", outfile_name
            if Ne == None:
                for i in post:
                    output.write(("%s\t%s\n")%(i,names[sim_count]))
                
            else:
                for i in post:
                    #print i
                    output.write(("%s\t%s\t%s\n")%(i[0],i[1],names[sim_count]))
            sim_count += 1

    if Sim:
        outfile_name = (("%s.sim")%(outfile))
        output = open(outfile_name, mode = "w")
        
        params = read_inparams(infile)
        
        Ne_L, Ng_L, n_L = params[0], params[1], params[2]
        
        print "Input Ne values: ", ",".join(Ne_L)
        print "Input Ng values: ", ",".join(Ng_L)
        print "Input n values: ", ",".join(n_L)
        
        sim_per_proc, remain = int(n_sims)/processes, n_sims%processes
        proc_L = [sim_per_proc] * (processes - 1)
        proc_L.append(sim_per_proc + remain)
        
        t_c=0
        processes_L=[]
            
        queues=[mp.Queue() for i in proc_L]
        
        for i in proc_L:
            processes_L.append(mp.Process(target=run_sims, args=(Ne_L, n_L, Ng_L, i, queues[t_c],t_c)))
            t_c+=1
    
        process_results=[]
        
        #Start the processes in the job queue
        for p in processes_L:
            p.start()
        
        #Return the output from each process in the queue and place in list
        for q in queues:
            process_results.append(q.get())
            
        #Coordinate and join each process
        for p in processes_L:
            p.join()
        print "Joining together final simulation results"
        final_results=[item for sublist in process_results for item in sublist]
        #print final_results

        for i in final_results:
            output.write(("%s\n")%("\t".join(i)))

    if Lklh:
        pass
    else:
        print "All simulations are completed and written to output"

if __name__ == '__main__':
    main()
