#!/usr/bin/env python

import sys, gzip, argparse
import numpy as np

from multiprocessing import Pool

def bincountarray(arr):
    return np.bincount(arr.reshape(np.prod(arr.shape)))

def maximise_sample_similarity(diploid_array, sample_idx, n_best_scores = 3):
    
    s = sample_idx
    
    nsites = diploid_array.shape[0]
    nsamples = diploid_array.shape[1]
    
    hapA = diploid_array[:,s,0].copy()
    hapB = diploid_array[:,s,1].copy()
    
    #temporary function to get best scores sum
    get_best_scores_sum = lambda scores: sum(sorted(scores.reshape((nsamples-1)*2))[-n_best_scores:])
    
    #similrity with ll other sequencs
    similarityA = hapA[:,None,None] == diploid_array[:,[_s_ for _s_ in range(nsamples) if _s_ != s],:]
    similarityB = hapB[:,None,None] == diploid_array[:,[_s_ for _s_ in range(nsamples) if _s_ != s],:]
    
    switch_made = False
    
    #from left
    #get all potential switch points
    switch_points_left = (i for i in range(1, nsites) if hapA[i] != hapB[i] and (hapA[i] != hapA[i-1] or hapB[i] != hapB[i-1]))
    #we don't have a previous switch point, so we start with -1. This is needed for summing scores
    previous_i = -1
    #current scores are zero
    current_scoresA = np.zeros((nsamples-1, 2), dtype="int")
    current_scoresB = np.zeros((nsamples-1, 2), dtype="int")
    
    for i in switch_points_left:
        #print(i, file=sys.stderr)
        
        intervening_scoresA = similarityA[previous_i+1:i,:,:].sum(axis=0)
        intervening_scoresB = similarityB[previous_i+1:i,:,:].sum(axis=0)
        
        current_scoresA += intervening_scoresA
        current_scoresB += intervening_scoresB
        
        #sum the best three haplotype matches, assuming that haplotypes should be present multiple times in the dataset
        top_sum_for_staying = get_best_scores_sum(current_scoresA + similarityA[i,:,:]) + get_best_scores_sum(current_scoresB + similarityB[i,:,:])
        
        top_sum_for_switching = get_best_scores_sum(current_scoresA + similarityB[i,:,:]) + get_best_scores_sum(current_scoresB + similarityA[i,:,:])
        
        if (top_sum_for_switching > top_sum_for_staying):
            hapA_old = hapA.copy()
            similarityA_old = similarityA.copy()
            
            #update scores
            current_scoresA += similarityB[i,:,:]
            current_scoresB += similarityA[i,:,:]
            
            #switch them from i onwards!
            hapA[i:] = hapB[i:]
            hapB[i:] = hapA_old[i:]
            
            #and switch the similarities
            similarityA[i:,:,:] = similarityB[i:,:,:]
            similarityB[i:,:,:] = similarityA_old[i:,:,:]
            switch_made = True
        
        else:
            #update scores
            current_scoresA += similarityA[i,:,:]
            current_scoresB += similarityB[i,:,:]
        
        previous_i = i
    
    #now from right
    switch_points_right = (i for i in reversed(range(nsites-1)) if hapA[i] != hapB[i] and (hapA[i] != hapA[i+1] or hapB[i] != hapB[i+1]))
    previous_i = nsites
    current_scoresA = np.zeros((nsamples-1, 2), dtype="int")
    current_scoresB = np.zeros((nsamples-1, 2), dtype="int")
    
    for i in switch_points_right:
        #print(i, file=sys.stderr)
        
        intervening_scoresA = similarityA[i+1:previous_i,:,:].sum(axis=0)
        intervening_scoresB = similarityB[i+1:previous_i,:,:].sum(axis=0)
        
        current_scoresA += intervening_scoresA
        current_scoresB += intervening_scoresB
        
        #sum the best three haplotype matches, assuming that haplotypes should be present multiple times in the dataset
        top_sum_for_staying = get_best_scores_sum(current_scoresA + similarityA[i,:,:]) + get_best_scores_sum(current_scoresB + similarityB[i,:,:])
        
        top_sum_for_switching = get_best_scores_sum(current_scoresA + similarityB[i,:,:]) + get_best_scores_sum(current_scoresB + similarityA[i,:,:])
        
        if (top_sum_for_switching > top_sum_for_staying):
            hapA_old = hapA.copy()
            similarityA_old = similarityA.copy()
            
            #update scores
            current_scoresA += similarityB[i,:,:]
            current_scoresB += similarityA[i,:,:]
            
            #switch them from i onwards!
            hapA[:i+1] = hapB[:i+1]
            hapB[:i+1] = hapA_old[:i+1]
            
            #and switch the similarities
            similarityA[:i+1,:,:] = similarityB[:i+1,:,:]
            similarityB[:i+1,:,:] = similarityA_old[:i+1,:,:]
            switch_made = True
        
        else:
            #update scores
            current_scoresA += similarityA[i,:,:]
            current_scoresB += similarityB[i,:,:]
        
        previous_i = i
    
    return np.column_stack((hapA, hapB)) if switch_made else None


###############################################################################


if __name__ == '__main__':

    ### parse arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="Input file", action = "store", required = True)
    parser.add_argument("-o", "--output_file", help="output file", action = "store", required = True)
    
    parser.add_argument("-t", "--threads", help="Number of threads for parallel processing", type=int, action = "store", default = 1)
    parser.add_argument("--max_iterations", help="Maximum iterations for phase improvement", type=int, action = "store", default = 10)
    parser.add_argument("--n_best_scores", help="Consider the top N best matches for each haplotype", type=int, action = "store", default = 3)
    parser.add_argument("--ignore_first_n_columns", help="Number of columns before first sample", type=int, action = "store", default = 2)
    
    args = parser.parse_args()

    threads = args.threads
    max_iterations = args.max_iterations
    n_best_scores = args.n_best_scores
    ignore_n = args.ignore_first_n_columns
    
    anc = np.loadtxt(args.input_file, skiprows=1, dtype='U10')
    
    with gzip.open(args.input_file, "rt") if args.input_file.endswith(".gz") else open(args.input_file, "rt") as ancfile:
        header = ancfile.readline()
    
    ignore_columns = anc[:,:ignore_n]
    anc = np.delete(anc, tuple(range(ignore_n)), axis=1)
    
    assert anc.shape[1] % 2 == 0, "Odd number of sample columns. Check that that all are diploid, and you have set the option --ignore_first_n_columns correctly."
    
    nsamples = int(anc.shape[1]/2)
    
    nregions = anc.shape[0]
    
    print(f"Preparing to phase {nregions} regions in {nsamples} samples.", file=sys.stderr, flush=True)

    anc = anc.astype(int)
    
    diploid_array = np.zeros(shape=(nregions, nsamples, 2), dtype="int")
    
    diploid_array[:,:,0] = anc[:,range(0,anc.shape[1], 2)]
    diploid_array[:,:,1] = anc[:,range(1,anc.shape[1], 2)]
    
    new_diploid_array = diploid_array.copy()
    
    #multiprocessing pool.map rquires a funcion with one input, so we make a wrapper for the phasing function
    def phase_sample(s):
        return maximise_sample_similarity(new_diploid_array, s, n_best_scores)
    
    switch_made = True
    
    for iteration in range(max_iterations):
        if not switch_made: break
        print("Starting iteration", iteration, file=sys.stderr, flush=True)
        
        switch_made = False
        
        #run phasing in prallel
        with Pool(threads) as p:
            phased_samples = p.map(phase_sample, range(nsamples))
        
        for s in range(nsamples):
            if phased_samples[s] is not None:
                switch_made = True
                new_diploid_array[:, s, :] = phased_samples[s]
    
    #convert back to flat
    output_array = np.column_stack((ignore_columns, new_diploid_array.reshape(nregions, nsamples*2).astype(str)))
    
    np.savetxt(args.output_file, output_array, fmt='%s', delimiter = "\t", header=header.rstrip(), comments="")
