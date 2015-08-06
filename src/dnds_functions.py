
ZERO = 1e-10
import numpy as np
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]



def codon_fitness_to_freqs(codon_fitness):
    ''' 
        Convert amino acid fitness values to stationary codon frequencies using Sella and Hirsh 2005 (Boltzmann).
        
    '''
 
    codon_freqs = np.zeros(61)
    count = 0
    for codon in codons:
        codon_freqs[count] = np.exp( codon_fitness[codon] )
        count += 1
    codon_freqs /= np.sum(codon_freqs)                   
    assert( abs(1. - np.sum(codon_freqs)) <= ZERO), "codon_freq doesn't sum to 1 in codon_fitness_to_freqs"
    return codon_freqs 


def derive_dnds(codon_freqs_dict, mu_dict):
    ''' By default, calculate dS. If no bias and symmetric mutation rates, it will be 1 anyways at virtually no computational cost... '''

    numer_dn = 0.; denom_dn = 0.;
    numer_ds = 0.; denom_ds = 0.;

    for codon in codon_freqs_dict:
        rate = 0.
        sites = 0.
        rate, sites = calc_paths(codon, codon_freqs_dict, mu_dict, 'nonsyn')
        numer_dn += rate
        denom_dn += sites

        rate = 0.
        sites = 0.
        rate, sites = calc_paths(codon, codon_freqs_dict, mu_dict, 'syn')
        numer_ds += rate
        denom_ds += sites

    assert( denom_dn != 0. and denom_ds != 0.), "dN/dS calc indicates no evolution, maybe????"
    return (numer_dn/denom_dn),(numer_ds/denom_ds)
    


     
def calc_paths(source, cfreqs, mu_dict, type):
    ''' type is either syn or nonsyn '''
    rate = 0.
    sites = 0.
    source_freq = cfreqs[source]
    for target in codons:
        diff = get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
        if (type == 'nonsyn' and codon_dict[source] != codon_dict[target]) or (type == 'syn' and codon_dict[source] == codon_dict[target]):
            if len(diff) == 2:
                rate  += calc_subst_prob( source_freq, cfreqs[target], mu_dict[diff], mu_dict[diff[1]+diff[0]] )
                sites += mu_dict[diff]
    rate  *= source_freq
    sites *= source_freq
    return rate, sites



def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff

 
    
def calc_subst_prob(pi, pj, mu_ij, mu_ji):
    if abs(pi) <= ZERO or abs(pj) <= ZERO:
        fixation_rate = 0.
    else:
        p_mu = (mu_ji*pj)/(mu_ij*pi)
        
        # If p_mu == 1, L'Hopitals gives fixation rate of 1 (substitution probability is the forward mutation rate) 
        if abs(1. - p_mu) <= ZERO:
            fixation_rate = 1. 
        else:
            fixation_rate =  (np.log(p_mu))/(1. - (1./p_mu))
    return fixation_rate * mu_ij            










