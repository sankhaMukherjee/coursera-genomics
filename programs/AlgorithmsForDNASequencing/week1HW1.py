import simpleAlgo as sA
import numpy as np

def naive_rc(p, t):
    '''naive with reverse complement
    
    [description]
    
    Arguments:
        p {str} -- The pattern that we are looking to check
        t {str} -- The test string within which the pattern
            might occur
    
    Returns:
        list of ints -- The list of locations where the pattern
            is found in the test string
    '''

    p_rc = sA.reverseComplement(p)

    occurences =  sA.naive(p, t)
    occurences += sA.naive(p_rc, t)
    return sorted(list(set(occurences)))

def naive_nErr(p, t, n):
    occurrences = []

    for i in range(len(t) - len(p) + 1):  # loop over alignments
        cErr = 0
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                cErr += 1
                if cErr > n:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def fastqQual(fileName):

    with open(fileName) as f:
        while True:
            f.readline()
            seq = f.readline()
            if len(seq) < 10: break
            f.readline()
            qual = f.readline().rstrip()
            qual = [ (ord(q)-33) for q in qual]
            yield qual

    return

def main():
    lambdaVirus = '../../data/lambda_virus.fa'
    humanGenome = '../../data/ERR037900_1.first1000.fastq'
    humanData   = fastqQual(humanGenome)
    virusData = sA.readGenome(lambdaVirus)

    'Question 1'
    l = len(naive_rc(  'AGGT', virusData))
    print('Question 1: {}'.format(l))
    
    'Question 2'
    l = len(naive_rc('TTAA', virusData))
    print('Question 2: {}'.format(l))

    'Question 3'
    p = min(naive_rc('ACTAAGT', virusData))
    print('Question 3: {}'.format(p))

    'Question 4'
    p = min(naive_rc('AGTCGA', virusData))
    print('Question 4: {}'.format(p))

    'Question 5'
    t = naive_nErr('ACTTTA', 'ACTTACTTGATAAAGT', 2)
    v = naive_nErr('TTCAAGCC', virusData, 2)
    print('Question 5 test: {}'.format(t))
    print('Question 5: {}'.format(len(v)))

    'Question 6'
    v = naive_nErr('AGGAGGTT', virusData, 2)
    print('Question 6: {}'.format(min(v)))

    humanData = np.array(list(humanData))
    dataMean  = humanData.mean(axis=0)
    print('The mean quality scores for each sequence are: ')
    print(dataMean)
    print('Location of the lowest mean:')
    print( np.arange(len(dataMean))[  dataMean == dataMean.min() ]  )
    
    return

if __name__ == '__main__':
    main()