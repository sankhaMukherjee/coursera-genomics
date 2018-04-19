import simpleAlgo as sA 
from tqdm import tqdm

def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix

    return D[-1][-1]

def bestMatch(p, t):

    L  = len(p)
    t1 = t[:L]

    minDist = editDistance(t1,p)
    print(minDist)
    for i, m in enumerate(tqdm(t[L:])):
        t1   = t1[1:] + m
        temp = editDistance(t1,p)


        if temp < minDist:
            tqdm.write('{} [{}] --> {:5}'.format(t1, i, temp))
            minDist = temp

        if minDist == 0:
            break

    return minDist, i

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def generateKmerDict(seq, k=3):

    kmer = {}
    for i, s in enumerate(seq):
        
        for j in range(len(s)-k+1):
            t = s[j:k+j]
            if t not in kmer:
                kmer[t] = set([i])
            else:
                kmer[t].add(i)

    return kmer

def main():

    fileName = '../../data/chr1.GRCh38.excerpt.fasta'
    t = sA.readGenome(fileName)


    p  = 'GCTGATCGATCGTACG'
    # Q1 = bestMatch(p, t) # ans = 3
    Q1 = 3

    p  = 'GATTTACCAGATTGAG'
    # Q2 = bestMatch(p, t) # ans = 2
    # Q2 = 2

    fileName = '../../data/ERR266411_1.for_asm.fastq'
    seq, _ = sA.readFastq(fileName)
    print(seq[:5])


    k, mOverlap = 20, 30
    # seq = ['ABCDABCD', 'ABCDEFGH']
    kmer = generateKmerDict(seq, k)
    result = []
    for i, s in enumerate(tqdm(seq)):
        suffix = s[-k:]
        kList  = [m for m in kmer[suffix] if m != i]
        
        for kVal in kList:
            temp = overlap(s, seq[kVal], mOverlap)
            if temp > 0:
                result.append((s, seq[kVal], temp))

    print(result[:10])
    a, _, _ = zip(*result)


    print('Questin 1: {}'.format( Q1 ))
    print('Questin 2: {}'.format( Q2 ))
    print('Questin 3: {}'.format( len(result) ))
    print('Questin 4: {}'.format( len(list(set(a))) ))
    




    return

if __name__ == '__main__':
    main()

