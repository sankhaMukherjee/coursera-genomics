
# ------- Already provided ---------------------------
# we shall change this function at specific locations

'''
    Third, implement versions of the naive exact matching and 
    Boyer-Moore algorithms that additionally count and return 
    
        (a) the number of character comparisons performed and 
        (b) the number of alignments tried. Roughly 

    speaking, these measure how much work the two different 
    algorithms are doing.
'''

import bm_preproc  as bmPpc
import kmer_index  as kmI
import week1HW1    as w1
import colorOutput as cO

def naive(p, t):
    charCompCount = 0
    alignCount    = 0

    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        alignCount += 1

        match = True
        for j in range(len(p)):  # loop over characters
            charCompCount += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, charCompCount, alignCount

def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p 
    """
    charCompCount = 0
    alignCount    = 0

    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                charCompCount += 1
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
        alignCount += 1
    return occurrences, charCompCount, alignCount

def readFasta(fileName):

    with open(fileName) as f:
        f.readline()
        data = ''.join([l.strip() for l in f])

    return data

def confirm(s1, s2, nMismatch):
    n = sum([ a!=b  for a, b in zip(s1, s2)])
    return n <= nMismatch

def main():


    __VERBOSE__ = False
    fileName = '../../data/chr1.GRCh38.excerpt.fasta'
    t = readFasta(fileName)
    kMer = 8 # kmer index
    index = kmI.Index(t, kMer)
    
    p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    p_bm = bmPpc.BoyerMoore(p)
    
    o, cC, aC = naive(p, t)
    print('Question 1: {}'.format(aC))
    print('Question 2: {}'.format(cC))

    o, cC, aC = boyer_moore(p, p_bm, t)
    print('Question 3: {}'.format(aC))

    # print(index.index)
    p = 'GGCGCGGTGGCTCACGCCTGTAAT'

    o, iHits = [], 0
    # Three parts of the pigeonhole problem ...
    for j in range(3):
        
        inner = p[j*kMer : (j+1)*kMer]

        if __VERBOSE__:
            print('This is pegionhole step: {}, {}'.format(j, inner))
            toPrint = [' ']
            for l in range(3):
                toPrint.append( ' '*kMer if l != j else inner )
            print('{:10}   {}'.format( ' ', p ))
            print('{:10}   {}{}{}'.format( *toPrint ))

        hits = index.query( inner )
        iHits += len(hits)

        for m in hits:
            if not confirm(t[m - j*kMer : m  - j*kMer + len(p)], p, 2):
                continue

            if (m - j*kMer) not in o:
                o.append( m - j*kMer )

            if __VERBOSE__:
                print( '{:10} : {}'.format(
                        m, cO.cStr(t[m - j*kMer : m - j*kMer + len(p)], p )) )


    o = sorted(o)
    print('Question 3: {}'.format(len(o)))
    print('Question 4: {}'.format( iHits ))
    if __VERBOSE__:
        print('This is the Indexed match')
        for m in o:
            print( '{:10} : {}'.format(
                m, cO.cStr(t[m:m+len(p)], p )) )


    # ---- Rest of the examples for confirmation --------
    if __VERBOSE__:
        p_bm = bmPpc.BoyerMoore(p)
        o, cC, aC = boyer_moore(p, p_bm, t)
        print(o)
        
        print('This is exact matching with Boyer More')
        for m in o:
            print( '{:10} : {}'.format(
                m, cO.cStr(t[m:m+len(p)], p )) )

            
        print('This is the naive match')
        o = w1.naive_nErr(p,t,2)
        for m in o:
            print( '{:10} : {}'.format(
                m, cO.cStr(t[m:m+len(p)], p )) )

    return

if __name__ == '__main__':
    main()