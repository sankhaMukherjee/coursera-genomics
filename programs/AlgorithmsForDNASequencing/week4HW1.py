import itertools
import simpleAlgo as sA
from tqdm import tqdm 

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def scsAll(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup == []:
            shortest_sup = [sup]
        else:
            if len(sup) == len(shortest_sup[0]):
                shortest_sup.append(sup)  

            if len(sup) < len(shortest_sup[0]):
                shortest_sup = [sup]  # found shorter superstring


    return shortest_sup  # return shortest

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def getGraph(substrings, mOverlap=1):

    graph = {
        'nodes': list(substrings),
        'edges': []
    }

    for i in tqdm(graph['nodes']):
        for j in graph['nodes']:
            if i == j: continue
            oLen = overlap(i, j, mOverlap)

            if oLen > 0:
                graph['edges'].append((oLen, i, j))

    return graph

def updateGraph(graph, oldNodes, newNode, mOverlap=1):

    # drop the old nodes
    for n in oldNodes:
        graph['nodes'].remove(n)

    # drop the old edges
    for e in graph['edges'].copy():
        if (e[1] in oldNodes) or (e[2] in oldNodes):
            graph['edges'].remove( e )
    

    # insert new edges
    for n in graph['nodes']:
        oLen = overlap(n, newNode, mOverlap)
        if oLen > 0:
            graph['edges'].append((oLen, n, newNode))

        oLen = overlap(newNode, n, mOverlap)
        if oLen > 0:
            graph['edges'].append((oLen, newNode, n))

    # insert the new node
    graph['nodes'].append(newNode)

    return graph

def main():


    fileName   = '../../data/ads1_week4_reads.fq'
    substrings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']



    # result = scsAll(substrings)

    # print('Question 1: {}'.format(len(result[0]))) #  Ans = 11
    # print('Question 1: {}'.format(len(result)))    #  Ans = 4

    substrings = ['AAA', 'AAB', 'BBB', 'ABB', 'BBA']
    substrings, _ = sA.readFastq(fileName)

    mLen = 80
    graph = getGraph(substrings, mLen)

    print( '{} -> {}'.format(len(graph['nodes']), len(graph['edges'])))

    while len(graph['edges']) > 0:
        
        oLen, v1, v2 = sorted(graph['edges'], reverse=True)[0]
        v3 = v1[:-oLen] + v2
        
        graph = updateGraph(graph, [v1, v2], v3, mLen)
        tqdm.write( '[{:4}] -> {} -> {}'.format(oLen, len(graph['nodes']), len(graph['edges'])))

        if (mLen > 5) and len(graph['edges']) < 10:
            mLen -= 3

    finalResult = ''.join(graph['nodes'])
    # print(finalResult)
    print(len(finalResult))

    nAs = sum([1 for m in finalResult if m == 'A'])
    nTs = sum([1 for m in finalResult if m == 'T'])

    print(nAs, nTs)

    return

if __name__ == '__main__':
    main()