import sys
import re

def similarity(clusters1, clusters2):

    print(f'Filling first set')
    items1 = {}
    for i, clust in enumerate(clusters1):
        for item in clust:
            items1[item] = i

    print(f'Filling second set')
    items2 = {}
    for i, clust in enumerate(clusters2):
        for item in clust:
            items2[item] = i

    print(f'Iterating over items')
    score = 0
    diffs = set()
    
    num_rows = max(len(clusters1), len(clusters2))
    jaccards = [{} for i in range(num_rows)]
    
    
    for i,item in enumerate(items1):
        if (i > 0 and i % 100000 == 0):
            print(f'\r{i}', end='')
        
        idx1 = items1[item]
        idx2 = items2[item]
                 
        if idx2 in jaccards[idx1]:
            j = jaccards[idx1][idx2]
        else:
            s1 = clusters1[idx1]
            s2 = clusters2[idx2]
            j = len(s1.intersection(s2)) / len(s1.union(s2))
            jaccards[idx1][idx2] = j
            jaccards[idx2][idx1] = j
           
        score += j
        if j < 1:
            diffs.add(idx1)
            diffs.add(idx2)
            #print(jaccard, s1, s2)
    
    print(f'\r{i} done')
      
    return (score / len(items1), len(diffs))

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("USAGE: cmp.py <file_1> <file_2>")
        sys.exit(-1)

    files = [sys.argv[1], sys.argv[2]]

    print(f'File 1: {files[0]}, File 2: {files[1]}\n')

    clusters = []
    for fid, file in enumerate(files):
        print(f'Loading file: {file}')
        d = {}
        fh = open(file)
        fh.readline()
        for i,line in enumerate(fh):
            if (i > 0 and i % 100000 == 0):
                print(f'\r{i}', end='')
            cols = re.split('[ ,\t]', line)
            genome_id = cols[0]
            cluster_id = cols[1]
            if cluster_id not in d:
                d[cluster_id] = set()
            d[cluster_id].add(genome_id)
        fh.close()
        clusters.append(list(d.values()))
        print(f'\r{i} done')

    clusters1, clusters2 = clusters
    sim, n_diffs = similarity(clusters1, clusters2)

    print(f'Similarity: {sim}, #differences: {n_diffs}')

    if (sim > 0.999999):
        sys.exit(0)
    else:
        sys.exit(-1)
