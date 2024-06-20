import sys
import re

def similarity(clusters1, clusters2):

    items1 = {}
    for i, clust in enumerate(clusters1):
        for item in clust:
            items1[item] = i

    items2 = {}
    for i, clust in enumerate(clusters2):
        for item in clust:
            items2[item] = i

    score = 0
    diffs = set()
    for item in items1:
        idx1 = items1[item]
        s1 = clusters1[idx1]
        idx2 = items2[item]
        s2 = clusters2[idx2]
        jaccard = len(s1.intersection(s2)) / len(s1.union(s2))
           
        score += jaccard
        if jaccard < 1:
            diffs.add(idx1)
            diffs.add(idx2)
            #print(jaccard, s1, s2)


    return (score / len(items1), len(diffs))

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("USAGE: cmp.py <file_1> <file_2>")
        sys.exit(-1)

    files = [sys.argv[1], sys.argv[2]]

    print(f'File 1: {files[0]}, File 2: {files[1]}\n')

    clusters = []
    for fid, file in enumerate(files):
        d = {}
        fh = open(file)
        fh.readline()
        for line in fh:
            cols = re.split('[ ,\t]', line)
            genome_id = cols[0]
            cluster_id = cols[1]
            if cluster_id not in d:
                d[cluster_id] = set()
            d[cluster_id].add(genome_id)
        fh.close()
        clusters.append(list(d.values()))

    clusters1, clusters2 = clusters
    sim, n_diffs = similarity(clusters1, clusters2)

    print(f'Similarity: {sim}, #differences: {n_diffs}')

    if (sim > 0.999999):
        sys.exit(0)
    else:
        sys.exit(-1)
