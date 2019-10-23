def read_fastq(fasta_q):
  with open(fasta_q, 'r') as f:
    for count, line in enumerate(f, start=1):
      if (count/2)%2 == 1: #selectionner seulement les lignes du style 2,6,10 et pas 4 et 8
        yield (line.strip("\n"))


def cut_kmer(seq,k_mer):
    for i in range(len(seq)-k_mer+1):
        yield seq[i:i+k_mer]


def build_kmer_dict(fasta_q,k_mer): # finalement il renvoie ça que pour la derniere seq ! faut t'il conctner les seq?
        dic = {}

        it = read_fastq(fasta_q)

        for sequence in it:
            it_tmp = cut_kmer(sequence,k_mer) #renvoie un itér des seq de taille k_mer
            list_tmp_kmer = list(it_tmp)
            for sous_k_mer in list_tmp_kmer:

                if sous_k_mer in dic:
                    dic[sous_k_mer] += 1
                else:
                      dic[sous_k_mer] = 1


        return dic

def build_graph(dic_kmer) :

    Grph = nx.DiGraph()

    for kmer, valeur in dic_kmer.items():

        nd1 = kmer[:-1]

        nd2 = kmer[1:]

        Grph.add_edge(nd1 , nd2 , weight = valeur)

    return Grph


def get_starting_nodes(graph):

    node_entry = []

    for node in graph.nodes:

        if len(list(graph.predecessors(node))) == 0:

            node_entry.append(node)

    return node_entry

