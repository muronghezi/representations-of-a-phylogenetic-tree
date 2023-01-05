#The file titled "representations.py" contains the python code for the four representations of a phylogenetic tree, plus the conversions in between them. For more details see preprint arXiv:2011.11774v2. In this file after the defined functions, there are also examples on how to call the functions to use them, notated out. 




import itertools
from collections import Counter
import networkx as nx
from ast import literal_eval


def is_partition(partition,n):# "partition" is a list of tuples of numbers. Output is 'yes' if the sets of those tuples form a partition of the set {1,2,...,n}, n is a natural number bigger than 2. 
    N = range(1,n+1)
    S = []
    for part in partition:
        S.append(set(part))
    if (set.union(*S) != set(N) ):
        return 'no'
    while (S != []):
        p1 = S[0]
        i = 1
        while (i < len(S)):
            p2 = S[i]
            i = i+1
            if (p1.intersection(p2) != set()):
                return 'no'
        S.remove(p1)
    return 'yes'    
        



def is_phylogenetic_partition_set(P,n):# P is a list of lists of tuples, representing a set of partitions of the set {1,2,...,n}; n is a natural number bigger than 2. Output 'yes' if P is a phylogenetic partition set on {1,2,...,n}. For detailed definition see Definition 3 in Section 1.2, in the preprint arXiv:2011.11774v2.
    l = []
    l1 = []
    singletons = []
    N = range(1,n+1)
    for partition in P:
        if (is_partition(partition,n) != 'yes'):
            return 'No, since the following set is not a partition on '+str(N)+': '+str(partition)+'.'
        if (len(partition) < 3):
            return 'No, since the following partition does not have at least three parts: '+str(partition)+'.'
        for part in partition:
            l.append(set(part))        
    c = Counter(frozenset(s) for s in l)
    for e in l:
        if c[frozenset(e)] > 1:
            return 'No, since the following part appears more than once: '+str(e)+'.'
        elif (len(frozenset(e)) == 1):
            singletons.append(e)
        else:    
            l1.append(e)    
    if ( set.union(*singletons) != set(N)):
        return 'No, since some one-element subset is not a part of any partition.'
    for e in l1:
        b = 0
        e1 = set(N).difference(e)
        if (e1 not in l1):
            return 'No, since the following part with cardinality bigger than one does not have a complement in some other partition: '+str(e)+'.'
        else: 
            l1.remove(e1)
            l1.remove(e) 
    return 'yes'
            
        
#P = [ [(1,), (2,), (3,), (4,5,6,7,8,9) ], [(1,2,3), (4,), (5,), (6,7,8,9)], [(1,2,3,4,5), (6,7), (8,9)], [(1,2,3,4,5,8,9), (6,), (7,)], [(1,2,3,4,5,6,7), (8,), (9,)] ]

#print is_phylogenetic_partition_set(P,9)      


def is_a_cut(B,n):# B is a list of two tuples, representing a set of two subsets of {1,2,...,n}; n is a natural number bigger than 2. Output 'yes' if the set is a cut. See Definition 7 (for "cut") in Section 1.3 of the preprint arXiv:2011.11774v2.
    N = range(1,n+1)
    if ( ( len(B[0]) <= 1 ) or ( len(B[1]) <= 1 ) ):
        return 'no'
    if ((set(B[0]).intersection(set(B[1])) != set()) or ( set(B[0]).union(set(B[1])) != set(N) )):
        return 'no'
    return 'yes'
    




def is_phylogenetic_cut_set(C,n):# C is a list of lists of two tuples, representing a collection of cuts on {1,2,...,n}; n is a natural number bigger than two. Output 'yes' if the set C is a phylogenetic cut set on {1,2,...,n}. See Definition 7 in Section 1.3 in the preprint arXiv:2011.11774v2 for a "phylogenetic cut set on {1,2,...,n}".
    while (len(C) > 1):
        c = C[0]
        if (is_a_cut(c,n) != 'yes'):
            return 'No, since the following set is not a cut on '+str(N)+': '+str(c)+'.'
        i = 1
        while ( i < len(C) ):
            c1 = C[i]
            i = i+1
            if (c1 != c):
                if ( ( set(c[0]).intersection(set(c1[0])) != set() ) and ( set(c[0]).intersection(set(c1[1])) != set() ) and ( set(c[1]).intersection(set(c1[0])) != set() ) and ( set(c[1]).intersection(set(c1[1])) != set() ) ):
                    return 'no'
        C.remove(c)
    return 'yes'    

#C = [ [(1,2,3), (4,5,6,7,8,9)], [(1,2,3,4,5), (6,7,8,9)], [(1,2,3,4,5,8,9), (6,7)], [(1,2,3,4,5,6,7), (8,9)] ]

#print is_phylogenetic_cut_set(C,9)



def is_phylogenetic_tree(V,E,n):# V is a list of tuples containing information on the label and singletons of nodes, E is a list with cardinality two of tuples representing the edges information. Output 'yes' if the graph (V,E) is a phylogenetic tree with leaf set {1,2,...,n}; n is a natural number bigger than 2. See Definition 2 in Section 1.1 of the preprint arXiv:2011.11774v2 for this concept.
    G = nx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    leaves = []
    N = range(1,n+1)
    if ( (not nx.is_connected(G)) or (nx.number_of_nodes(G) != nx.number_of_edges(G) + 1) ):
        return 'No, it is not a tree.'
    for v in list(G.nodes):
        if (len(G.nodes[v]['singletons']) + G.degree[v] < 3 ):
            return 'no'
        leaves.append(G.nodes[v]['singletons'])
    if ( set.union(*leaves) != set(N) ):
        return 'No, since the union of leaves is not the set '+str(N)+'.'
    return 'yes'
    
#V = [ ('a',{'singletons': {1,2,3} }), ('b',{'singletons': {4,5} }), ('c',{'singletons': {} }), ('d',{'singletons': {6,7} }), ('e',{'singletons': {8,9} }) ]      
#E = [ ('a','b'), ('b','c'), ('c','d'), ('c','e') ]

#print is_phylogenetic_tree(V1,E1,9)

    
def is_partition_on_triples(E,n):# E is a list of lists of tuples with cardinality three, representing a collection of sets of triples, representing a collection of sets of triples; n is a natural number bigger than 2. Output 'yes' if the triples that appear in the input data form a partition of all the triples with labels in {1,2,...,n}.
        N = range(1, n+1)
        G = []
        A = [(x,y,z) for x,y,z in sorted(itertools.combinations(N,3))]
        for S in E:
            for s in S:
                if (len(s) != 3):
                    return 'No, '+str(s)+' is not a triple.'
                elif ( not set(s).issubset(N) ):
                    return 'No, since '+str(s)+' is not a triple on '+str(N)+'.'
                s1 = tuple(sorted(s))
                A.remove(s1)
        if (A != []):
            return 'No, since the following triple(s) is (are) not included in the input data: '+str(A)+'.'
        return 'yes'
        
                
#E = [[(1,2,3),(1,2,4),(1,2,5)], [(1,4,5),(4,5,3)], [ (1,3,4),(1,3,5), (2,3,4), (2,3,5)]]                
        
#print is_partition_on_triples(E,5)  



def is_diverse_triple_set(S,n):# S is a list of tuples with cardinality three, representing a collection of triples; ; n is a natural number bigger than 2. Output 'yes' if the set S is a diverse triple set; output the reason on which axiom being violated so that S is not a diverse triples set. See Section 1.5 of the preprint arXiv:2011.11774v2 for the definition of a "diverse triple set".
    mark = [0,0,0]
    S1 = []
    A = []
    for s in S:
        s1 = tuple(sorted(s))
        S1.append(s1) 
    for a,b,c,d in itertools.combinations(range(1,n+1), 4):
        for x,y,z in itertools.combinations([a,b,c,d],3):
            A.append(tuple(sorted((x,y,z))))
        K = list(set(A) & set(S1))   
        if (len(K) == 1):
            return 'No, since Axiom D1 is not fulfilled for the quadruple '+str({a,b,c,d})+': only one triple '+str(K[0])+' out of the four is in '+str(S)+'.'
        A = []
    for x,y,z in itertools.combinations(range(1,n+1),3):
        t = tuple(sorted((x,y,z)))
        if (t not in S1):
            for s in S1:
                if (set([x,y]).issubset(set(s))):
                    mark[0] = 1
                elif (set([x,z]).issubset(set(s))):
                    mark[1] = 1
                elif (set([y,z]).issubset(set(s))):
                    mark[2] = 1
                if (mark == [1,1,1]):
                    return 'No, since Axiom D2 is not fulfilled at the triple '+str(t)+'.'
        mark = [0,0,0]        
    return 'yes'            
                
#S = [(1,3,4),(1,3,5),(2,3,4),(2,3,5)]
#print is_diverse_triple_set(S,5)
    
      
    


def is_phylogenetic_triple_equivalences(E,n):# E is a list of lists of tuples with cardinality three, representing a collection of sets of triples; n is a natural number bigger than 2. Output 'yes' if E is a phylogenetic equivalence relation on the triples on {1,2,...,n}; otherwise, output the reason why it is not. See Section 1.5 of the preprint arXiv:2011.11774v2 for more details on this concept.
    N = range(1, n+1)
    if (is_partition_on_triples(E,n) != 'yes'):
        return is_partition_on_triples(E,n)
    for S in E:
        if ( is_diverse_triple_set(S,n) != 'yes' ):
            return is_diverse_triple_set(S,n)
    return 'yes'    
    
#E = [[(1,2,3),(1,2,4),(1,2,5)], [(1,4,5),(2,4,5),(3,4,5)], [ (1,3,4), (2,3,4), (2,3,5), (1,3,5)]]

#print is_phylogenetic_triple_equivalences(E,5)


###conversions starts from here###
###conversions starts from here###

def partitions_to_triples(P,n):# P is a list of lists of tuples on {1,2,...,n},representing a set of partitions of the set {1,2,...,n}; n is a natural number bigger than 2. Output is the corresponding phylogenetic equivalent relation on {1,2,...,n} of P, represented as a list of lists of tuples with cardinality three, if P is phylogenetic; otherwise, output "please input a phylogenetic set of partitions". See Definition 3 in Section 1.2 for "phylogenetic set of partitions", see Section 1.5 for "phylogenetic equivalence relation on the triples", and see Definition 41 in Section 2.4 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_partition_set(P,n) != 'yes'):
        return "Please input a phylogenetic set of partitions."
    A = []
    E = []
    for p in P:
        for p1,p2,p3 in itertools.combinations(p,3):
            for x in p1:
                for y in p2: 
                    for z in p3:
                        a = tuple(sorted((x,y,z)))
                        A.append(a)
        E.append(A)  
        A = []
    return E        
                
#P = [ [(1,), (2,), (3,4) ], [(3,), (4,), (1,2)] ]    

#E = partitions_to_triples(P)
#print E

#print is_phylogenetic_triple_equivalences(E,4)
#P = [ [(1,), (2,), (3,), (4,5,6,7,8,9) ], [(1,2,3), (4,), (5,), (6,7,8,9)], [(1,2,3,4,5), (6,7), (8,9)], [(1,2,3,4,5,8,9), (6,), (7,)], [(1,2,3,4,5,6,7), (8,), (9,)] ]           
        
#E2 = partitions_to_triples(P,9)
#print E2
#print is_phylogenetic_triple_equivalences(E2,9)  


def graph_of_a_triple_set(S,n):# S is a list of tuples with cardinality three, representing a collection of triples on the set {1,2,...,n}; n is a natural number bigger than 2. Output its graph. See Definition 41 in Section 2.4 for "the graph of a triple set of {1,2,...,n}", in the preprint arXiv:2011.11774v2. 
    a = 0
    G = nx.Graph()
    for i,j in itertools.combinations(range(1,n+1),2):
        for s in S:
            if ({i,j}.issubset(set(s))):
                a = 1
                break
        if ( a == 0 ):    
            G.add_edge(i,j) 
        a = 0
    return G      
                
    
def triples_to_partitions(E,n):# E is a list of lists of tuples with cardinality three, representing a collection of sets of triples, representing a collection of sets of triples; n is a natural number bigger than 2. Output is the corresponding phylogenetic set of partitions on {1,2,...,n} of P, represented by a list of lists of tuples, if E is phylogenetic; otherwise output "please input a phylogenetic equivalence relation on triples". See Definition 3 in Section 1.2 for "phylogenetic set of partitions", see Section 1.5 for "phylogenetic equivalence relation on the triples", and see Definition 41 in Section 2.4 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_triple_equivalences(E,n) != 'yes'):
        return 'Please input a phylogenetic equivalence relation on triples.'
    M = []
    P = []
    N = range(1,n+1)
    for S in E:
        G = graph_of_a_triple_set(S,n)
        A = sorted(nx.connected_components(G))
        for a in A:
            M.append(tuple(a))
        B = set.union(*A)
        C = list(set(N).difference(B))
        while (len(C)>0):
            C1 = (C[0],)
            M.append(C1)
            C.pop(0)   
        P.append(M)
        M = []
    return P   
            
    
        

#E = [[(1,2,3),(1,2,4),(1,2,5)], [(1,4,5),(2,4,5),(3,4,5)], [ (1,3,4), (2,3,4), (2,3,5), (1,3,5)]]
#P = triples_to_partitions(E, 5)


#print partitions_to_triples(P,5)


def partitions_to_cuts(P,n):# P is a list of lists of tuples, representing a set of partitions on the set {1,2,...,n}; n is a natural number bigger than 2. Output is the corresponding phylogenetic set of cuts on {1,2,...,n} of P, represented as a list of lists with cardinality two of tuples, if P is phylogenetic; otherwise, output "please input a phylogenetic set of partitions". See Definition 3 in Section 1.2 for "phylogenetic set of partitions", see Definition 7 in Section 1.3 for a "phylogenetic cut set on {1,2,...,n}", and see Definition 10 in Section 1.3 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_partition_set(P,n) != 'yes'):
        return 'Please input a phylogenetic set of partitions.'
    A = []
    C = []
    N =range(1,n+1)
    for p in P:
        for x in p:
            if (len(x) != 1):
                A.append(x)
    while (len(A)>0):
        a = A[0]
        a1 = set(N).difference(set(a))
        a2 = tuple(a1)
        C.append([a,a2])
        A.remove(a)
        A.remove(a2)
    return C
 
 
#P = [ [(1,), (2,), (3,), (4,5,6,7,8,9) ], [(1,2,3), (4,), (5,), (6,7,8,9)], [(1,2,3,4,5), (6,7), (8,9)], [(1,2,3,4,5,8,9), (6,), (7,)], [(1,2,3,4,5,6,7), (8,), (9,)] ]

#print is_phylogenetic_cut_set(partitions_to_cuts(P,9),9)

def tree_to_partitions(V,E,n):# V is a list of tuples containing information on the label and singletons of nodes, E is a list with cardinality two of tuples representing the edges information; n is a natural number bigger than 2. Output is the corresponding phylogenetic set of partitions on {1,2,...,n} of (V,E), represented by a list of lists of tuples, if (V,E) represents a phylogenetic tree with leaf set {1,2,...,n}; otherwise output "please input a phylogenetic tree". See Definition 2 in Section 1.1 for "phylogenetic tree", see Definition 3 in Section 1.2 for "phylogenetic set of partitions", and see Definition 5 in Section 1.1 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_tree(V,E,n) != 'yes'):
        return 'Please input a phylogenetic tree.'
    A = set()
    B = []
    P = []
    G = nx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    for v in V:
        sin = G.nodes[str(v[0])]['singletons']
        if (sin != []):
            for s in sin:
                B.append((s,))
        L = G.edges(str(v[0]))
        G.remove_edges_from(L)
        C = sorted(nx.connected_components(G))
        C.remove(set([str(v[0])]))       
        for c in C:
            for x in c:
                a = G.nodes[x]['singletons']
                if (a != set()):
                    A = A.union(a) 
            A1 = tuple(A)
            B.append(A1)
            A = set()
        P.append(B)    
        B = []
        G.add_edges_from(E)
    return P        

#V1 = [('(8, 9)', {'singletons': {8,9}}), ('(4, 5, 6, 7, 8, 9)', {'singletons': {4,5}}), ('(6, 7)', {'singletons': {6,7}}), ('(6, 7, 8, 9)', {'singletons': {}}), ('(1, 2, 3)', {'singletons': {1,2,3}})]
#E1 = [('(8, 9)', '(6, 7, 8, 9)'), ('(4, 5, 6, 7, 8, 9)', '(1, 2, 3)'), ('(4, 5, 6, 7, 8, 9)', '(6, 7, 8, 9)'), ('(6, 7)', '(6, 7, 8, 9)')]

#V = [ ('a',{'singletons': {1,2,3} }), ('b',{'singletons': {4,5} }), ('c',{'singletons': {} }), ('d',{'singletons': {6,7} }), ('e',{'singletons': {8,9} }) ]     
#E = [ ('a','b'), ('b','c'), ('c','d'), ('c','e') ]   

#print tree_to_partitions(V1,E1,9)  


def tree_to_triples(V,E,n):# V is a list of tuples containing information on the label and singletons of nodes, E is a list with cardinality two of tuples representing the edges information; n is a natural number bigger than 2. Output is the corresponding phylogenetic set of partitions on {1,2,...,n} of (V,E), represented by a list of lists of tuples, if (V,E) is a phylogenetic tree with leaf set {1,2,...,n}; otherwies, output "please input a phylogenetic tree". See Definition 2 in Section 1.1 for "phylogenetic tree", see Section 1.5 for "phylogenetic equivalence relation on the triples", and see Definition 5 in Section 1.1 for the conversion from a phylogenetic tree to a set of partitions and see Definition 41 in Section 2.4 for the conversion from a phylogenetic partition set to an equivalent relation on triples --- this function implemented the conversion that goes through these two conversions --- in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_tree(V,E,n) != 'yes'):
        return 'Please input a phylogenetic tree.'
    G = nx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    P = tree_to_partitions(V,E,n)
    T = partitions_to_triples(P,n)
    return T
    
#V = [ ('a',{'singletons': {1,2} }), ('b',{'singletons': {4} }), ('c',{'singletons': {3,5} }) ]     
#E = [ ('a','b'), ('b','c') ]    

#print tree_to_triples(V,E,5) 


def tree_to_cuts(V,E,n):# V is a list of tuples containing information on the label and singletons of nodes, E is a list with cardinality two of tuples representing the edges information; n is a natural number bigger than 2. Output is the corresponding phylogenetic set of cuts on {1,2,...,n} of (V,E), represented by a list of lists with cardinality two of tuples, if (V,E) is a phylogenetic tree with leaf set {1,2,...,n}; otherwise, output "please input a phylogenetic tree". See Definition 2 in Section 1.1 for "phylogenetic tree",see Definition 7 in Section 1.3 for a "phylogenetic cut set on {1,2,...,n}", and see Definition 10 in Section 1.3 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_tree(V,E,n) != 'yes'):
        return 'Please input a phylogenetic tree.'
    N = range(1,n+1)
    A = set()
    B = []
    C = []
    G = nx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    for e in E:
        G.remove_edge(str(e[0]),str(e[1]))
        D = sorted(nx.connected_components(G))
        d = D[0]
        for x in d:
            a = G.nodes[str(x)]['singletons']
            if (a != set()):
                A = A.union(a) 
        A1 = set(N).difference(A)
        B = [tuple(A), tuple(A1)]
        C.append(B) 
        A = set()
        B = []
        G.add_edge(e[0],e[1])
    return C    
    
#V = [ ('a',{'singletons': {1,2,3} }), ('b',{'singletons': {4,5} }), ('c',{'singletons': {} }), ('d',{'singletons': {6,7} }), ('e',{'singletons': {8,9} }) ]     
#E = [ ('a','b'), ('b','c'), ('c','d'), ('c','e') ]

#print is_phylogenetic_cut_set(tree_to_cuts(V,E,9),9)
#print tree_to_cuts(V,E,9)

def cuts_to_tree(C,n):# C is a list of lists with cardinality two of two tuples, representing a collection of cuts on {1,2,...,n}; n is a natural number bigger than two. Output its corresponding phylogenetic tree if the set C is a phylogenetic cut set on {1,2,...,n}; otherwise, output "please input a phylogenetic set of cuts".. See Definition 2 in Section 1.1 for "phylogenetic tree", see Definition 7 in Section 1.3 for a "phylogenetic cut set on {1,2,...,n}", and see Definition 23 in Section 2.2 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_cut_set(C,n) != 'yes'):
        return 'Please input a phylogenetic set of cuts.'
    G = nx.DiGraph()
    G1 = nx.Graph()
    P = []
    x = C[0][0]
    y = C[0][1]
    del C[0]
    for i in C:
        P.append(i[0])
        P.append(i[1])
    for j in P:
        if (set(j).issubset(set(x))):
            G.add_edge(str(j),str(x))
        else:
            if (set(j).issubset(set(y))):
                G.add_edge(str(j),str(y))
    V = list(G.nodes)
    if (str(x) in V):
        V.remove(str(x))
    else:
        if (str(y) in V):
            V.remove(str(y))
    if (len(V) > 0):
        for c in V:
            for d in V:
                c1 = literal_eval(c)
                d1 = literal_eval(d)
                if ( d1 != c1 and set(d1).issubset(set(c1))):
                    G.add_edge(d, c)
    L = list(G.edges)
    if (len(L) > 0):
        for e in L:
            for f in L:
                if ( f != e and f[0] == e[0]):
                    f1 = literal_eval(f[1])
                    e1 = literal_eval(e[1])
                    if set(f1).issubset(set(e1)):
                        L.remove(e)
                    else:
                        if set(e1).issubset(set(f1)):
                            L.remove(f)                
    G1.add_edges_from(L)
    G1.add_edge(str(x),str(y))
    V1 = G1.nodes
    D = [x for x in V1 if G1.degree(x) == 1]
    for v in D:
        G1.nodes[v]['singletons'] = set(literal_eval(v))
    V2 = list( set(V1).difference(set(D)) )    
    for v in V2:
        v1 = v
        v0 = literal_eval(v)
        for x in G1[v]:
            if (set(literal_eval(x)).issubset(set(v0))):
                v0 = tuple(set(v0).difference(set(literal_eval(x))))
        G1.nodes[v1]['singletons'] = set(v0)
    return G1


#C = [ [(1,2,3), (4,5,6,7,8,9)], [(1,2,3,4,5), (6,7,8,9)], [(1,2,3,4,5,8,9), (6,7)], [(1,2,3,4,5,6,7), (8,9)] ]

#G = cuts_to_tree(C,9)


#print tree_to_cuts(G.nodes.data(),G.edges,9)


def partitions_to_tree(P,n):# P is a list of lists of tuples, representing a set of partitions on the set {1,2,...,n}; n is a natural number bigger than 2. Output is the corresponding phylogenetic tree of P, if P is phylogenetic; otherwise, output "please input a phylogenetic set of partitions". See Definition 3 in Section 1.2 for "phylogenetic set of partitions", see Definition 2 in Section 1.1 for "phylogenetic tree", and see Definition 20 in Section 2.1 for this conversion, in the preprint arXiv:2011.11774v2.
    if (is_phylogenetic_partition_set(P,n) != 'yes'):
        return 'Please input a phylogenetic set of partitions.'
    N = range(1, n+1)
    G = nx.Graph()
    A = []
    for p in P:
        G.add_node(str(p))
        for x in p:
            if (len(x) == 1):
                A.append(x[0])
            else:
                y = tuple(set(N).difference(set(x)))
                for q in P:
                    if ((q != p) and (y in q)):
                        for y in q:
                            G.add_edge(str(p),str(q))
        G.node[str(p)]['singletons'] = set(A)
        A = []
    return G
            
            
#P = [ [(1,), (2,), (3,), (4,5,6,7,8,9) ], [(1,2,3), (4,), (5,), (6,7,8,9)], [(1,2,3,4,5), (6,7), (8,9)], [(1,2,3,4,5,8,9), (6,), (7,)], [(1,2,3,4,5,6,7), (8,), (9,)] ]
#G = partitions_to_tree(P,9)

#V = G.nodes.data()
#E = G.edges
 
#C = [ [(1,2,3), (4,5,6,7,8,9)], [(1,2,3,4,5), (6,7,8,9)], [(1,2,3,4,5,8,9), (6,7)], [(1,2,3,4,5,6,7), (8,9)] ] 
 
#print is_phylogenetic_cut_set(tree_to_cuts(V,E,9),9)

#print is_phylogenetic_partition_set(tree_to_partitions(V,E,9),9)

#print tree_to_triples(V,E,9)
#print is_phylogenetic_triple_equivalences( tree_to_triples(V,E,9),9)

#E1 = tree_to_triples(V,E,9)


### some extra functions start from here###
### some extra functions start from here###

def two_same_lists_of_tuples(L,P):# L, P are both lists of tuples. Output 'yes' if they are the same, when considered as two sets of sets. 
    L1 = []
    P1 = []
    for l in L:
        L1.append(tuple(sorted(l)))
    L2 = sorted(L1)
    for p in P:
        P1.append(tuple(sorted(p)))
    P2 = sorted(P1)
    if (L2 !=  P2):
        return 'no'
    return 'yes'

#L = [(3,2,1),(6,7,6)]
#P = [(6,6,7),(1,2,3)]
#print two_same_lists_of_tuples(L,P)


def two_same_partitions(P,Q):# P, Q are both lists of lists of tuples, representing two sets of partitions respectively. Output 'yes' if they are the same, when considered as sets of partitions (the elements of which are also sets).
    A = []
    B = []
    P1 = []
    Q1 = []
    for p in P:
        for x in p:
            A.append(tuple(sorted(x)))
        P1.append(sorted(A))
        A = []
    P2 = sorted(P1)
    for q in Q:
        for y in q:
            B.append(tuple(sorted(y)))
        Q1.append(sorted(B))
        B = []
    Q2 = sorted(Q1)  
    if (P2 != Q2):
        return 'no'
    return 'yes'
            
#P = [ [(1,), (2,), (3,), (4,5,6,7,8,9) ], [(1,2,3), (4,), (5,), (6,7,8,9)], [(1,2,3,4,5), (6,7), (8,9)], [(1,2,3,4,5,8,9), (6,), (7,)], [(1,2,3,4,5,6,7), (8,), (9,)] ]

#Q = [  [(1,2,3), (4,), (5,), (6,7,8,9)],[(1,), (2,), (4,5,8,9,7,6), (3,) ], [(1,2,3,4,5), (6,7), (8,9)], [(1,2,3,4,5,8,9), (6,), (7,)], [(1,2,3,4,5,6,7), (8,), (9,)] ]
            
#print two_same_partitions(P,Q)            
      
    
def two_same_cuts(C,U):# C, U are both lists of lists with cardinality two of tuples, representing two sets of cuts respectively. Output 'yes' if they are the same, when considered as sets of cuts (the elements of which are also sets).
    return two_same_partitions(C,U)
    
#C = [ [(1,2,3), (4,5,6,7,8,9)], [(1,2,3,4,5), (6,7,8,9)], [(1,2,3,4,5,8,9), (6,7)], [(1,2,3,4,5,6,7), (8,9)] ]     
    
#U = [ [ (4,5,6,7,8,9),(3, 1, 2)],  [(1,2,3,4,5,8,9), (6,7)], [ (6,7,8,9),(1,2,5,4,3)],[(1,2,3,4,5,6,7), (8,9)] ]     
    
#print two_same_cuts(C,U) 
 
def two_same_triple_sets(T,R):# T, R are both lists of lists of tuples with cardinality three, representing two sets of collection of triples respectively. Output 'yes' if they are the same, when considered as sets of sets of tripled-sets.
    return two_same_partitions(T,R)

#T = [[(1,2,3),(1,2,4),(1,2,5)], [(1,4,5),(2,4,5),(3,4,5)], [ (1,3,4), (2,3,4), (2,3,5), (1,3,5)]]
#R = [ [(1,4,5),(2,4,5),(3,4,5)], [ (1,3,4), (5,3,2),(2,3,4),  (1,3,5)], [(1,2,3),(5,2,1),(1,2,4),]]

#print two_same_triple_sets(T,R)

def two_phylogenetic_trees_isomorphic(V1,E1,V2,E2,n):# V1, V2 are two lists of tuples containing information on the label and singletons of nodes respectively; E1, E2 are two lists with cardinality two of tuples representing the edges information; n is a natural number bigger than 2. Output 'yes' if the two phylogenetic trees are isomorphic; otherwise, output "no" if they are not isomorphic, output "please input two phylogenetic trees" when at least one of the input is not phylogenetic. See Definition 2 in Section 1.1 for "two phylogenetic trees being isomorphic".
    N = range(1,n+1)
    if ((is_phylogenetic_tree(V1,E1,n) != 'yes') or (is_phylogenetic_tree(V2,E2,n) != 'yes')):
        return 'Please input two phylogenetic trees on '+str(N)+'.'
    C1 = tree_to_cuts(V1,E1,n)
    C2 = tree_to_cuts(V2,E2,n)
    return two_same_cuts(C1,C2)
    
#V1 = [ ('a',{'singletons': {1,2,3} }), ('b',{'singletons': {4,5} }), ('c',{'singletons': {} }), ('d',{'singletons': {6,7} }), ('e',{'singletons': {8,9} }) ]     
#E1 = [ ('a','b'), ('b','c'), ('c','d'), ('c','e') ]

#V2 = [('(8, 9)', {'singletons': {8,9}}), ('(4, 5, 6, 7, 8, 9)', {'singletons': {4,5}}), ('(6, 7)', {'singletons': {6,7}}), ('(6, 7, 8, 9)', {'singletons': {}}), ('(1, 2, 3)', {'singletons': {1,2,3}})]
#E2 = [('(8, 9)', '(6, 7, 8, 9)'), ('(4, 5, 6, 7, 8, 9)', '(1, 2, 3)'), ('(4, 5, 6, 7, 8, 9)', '(6, 7, 8, 9)'), ('(6, 7)', '(6, 7, 8, 9)')]

#print two_phylogenetic_trees_isomorphic(V1,E1,V2,E2,9)
