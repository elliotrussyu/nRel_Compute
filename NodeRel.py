import networkx as nx
import itertools as it
import numpy as np
import igraph as ig
import time 


def combi(N,num):
    """
        Find all comibinations of 'num' elements for int in [0,N-1].
        
        Input: 
            N       total number of elements (numbering starts from zero!)
            num     number of elements in one combination

        Output:
            C       the list with all combinations
    """
    return it.combinations(tuple(range(N)),num)


def combi_all(N):
    """
    Find all comibinations for int in [0,N-1].
    
    Input: 
        N       total number of elements (numbering starts from zero!)

    Output:
        AC       the list with all combinations, arranged from the empty set to N-element sets
    """
    return tuple(map(combi,it.repeat(N),range(N+1)))


def examine_connected_sets(G,N_minus_f):
    c_N_minus_f = 0
    for failed_vertice_set in combi(G.vcount(),N_minus_f):
        GG = G.copy()
        GG.delete_vertices(failed_vertice_set)
        if GG.is_connected() == True:
            c_N_minus_f += 1
    print(c_N_minus_f)
    return c_N_minus_f

def to_ig(H):
    if type(H) is nx.classes.graph.Graph:
        G0 = H.copy()
        N = G0.order()
        mapping = dict(zip(G0, range(N)))
        G0 = nx.relabel_nodes(G0,mapping)
        G = ig.Graph.from_networkx(G0)
    elif type(H) is ig.Graph:
        G = H.copy()
    elif type(H) is list:
        N = len(np.unique(H)) #number of nodes
        G = ig.Graph(N,H)
    elif type(H) is np.ndarray:
        G = ig.Graph.Adjacency(H)
        
    return G


def connected_set_poly(H):
    """
    This function generates the connected set node polynomial of a graph
    Input:
        H       H can be a list of tuples that describes all edges in a graph, 
                or a networkx graph.
    Output:
        s       The node connected set polynomial coefficients, 
                begins with the leading coefficient(highest order) and ends with the constant coefficient(zero order)
    """
    #Processing the input graph, convert input H to an igraph graph
    G = to_ig(H) 
    N = G.vcount()
    
    c_k = np.zeros(N+1) 
    # Initiallize the array of k-ordered connected set coefficient c_k. 
    # Where the numpy array c_k start with the leading coefficient when k = N. 
    # => c_k[0] is c_N, c_k[N+1] is c_(N-1), ...,
    # c_k[N-1] is c_1, and c_k[N] is c_0 which is always 0. 
    if N>0:
        c_k[N-1] = N 
        if N>1:
            c_k[N-2] = G.ecount()
            kappa_det = G.cohesion()
            # Kappa is the vertex connectivity of a graph, then for the first kappa highest order coefficients, the value can be directly calculated via the binomial coefficients. 
            c_k[0:kappa_det] = np.array(list(map(math.comb,it.repeat(N),range(N-kappa_det+1,N+1))))[::-1]
            c_k[kappa_det:N-2] = np.array(list(map(examine_connected_sets,[G.copy() for _ in range(kappa_det,N-2)],range(kappa_det,N-2))))
    return c_k.astype(int)