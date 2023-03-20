import networkx as nx
import itertools as it
import numpy as np
import igraph as ig
import time 
import math

class GraphClass():
    
    def __init__(self, N, L):
        self.N = N
        self.L = L
        self.attachedfile = None
        self.Gset = None
        self.size = None
        self.Gpoly = None
    
    def attach(self, filename):
        self.attachedfile = filename
    
    
    def load(self, fname = None):
        if not fname:
            if self.attachedfile:
                fname = self.attachedfile
            else:
                raise ValueError('Missing graph set file to load!')
        self.Gset = load_adjacency(fname)
        self.size = len(self.Gset)
        self.Gcoeff = tuple([nRelpoly(g) for g in self.Gset])
        self.Gpoly = tuple([coeff_to_poly(co) for co in self.Gcoeff])
        
    
    def analysis(self, specify = 'all',rec = False, sh_report = False, interact = False, sample = True):
        
        processed_poly = {'extrema': tuple(map(np.polyder,self.Gpoly, it.repeat(1))),
                         'inflection_point':tuple(map(np.polysub,self.Gpoly, it.repeat([1,0]))),
                         'fixed_point':tuple(map(np.polysub,self.Gpoly, it.repeat([1,0]))),
                         'intersection_point': self.Gpoly,
                          'all': ('extrema','inflection_point','fixed_point','intersection_point')
                         }
        
        if specify == 'intersection_point':
            if self.size >= 200:
                if interact:
                    warn = input("""
Warning! There are too many graphs in the graph set, computation could cost very long time. 
Randomly sampling 100 graphs in the graph set.
To override this feature, press 'o' and confirm by pressing enter. Otherwise continue by pressing enter.""")
                    if warn != 'o':
                        processed_poly[specify] = tuple(np.array(processed_poly[specify],dtype=object)[rand.sample(range(0, self.size), 200)])
                        print('Successfully sampled {} graphs from the original graph set.'.format(len(processed_poly[specify])))
                else:
                    if sample:
                        processed_poly[specify] = tuple(np.array(processed_poly[specify],dtype=object)[rand.sample(range(0, self.size), 200)])
            pairs = tuple(combi(processed_poly[specify], 2))
            if interact:
                print('There are in total {} different graph pairs.'.format(len(pairs)))
            processed_poly[specify] = tuple(map(np.polysub,[pair[0] for pair in pairs],[pair[1] for pair in pairs]))
        
        count_record = {}
        value_record = {}
        if specify == 'all':
            for spec in processed_poly[specify]:
                count, values = polyarray_realroot_statistics(processed_poly[spec], rec=False)
                count_record[spec] = count
                value_record[spec] = values
        else:
            count, tot = polyarray_realroot_statistics(processed_poly[spec], rec=False)
            count_record[specify] = count
            value_record[specify] = values
        if sh_report:
            pass
            # report_gen('extrema',count)
        if rec:
            return count_record, value_record
        else:
            return count_record
    
    def plotALL(self, x = np.linspace(0,1,101)):
        plt.figure()
        plt.xlim([0,1])
        plt.ylim([0,1])
        for g in self.Gpoly:
            plt.plot(x,g(x), pltpoly=True)
        plt.show()


def coeff_to_poly(coeff, pltpoly = False, x=np.linspace(0,1,101), plt_arg = None):
    poly1 = np.poly1d(coeff)
    if pltpoly:
        plt.plot(x,poly1(x))
    
    return poly1


def findrealroots(p, interval = (0.0001,0.9999)):
    roots = p.r
    realroots = roots.real[abs(roots.imag)<1e-5]
    return realroots[np.where(np.logical_and(interval[0]<= realroots ,realroots<= interval[1]))]


def polyarray_realroot_statistics(polyarray, rec=False):
    count = {}
    tot = {}

    for i,p in enumerate(polyarray):
        rr = findrealroots(p)
        key = str(len(rr))
        if key in count.keys():
            count[key] += 1
        else:
            count[key] = 1
        if rec:
            tot[str(i)] = rr
    return count, tot



def combi( N:int, num:int)->it.combinations:
    """
        Find all comibinations of 'num' elements for int in [0,N-1].
        
        Input: 
            N     total number of elements (numbering starts from zero!)
            num   number of elements in one combination

        Output:
            C     the list with all combinations
    """
    return it.combinations(tuple(range(N)),num)


def combi_all(N:int) -> tuple:
    """
    Find all comibinations for int in [0,N-1].
    
    Input: 
        N    total number of elements (numbering starts from zero!)

    Output:
        AC   the list with all combinations, arranged from the empty set to N-element sets
    """
    return tuple(map(combi,it.repeat(N),range(N+1)))


def examine_connected_sets(G: ig.Graph,f: int) -> int:
    """
    This function enumerates the connected sets of cardinality k = N-f. f represents the cardinality of the failure node set.
    Input: 
        G    a igraph.Graph object
        f    number of nodes to be removed

    Output:
        c_n_minus_f   the cardinality of the k-ordered connected sets.
    """
    
    c_N_minus_f = 0
    for failed_vertice_set in combi(G.vcount(),f):
        GG = G.copy()
        GG.delete_vertices(failed_vertice_set)
        if GG.is_connected() == True:
            c_N_minus_f += 1
    return c_N_minus_f

def to_ig(H) -> ig.Graph:
    """
    This function converts the input graph H to the correponding igraph Graph object.
    Supported format of H is described below.
    Input:
        H       H can be a list of tuples that describes all edges in a graph, 
                a numpy matrix that represents the adjacency matrix,
                an igraph graph object,
                or a networkx graph object.
    Output:
        G       An igraph graph object that correspond to the input graph H.
    """
    
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


def connected_set_poly(H) -> np.ndarray:
    """
    This function generates the connected set node polynomial of a graph
    Input:
        H       H can be a list of tuples that describes all edges in a graph, 
                a numpy matrix that represents the adjacency matrix,
                an igraph graph object,
                or a networkx graph object.
    Output:
        c_k       The node connected set polynomial coefficients, 
                begins with the leading coefficient(highest degree) and ends with the constant coefficient(zero degree)
    """
    #Processing the input graph, convert input H to an igraph graph
    G = to_ig(H) 
    N = G.vcount()
    
    c_k = np.zeros(N+1) 
    # Initiallize the array of k-cardinality connected set coefficient c_k. 
    # Where the numpy array c_k start with the leading coefficient when k = N. 
    # => c_k[0] is c_N, c_k[N+1] is c_(N-1), ...,
    # c_k[N-1] is c_1, and c_k[N] is c_0 which is always 0. 
    if N>0:
        c_k[N-1] = N 
        if N>1:
            c_k[N-2] = G.ecount()
            kappa_det = G.cohesion()
            # Kappa is the vertex connectivity of a graph, then for the first kappa-number terms of highest degree, the coefficients can be efficiently calculated using the binomial coefficients. 
            c_k[0:kappa_det] = np.array(list(map(math.comb,it.repeat(N),range(N-kappa_det+1,N+1))))[::-1]
            # For the rest of the coefficients, enumerate all possible sets and count the connected sets. 
            c_k[kappa_det:N-2] = np.array(list(map(examine_connected_sets,[G.copy() for _ in range(kappa_det,N-2)],range(kappa_det,N-2))))
    return c_k.astype(int)


def nRelpoly(H):
    """
    Calculates the coefficients of node reliability polynomial for the given graph H.
    Input:
        H    H can be a list of tuples that describes all edges in a graph, 
             or a networkx graph.
    Output:
        coeff    The array of coefficients of the node reliability polynomial for the graph.
    """
    G = to_ig(H)
    N = G.vcount()
    # Get the connected-set polynomial coefficients
    c_k = connected_set_poly(G)[::-1]
    coeff = np.zeros(N+1)
    # Find the coefficient from the definition of node reliability polynomial. 
    # c_k*(p^k)*(1-p)^(N-k)
    for k in range(N+1):
        coeff[k:] += np.array([c_k[k]*((-1)**i)*math.comb(N-k,i) for i in range(N-k+1)])
    
    return coeff[::-1]