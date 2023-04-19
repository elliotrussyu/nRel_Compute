import re
import os


def checkNL(N,L):
    """
    Checks whether the given graph order N and size L pair is valid.
    """
    return N-1 <= L <= N*(N-1) and 1 <= N <= 63


def checkfilename(filename):
    """
    Checks whether the user specified file is valid or not.
    Criteria: 
    1. txt extension or no extension name
    2. must contain {N}n{L}l to indicate the number of nodes and links of the graph set contained 
    """
    filename = filename.lower()
    #matching graph name
    graph = re.search('\d+n+\d+l',filename)
    #Checking file extension
    if re.search('\.+\S',filename):
        if filename.split('.')[-1] != 'txt':
            raise ValueError('Invalid file type, this program only supports .txt file or no extension')
    #If there is a valid match {N}n{L}l
    if graph: 
        N = re.search('\d+n',graph[0])[0][:-1]
        L = re.search('\d+l',graph[0])[0][:-1]
        if checkNL(int(N),int(L)): # If the node and link number is valid
            return int(N),int(L)
        else:
            raise ValueError('Invalid graph set! Check the order and size of the graph.')
    
    raise ValueError('Invalid filename. The standard filename for N-node-L-link graph set should be formatted as  {N}n{L}l')
    
    
def graph6_to_adjacency(g,N,L = '0'):
    """
    The main method of this module. This method converts a graph6 format graph to its adjacency matrix.
    Input:
        -g    graph6 format graph
        -N    graph order
        -L    graph size (only for verification purpose). When set to '0'(default), the checking step is ignored.
    Output:
        -A    The converted adjacency matrix.
    """
    upper_triangle_size = int(N*(N-1)/2)
    Aseed = [ord(i) for i in g]
    A_upper = []
    for i in Aseed[1:]:
        A_upper.extend(bin(i-63)[2:].zfill(6)) 
        #The true A_upper is in fact the first N*(N-1)/2 elements, zeros are padded at the end in graph6 encoding
    
    if type(L) == 'int' and sum(A_upper) != L:
        raise ValueError('Graph size mismatch error! The actual size of the decoded graph is {} however it was stated the graph should be of size {}. Check whether the filename matches the actual graph set stored'.format(sum(A_upper), L))
    
    A = [ [0]*N for _ in range(N) ]
    c=0
    for i in range(N):
        for j in range(i):
            A[i][j] = int(A_upper[c])
            A[j][i] = int(A_upper[c])
            c += 1
    return A

    
def read_and_convert(filename,flag=0):
    """
    This method reads the file given by the user. Checks the validity of the file, and converts all graphs inside to corresponding adjacency matrices.
    Input:
        -filename    The file contains the graph6 formatted graphs
        -flag        This flag is for printing, when set to 0(default), the program prints nothing.
    Output:
        -A_tot    A python list contains all of the adjacency matrices.
        -N        graph order
        -L        graph size
    """
    N, L = checkfilename(filename)
    with open(filename,'r') as g:
        if flag:
            print('***Start converting to adjacency matrices...')
        lines = g.readlines()
        lines = [x[:-1] for x in lines if x!='\n']
        A_tot = []
        for g in lines:
            A_tot.append(graph6_to_adjacency(g, N, L))
        if flag:
            print('***Convertion complete...')
            print('***Successfully converted {} {}-nodes-{}-links graphs'.format(len(lines),N,L))
    
    return A_tot, N, L


def write_A_to_file(A_tot, N, L, flag=0):
    """
    Writes the adjacency matrices to a text file named {order}_nodes_{size}_links_{number of graphs}_graphs.txt
    Input:
        -A_tot    The list of all adjacency matrices
        -N        graph order
        -L        graph size
        -flag     This flag is for printing, when set to 0(default), the program prints nothing.
    
    The file is stored in a new folder named: GraphSetsAdj
    """
    if flag:
        print('***Writing to file...')
    fAdj = 'GraphSetsAdj\{}_nodes_{}_links_{}_graphs.txt'.format(N,L,len(A_tot))
    with open(fAdj,'w') as f:
        f.write('\n')
        cc = 0
        for A in A_tot:
            cc += 1
            f.write('\n')
            f.write('Graph'+str(cc)+', order '+str(N)+', size '+str(L)+'\n')
            for r in A:
                f.write('\n')
                f.write(' '.join(map(repr, r)))
                f.write('\n')
    if flag:
        print('***Write to file successful! All {} non-isomorphic {}-nodes-{}-links graphs can be found in the text file: {}'.format(len(A_tot),N,L,fAdj))    

        
        
        
if __name__ == '__main__':
    import argparse
    
    ProgInitTXT = "This program transfers Graph6 format graph sets with certain number of nodes and links generated with 'geng' in nauty into corresponding adjacency matrices. The graph number starts with 1.\n|***NOTE: this program only supports the graph file:\n|***1)without header,\n|***2)using Graph6 format,\n|***3)currently only supports graph order smaller than 63,\n|***4)only supports plaintext files without file extension or with .txt extension name,\n|***5)each line in the file contains only one graph."
    
    parser = argparse.ArgumentParser(description=ProgInitTXT)
    parser.add_argument('-f',"--file", help="Specifies the file contains graph6 format graphs",default='')
    parser.add_argument("-i", "--interactive", help="Run in interactive mode.", action="store_true")
    args = parser.parse_args()
    
    print_flag = 1
    recurr_flag = 0
    
    if args.interactive or args.file == '':
        print('-- No valid filename detected, automatically entering the interactive mode.--\n')
        print(ProgInitTXT)
    
    while 1:
        if args.interactive or args.file == '':
            recurr_flag = 1       
            fname = input('\nPlease provide the file contains the graph set below. (Type exit to exit the program)\n File:  ')
            if fname == 'exit':
                break
        else:
            fname = args.file
            
        ### Main function ###
        try:
            A_all,N,L = read_and_convert(fname, print_flag)
            write_A_to_file(A_all,N,L, print_flag)
        except ValueError:
            print('Invalid file provided, try again!')
            
        if recurr_flag == 0:
            break