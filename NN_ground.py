import numpy as np
import random
import tfim
import itertools as it
import argparse
import networkx as nx
from line_profiler import LineProfiler
import json

import sys, os
from itertools import combinations


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('yheight', type=int, help='Height of grid')
    parser.add_argument('xwidth', type=int, help='Width of grid')
    parser.add_argument('initial_seed', type=int, help='First Jij seed')
    parser.add_argument('seed_range', type=int, help='Number of seeds')
    args = parser.parse_args()

    PBC = True
    yheight = args.yheight
    xwidth = args.xwidth
    L = [yheight, xwidth]
    lattice = tfim.Lattice(L, PBC)
    N = lattice.N
    initial = args.initial_seed
    num_seeds = args.seed_range
    center = (500,375)
    ground_states = {}
    
    f = open("ground_states.txt", "w+")

    for seed in range(initial, initial + num_seeds):

        bonds = bond_list(seed, N, PBC, xwidth, yheight)

        Jij = make_Jij(N, bonds, lattice)
        
        coordList = spinCoords(center, xwidth, yheight)
        
        plaq = make_plaquettes(PBC, lattice, N, xwidth, yheight)
        
        f_plaq = frustrated(Jij, plaq)
        
        node_pairs = plaq_pairing(f_plaq, coordList, PBC, xwidth, yheight)
        
        init_ground = initial_ground(node_pairs, xwidth, yheight)
        
        p_pairings = init_ground[0]
        
        ground_distance = init_ground[1]
        
        edges = viable_edges(node_pairs, p_pairings, ground_distance, f_plaq, xwidth, yheight)
        
        matchings = plaq_groups(edges, f_plaq)
        
        string_groups = add_all_strings(matchings, lattice, coordList)

        b_bonds = broken_bonds(string_groups, N, coordList, xwidth, yheight)
        
        true_ground = make_config(b_bonds, Jij, N, xwidth, lattice, string_groups)
        
        ground_config = true_ground[0]
        
        true_ground_strings = true_ground[1]
        
        number_ground_states = len(true_ground_strings)
        
        if number_ground_states == 0:
            found_ground = False
        else:
            found_ground = True
        
        incremented = 0
        while found_ground == False:
            incremented += 1
            ground_distance += 1
            edges = viable_edges(node_pairs, p_pairings, ground_distance, f_plaq, xwidth, yheight)
            
            if len(edges) != 0:
                
                matchings = plaq_groups(edges, f_plaq)
                
                string_groups = add_all_strings(matchings, lattice, coordList)

                b_bonds = broken_bonds(string_groups, N, coordList, xwidth, yheight)
                
                true_ground = make_config(b_bonds, Jij, N, xwidth, lattice, string_groups)
                
                ground_config = true_ground[0]
                
                true_ground_strings = true_ground[1]
                
                number_ground_states = len(true_ground_strings)
                
                if number_ground_states != 0:
                    found_ground = True
        ground_states.update({seed:ground_config})
        print('ground config: ', ground_config)
        print('seed: ', seed)
        print('incremented: ', incremented)
    f.write(json.dumps(ground_states))
    f.close()



def bond_list(seed, N, PBC, xwidth, yheight):
    np.random.seed(seed)
    # Generates a random list of bonds with equal numbers of ferromagnetic and antiferromagnetic bonds
    if PBC == True:
        num_of_bonds = 2*N
    else:
        num_of_bonds = (xwidth - 1)*(yheight) + (xwidth)*(yheight - 1)
    if num_of_bonds%2 == 0:
        a1 = [-1 for i in range(num_of_bonds//2)]
    else:
        a1 = [-1 for i in range((num_of_bonds//2) + 1)]
    a2 = [1 for i in range(num_of_bonds//2)]
    a = list(np.random.permutation(a1+a2))
    return a

def make_Jij(N, b_list, lattice):
    bond_index = 0
    Jij = np.zeros((N,N))
    for i in range(0,N):
        NNs = lattice.NN(i)
        for j in NNs:
            if Jij[i][j] == 0:
                Jij[i][j] = b_list[bond_index]
                Jij[j][i] = b_list[bond_index]
                bond_index += 1
    return Jij
    
    
def make_plaquettes(PBC, lattice, N, xwidth, yheight):
    p_list = []
    if PBC:
        for i in range(0, N):
            NNs = lattice.NN(i)
            plaq = [i]
            plaq.append(NNs[3])
            NNs2 = lattice.NN(NNs[3])
            plaq.append(NNs2[1])
            plaq.append(NNs[1])
            p_list.append(plaq)
    else:
        for y in range(0,yheight):
            for x in range(0, xwidth):
                if y == yheight-1 or x == xwidth-1: #This part adds empty plaquettes so the first number is also the index of each plaquette
                    #if i want to take it out, just need to subtract 1 from x and y range
                    p_list.append([])
                else:
                    plaq = []
                    i = y*xwidth + x
                    plaq.append(i)
                    plaq.append(i+1)
                    plaq.append(i+xwidth+1)
                    plaq.append(i+xwidth)
                    p_list.append(plaq)
    return p_list
    
    
def frustrated(Jij, plaqs):
    f_plaq = []
    for plaq in plaqs:
        count = 0
        if len(plaq)!=0:
            if Jij[plaq[0]][plaq[1]] == -1:
                count += 1
            if Jij[plaq[1]][plaq[2]] == -1:
                count += 1
            if Jij[plaq[2]][plaq[3]] == -1:
                count += 1
            if Jij[plaq[0]][plaq[3]] == -1:
                count += 1
            if count == 1 or count == 3:
                f_plaq.append(plaq[0])
    return f_plaq
        
        
def plaq_pairing(f_plaq, coordList, PBC, xwidth, yheight):
    '''Function returns a list of all possible pairs between frustrated plaquettes with the distances between them–The distance is stored as the maximum possible distance between two plaquettes minus the actual distance'''
    pair_list = []
    for index, p1 in enumerate(f_plaq):
        coord1 = coordList[p1]
        for p2 in f_plaq[index+1:]:
            coord2 = coordList[p2]
            x1 = coord1[0]
            x2 = coord2[0]
            
            y1 = coord1[1]
            y2 = coord2[1]
            
            xdiff = abs((x1 - x2))
            ydiff = abs((y2-y1))
            if PBC:
                if xdiff > (xwidth)//2:
                    xdiff = (xwidth) - xdiff
                if ydiff > (yheight)//2:
                    ydiff = (yheight) - ydiff
            else:
                if xdiff > (xwidth-1)//2:
                    xdiff = (xwidth - 1) - xdiff
                if ydiff > (yheight-1)//2:
                    ydiff = (yheight) - ydiff
            tot_dist = int(xdiff + ydiff)
            #Here we build a list of pairs with the distance between them
            max = xwidth//2 + yheight//2
            op_dist = (max - tot_dist)
            if p1 > p2:
                pair_list.append((p2, p1, op_dist))
            else:
                pair_list.append((p1, p2, op_dist)) #Pair_list does not have true distance between pairs, the true distance is subtracted from the max possible distance to help with node matching later
    return pair_list



def initial_ground(pair_list, xwidth, yheight):

  G = nx.Graph()
  G.add_weighted_edges_from(pair_list) #makes graph of all node pairs
  matching = nx.max_weight_matching(G) #gives one solution
  ground_dist = 0
  p_pairs = []
  for pair in matching:
      edge = G.get_edge_data(pair[0], pair[1])
      pair_dist = (xwidth//2 + yheight//2)-edge['weight']
      ground_dist += pair_dist #total string length for a ground state soln
      if pair[0] > pair[1]:
          p0 = pair[1]
          p1 = pair[0]
          pair = (p0, p1)
      p_pairs.append([pair, pair_dist]) #adds solution from above to list with pairs and pair distance
  return p_pairs, ground_dist



def viable_edges(pair_list, p_pairs, ground_dist, f_plaq, xwidth, yheight):
  '''Function takes the list of all possible pairings of nodes and returns a list of all combinations of those nodes that would result in a ground state energy'''
  
  edge_lst = []
  plaq_dict = {}
  
  #Make the list for the edges grouped by plaquette
  for index, plaq in enumerate(f_plaq):
      edge_lst.append([])
      plaq_dict[plaq] = index
  
  G = nx.Graph()
  G.add_weighted_edges_from(pair_list)
 
  first = False
  p_dist = 0
  for pair in p_pairs:
      dist = pair[1]
      p_dist += dist
  if p_dist == ground_dist:
      first = True #So we only remove certain edges if p_pairs is relevant
  
  if first:
      for pair in p_pairs:
          plaq = pair[0][0]
          ind = plaq_dict.get(plaq)
          edge_lst[ind].append(pair)
      
  for plaq in f_plaq:
      G2 = G.copy()
      if first:
          for pair in p_pairs:
              if pair[0][1] == plaq or pair[0][0] == plaq:
                  G2.remove_edge(*pair[0])
                  break
              
      ground_energy = True
      while ground_energy == True:
          matching = nx.max_weight_matching(G2)
          if len(matching) != len(f_plaq)/2: #This would happen if we have taken out all edges for a particular plaquette
              ground_energy = False
              break
          new_length = 0
          new_group = []
          for pair in matching:
              edge = G2.get_edge_data(pair[0], pair[1])
              if pair[0] == plaq or pair[1] == plaq:
                  rem_edge = (pair[0], pair[1])
              pair_dist = (xwidth//2 + yheight//2)-edge['weight']
              new_length += pair_dist
              if pair[0] > pair[1]:
                  p0 = pair[1]
                  p1 = pair[0]
                  pair = (p0, p1)
              new_group.append([pair, pair_dist])
          if new_length == ground_dist:
              G2.remove_edge(*rem_edge)
              for pair in new_group:
                  plaq = pair [0][0]
                  ind = plaq_dict.get(plaq)
                  
                  if pair not in edge_lst[ind]:
                      edge_lst[ind].append(pair)
          elif new_length < ground_dist:
              G2.remove_edge(*rem_edge)
          else:
              ground_energy = False
              ind = plaq_dict[plaq]
  zeroes = True
  for plaq in edge_lst:
      if len(plaq) != 0:
          zeroes = False
          break
  if zeroes:
      edge_lst = []
  return edge_lst


def plaq_groups(edges, f_plaq):
    
    group = []
    used_plaquettes = []
    all_groups = []
    current_plaq = 0
    p_ind = 0
    index = 0
    loop_count = 0
    
    new = False
    running = True
    
    plaq_dict = {}
    for index, plaq in enumerate(f_plaq): #Allows me to find the index of the frustrated plaquettes
        plaq_dict[plaq] = index
    
    #The main piece of the function
    while running:
        for group_index, p_edges in enumerate(edges[current_plaq:]):
            if new:
                new = False #new allows me to restart the for loop when i change current_plaq
                break
                
            if group_index + current_plaq == len(edges) - 1: #if we get to the last group without having a full ground state, we need to use a different combination of edges
                try_new = False
                for_loop = False
                for e_ind, edge in enumerate(group[::-1]): #going through the ground state group backwards to see if other edge choices will work
                    loop_count +=1
                    if loop_count > 1000000:
                        running = False
                        new = True
                    for_loop = True
                    if try_new == True:
                        break
                    else:
                        plaq_ind = plaq_dict.get(edge[0][0])
                        for e_index, edge2 in enumerate(edges[plaq_ind]):
                            if edge2 == edge:
                                if e_index == len(edges[plaq_ind])-1 and plaq_ind == 0: #end of program, we've reached the last entry of the first plaquette
                                    running = False
                                    new = True
                                    try_new = True
                                    break
                                elif e_index == len(edges[plaq_ind])-1: #Move to the previous plaquette list to find a viable edge
                                    break
                                else:
                                    current_plaq = plaq_ind #move to the next edge in the list for the plaquette
                                    p_ind = e_index + 1
                                    try_new = True
                                    new = True
                                    break
                if for_loop:
                    group = group[:-e_ind]
                    used_plaquettes = used_plaquettes[:-e_ind]
                    if len(group) == 1 and group[0][0][0] == plaq_dict.get(current_plaq):
                        group = []
                        used_plaquettes = []
            for pair in p_edges[p_ind:]: #This is a single edge from edge list
                p_ind = 0 #Resets p_ind for the next loop through
                if pair[0][0] in used_plaquettes:
                    break #Can move to next plaquette
                elif pair[0][1] in used_plaquettes:
                    continue #Need to go through to the next pair
                else:
                    group.append(pair)
                    used_plaquettes.append(pair[0][1])#maybe add the first element of pair if this is the first plaquette we visit
                    if len(group) == len(f_plaq)//2: #Group is full
                    
                        all_groups.append(group)
                        last_pair = group[-2] #This is the pair that we remove and replace before cycling through other options
                        ind = plaq_dict.get(last_pair[0][0]) #The plaquette index
                        group = group[:-2]
                        used_plaquettes = used_plaquettes[:-2]
                        found = False
                        while found == False:
                            loop_count += 1
                            if loop_count > 1000000:
                                running = False
                                new = True
                                found = True
                                break
                            for index, pairing in enumerate(edges[ind]):
                                if pairing == last_pair and index == len(edges[ind])-1: #This happens if we are at the last pair of a particular plaquette
                                    if len(group) == 0: #This happens if we have gotten through the last edge of the first plaquette, function is done
                                        running = False
                                        found = True
                                        break
                                    last_pair = group[-1] #Take off the last pair and go to that plaquette to see if there are further pairs to use
                                    ind = plaq_dict.get(last_pair[0][0])
                                    group = group[:-1]
                                    used_plaquettes = used_plaquettes[:-1]
                                elif pairing == last_pair: #This means there are more pairs for the plaquette in question, so we adjust current_plaq and p_ind, and go through the for loops again from there
                                    current_plaq = ind
                                    p_ind = index + 1
                                    found = True
                                    new = True
                                    break
                    break
    return all_groups
    

def add_all_strings(groups, lattice, coordList):
    edges = []
    for i in range(len(coordList)):
        NNs = lattice.NN(i)
        for j in NNs:
            if i < j:
                edges.append((i,j))
    G = nx.Graph()
    G.add_edges_from(edges) #G has edges connecting all points in a lattice with PBC
    
    all_groups = []
    index = 0
    for group in groups:
    
        if index > 20000:
            print('Not all ground states found')
            break
    
        single_pairing = []
        index += 1
        for pairing in group:
            paths = nx.all_shortest_paths(G, pairing[0][0], pairing[0][1]) #finds all possible paths between two points
            paths_list = []
            for path in paths:
                paths_list.append(path)
            single_pairing.append(paths_list)
        path_combos = it.product(*single_pairing)
        for combo in path_combos:
            all_groups.append(combo)
    return all_groups
    

def broken_bonds(string_groups, N, coordList, xwidth, yheight):
    '''Returns an NxN matrix with 1's where there are broken bonds between two spins'''
    config_Jij_list = []
    for str_index, state in enumerate(string_groups):
        config_Jij = np.zeros((N,N))
        for string in state:
            for index in range(0, len(string)-1):
                p1 = string[index]
                p2 = string[index + 1]
                if p1>p2:
                    hold = p1
                    p1 = p2
                    p2 = hold
                c1x = coordList[p1][0]
                c2x = coordList[p2][0]
                if c1x == c2x:
                    if p2 + xwidth > N - 1 and p1 < xwidth:
                        sp1 = p1
                        if (p1+1) % xwidth == 0:
                            sp2 = p1 - xwidth + 1
                        else:
                            sp2 = p1 + 1
                    else:
                        sp1 = p2
                        if (p1+1) % xwidth == 0: #on the far right
                            sp2 = p2 - xwidth + 1
                        else:
                            sp2 = p2 + 1
                else:
                    if p2 + xwidth > N - 1:
                        if (p2+1) % xwidth == 0:
                            if p1 % xwidth == 0:
                                sp1 = p1
                            else:
                                sp1 = p2
                        sp2 = sp1 - (xwidth * (yheight - 1))
                    elif (p2+1) % xwidth == 0:
                        if p1 % xwidth == 0:
                            sp1 = p1
                        else:
                            sp1 = p2
                        sp2 = sp1 + xwidth
                    else:
                        sp1 = p2
                        sp2 = p2 + xwidth
                bond = (sp1, sp2)
                config_Jij[sp1][sp2] = 1
                config_Jij[sp2][sp1] = 1
        config_Jij_list.append([config_Jij, str_index])
    return config_Jij_list
    
def make_config(b_bonds, bonds, N, xwidth, lattice, string_groups):
    ground_states = []
    true_strings = []
    for Jij in b_bonds:
        broken = Jij[0]
        spin_list = []
        spin_list.append(0) #Set the first spin as down
        valid = True
        
        #Loop through all other spins
        for sp1 in range(1, N):
            if valid == False:
                break
            if sp1 % xwidth == 0:
                sp2 = sp1 - xwidth
            else:
                sp2 = sp1 - 1
            
            spin2 = spin_list[sp2]
            bond = bonds[sp1][sp2]
            status = broken[sp1][sp2]
            
            #Set spin
            if bond == 1:
                #Spins want to be the same
                if status == 1: #broken
                    spin1 = abs(spin2 - 1)
                else:
                    spin1 = spin2
            else:
                #Spins want to be opposite
                if status == 1: #broken
                    spin1 = spin2
                else:
                    spin1 = abs(spin2 - 1)
            spin_list.append(spin1)
            
            #Check bonds to all lower spins
            NNs = lattice.NN(sp1)
            for i in NNs:
                if i < sp1:
                    spini = spin_list[i]
                    bond = bonds[sp1][i]
                    status = broken[sp1][i]
                    
                    if bond == 1: #Spins want to be same
                        if status == 1: #Spins should be opposite
                            if spin1 == spini:
                                valid = False
                                break
                        else: #Spins should be same
                            if spin1 != spini:
                                valid = False
                                break
                    else: #Spins want to be opposite
                        if status == 1: #spins should be same
                            if spini != spin1:
                                valid = False
                                break
                        else: #Spins should be opposite
                            if spini == spin1:
                                valid = False
                                break
        if valid:
            index = 0
            for i in range(0, N):
                if spin_list[i] == 1:
                    index += 2**i
            spin_list.reverse()
            ground_states.append(index)
            true_strings.append(string_groups[Jij[1]])
            
    
    return ground_states, true_strings
    
    
def spinCoords(center, xwidth, yheight):
    coords = []
    y_init = center[1] - ((yheight)/2)
    for j in range(0, yheight):
        y = y_init + (j)
        x_init = center[0] - ((xwidth)/2)
        for i in range(0, xwidth):
            x = x_init + i
            c = (x,y)
            coords.append(c)
    return coords


if __name__ == "__main__":
    main()
    
