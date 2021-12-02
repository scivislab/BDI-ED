import numpy as np
import math

from contourMergeTrees_helpers import *
            
def editDistance_constrained(nodes1,topo1,rootID1,nodes2,topo2,rootID2,editCost,traceback=False):
    memT = dict()
    memF = dict()
    
    #===================================================================
    # Recursive helper function that computes edit distance between two subtrees rooted in curr1,curr2
    def editDistance_constrained_tree(curr1,curr2):
        # If both trees empty, return distance 0
        if(curr1<0 and curr2<0):
            return 0
        # If first tree empty, return deletion cost of curr2 plus deletion cost of its subforest
        if(curr1<0):
            if((curr1,curr2) not in memT):
                d = (editDistance_constrained_forest(curr1,curr2)
                     + editCost(None,nodes2[curr2]))
                memT[(curr1,curr2)] = d
            return memT[(curr1,curr2)]
        # If second tree empty, return deletion cost of curr1 plus deletion cost of its subforest
        if(curr2<0):
            if((curr1,curr2) not in memT):
                d = (editDistance_constrained_forest(curr1,curr2)
                     + editCost(nodes1[curr1],None))
                memT[(curr1,curr2)] = d
            return memT[(curr1,curr2)]
        # If both trees not empty, find optimal edit operation
        if((curr1,curr2) not in memT):
            # Case for matching curr1 to curr2 and recurse to subforests
            d = (editDistance_constrained_forest(curr1,curr2)
                 + editCost(nodes1[curr1],nodes2[curr2]))
            # Cases for deleting curr1 and mapping one child to curr2
            for child1 in topo1[curr1]:
                d_ = (editDistance_constrained_tree(curr1,-1)
                      + editDistance_constrained_tree(child1,curr2)
                      - editDistance_constrained_tree(child1,-1))
                d = min(d,d_)
            # Cases for deleting curr2 and mapping one child to curr1
            for child2 in topo2[curr2]:
                d_ = (editDistance_constrained_tree(-1,curr2)
                      + editDistance_constrained_tree(curr1,child2)
                      - editDistance_constrained_tree(-1,child2))
                d = min(d,d_)
            memT[(curr1,curr2)] = d
        return memT[(curr1,curr2)]

    #===================================================================
    # Recursive helper function that computes edit distance between two subforests rooted in curr1,curr2
    def editDistance_constrained_forest(curr1,curr2):
        # If both trees empty, return distance 0
        if(curr1<0 and curr2<0):
            return 0
        # If first forest empty, return deletion cost of all trees in second forest
        if(curr1<0):
            if((curr1,curr2) not in memF):
                d=0.0
                for child2 in topo2[curr2]:
                    d += editDistance_constrained_tree(curr1,child2)
                memF[(curr1,curr2)] = d
            return memF[(curr1,curr2)]
        # If second forest empty, return deletion cost of all trees in first forest
        if(curr2<0):
            if((curr1,curr2) not in memF):
                d=0.0
                for child1 in topo1[curr1]:
                    d += editDistance_constrained_tree(child1,curr2)
                memF[(curr1,curr2)] = d
            return memF[(curr1,curr2)]
        # If both forests not empty, find optimal edit operations
        if((curr1,curr2) not in memF):
            d = float("inf")
            childMatchings = getAllMatchings(topo1[curr1],topo2[curr2])
            # Cases for mapping some trees in first forest to some trees in second forest and deleting all other trees in both forests
            for m in childMatchings:
                d_ = 0
                # trees not mapped in first forest
                for child1 in m[0]:
                    d_ += editDistance_constrained_tree(child1,-1)
                # trees not mapped in second forest
                for child2 in m[2]:
                    d_ += editDistance_constrained_tree(-1,child2)
                # trees mapped
                for child1,child2 in m[1]:
                    d_ += editDistance_constrained_tree(child1,child2)
                d = min(d,d_)
            # Cases for mapping one subtree in first forest to second forest and deleting all other trees in first forest 
            for child1 in topo1[curr1]:
                if(len(topo1[child1])==0):
                    continue
                d_ = (editDistance_constrained_forest(curr1,-1)
                      + editDistance_constrained_forest(child1,curr2)
                      - editDistance_constrained_forest(child1,-1))
                d = min(d,d_)
            # Cases for mapping one subtree in second forest to first forest and deleting all other trees in second forest 
            for child2 in topo2[curr2]:
                if(len(topo2[child2])==0):
                    continue
                d_ = (editDistance_constrained_forest(-1,curr2)
                      + editDistance_constrained_forest(curr1,child2)
                      - editDistance_constrained_forest(-1,child2))
                d = min(d,d_)
            memF[(curr1,curr2)] = d
        return memF[(curr1,curr2)]
    
    #===================================================================
    # Recursive helper function that computes the optimal edit mapping between two trees rooted in curr1,curr2 given the memoization table from distance computation
    def editDistance_constrained_tree_traceback(curr1,curr2):
        if(curr1<0 and curr2<0):
            return []
        if(curr1<0):
            match = editDistance_constrained_forest_traceback(curr1,curr2)
            match.append((-1,curr2))
            return match
        if(curr2<0):
            match = editDistance_constrained_forest_traceback(curr1,curr2)
            match.append((curr1,-1))
            return match
        d = (editDistance_constrained_forest(curr1,curr2)
             + editCost(nodes1[curr1],nodes2[curr2]))
        if(memT[(curr1,curr2)]==d):
            match = editDistance_constrained_forest_traceback(curr1,curr2)
            match.append((curr1,curr2))
            return match
        for child1 in topo1[curr1]:
            d = (editDistance_constrained_tree(curr1,-1)
                 + editDistance_constrained_tree(child1,curr2)
                 - editDistance_constrained_tree(child1,-1))
            if(memT[(curr1,curr2)]==d):
                match = (editDistance_constrained_tree_traceback(curr1,-1)
                         + editDistance_constrained_tree_traceback(child1,curr2))
                match_ = editDistance_constrained_tree_traceback(child1,-1)
                for m in match_:
                    match.remove(m)
                return match
        for child2 in topo2[curr2]:
            d = (editDistance_constrained_tree(-1,curr2)
                 + editDistance_constrained_tree(curr1,child2)
                 - editDistance_constrained_tree(-1,child2))
            if(memT[(curr1,curr2)]==d):
                match = (editDistance_constrained_tree_traceback(-1,curr2)
                         + editDistance_constrained_tree_traceback(curr1,child2))
                match_ = editDistance_constrained_tree_traceback(-1,child2)
                for m in match_:
                    match.remove(m)
                return match

    #===================================================================
    # Recursive helper function that computes the optimal edit mapping between two subforests rooted in curr1,curr2 given the memoization table from distance computation
    def editDistance_constrained_forest_traceback(curr1,curr2):
        if(curr1<0 and curr2<0):
            return []
        if(curr1<0):
            match = []
            for child2 in topo2[curr2]:
                match += editDistance_constrained_tree_traceback(curr1,child2)
            return match
        if(curr2<0):
            match = []
            for child1 in topo1[curr1]:
                match += editDistance_constrained_tree_traceback(child1,curr2)
            return match
        childMatchings = getAllMatchings(topo1[curr1],topo2[curr2])
        for m in childMatchings:
            d_ = 0
            match = []
            for child1 in m[0]:
                d_ += editDistance_constrained_tree(child1,-1)
                match += editDistance_constrained_tree_traceback(child1,-1)
            for child2 in m[2]:
                d_ += editDistance_constrained_tree(-1,child2)
                match += editDistance_constrained_tree_traceback(-1,child2)
            for child1,child2 in m[1]:
                d_ += editDistance_constrained_tree(child1,child2)
                match += editDistance_constrained_tree_traceback(child1,child2)
            if(memF[(curr1,curr2)] == d_):
                return match
        for child1 in topo1[curr1]:
            if(len(topo1[child1])==0):
                continue
            d_ = (editDistance_constrained_forest(curr1,-1)
                  + editDistance_constrained_forest(child1,curr2)
                  - editDistance_constrained_forest(child1,-1))
            if(memF[(curr1,curr2)] == d_):
                match = (editDistance_constrained_forest_traceback(curr1,-1)
                         + editDistance_constrained_forest_traceback(child1,curr2))
                match_ = editDistance_constrained_forest_traceback(child1,-1)
                for m in match_:
                    match.remove(m)
                return match
        for child2 in topo2[curr2]:
            if(len(topo2[child2])==0):
                continue
            d_ = (editDistance_constrained_forest(-1,curr2)
                  + editDistance_constrained_forest(curr1,child2)
                  - editDistance_constrained_forest(-1,child2))
            if(memF[(curr1,curr2)] == d_):
                match = (editDistance_constrained_forest_traceback(-1,curr2)
                         + editDistance_constrained_forest_traceback(curr1,child2))
                match_ = editDistance_constrained_forest_traceback(-1,child2)
                for m in match_:
                    match.remove(m)
                return match
    
    #===================================================================
    # if traceback flag set, return distance and mapping, otherwise only distance
    if(traceback):
        return editDistance_constrained_tree(rootID1,rootID2),editDistance_constrained_tree_traceback(rootID1,rootID2)
    else:
        return editDistance_constrained_tree(rootID1,rootID2)
