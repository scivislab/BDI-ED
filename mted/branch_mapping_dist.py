import math
import numpy as np
from scipy.optimize import linear_sum_assignment
        
def branchMappingDistance(tree1,rootID1,tree2,rootID2,editCost,sqrt=False,traceback=False):
    memT = dict()

    nodes1 = []
    for i in range(len(list(tree1.nodes))):
        nodes1.append(tree1.nodes(data=True)[i]["scalar"])
    topo1 = []
    for i in range(len(list(tree1.nodes))):
        topo1.append(list(tree1.neighbors(i)))

    nodes2 = []
    for i in range(len(list(tree2.nodes))):
        nodes2.append(tree2.nodes(data=True)[i]["scalar"])
    topo2 = []
    for i in range(len(list(tree2.nodes))):
        topo2.append(list(tree2.neighbors(i)))

    #===================================================================================
    # Recursive helper function that computes edit distance between two subtrees rooted in (parent1,curr1),(parent2,curr2)
    def editDistance_branch(curr1,parent1,curr2,parent2):
        #===============================================================================
        # if both trees are empty, return 0
        if(curr1<0 and curr2<0):
            return 0
        #===============================================================================
        # If first tree empty, delete entire second subtree
        if(curr1<0):
            if((curr1,parent1,curr2,parent2) not in memT):
                #-----------------------------------------------------------------------
                # If second subtree has only one branch, return deletion cost of this branch
                if(len(topo2[curr2])==0):
                    memT[(curr1,parent1,curr2,parent2)] = editCost(None,None,nodes2[curr2],nodes2[parent2])
                #-----------------------------------------------------------------------
                # If second subtree has more than one branch, try all decompositions
                else:
                    c = float("inf")
                    for child2_mb in topo2[curr2]:
                        c_ = editDistance_branch(curr1,parent1,child2_mb,parent2)
                        for child2 in topo2[curr2]:
                            if(child2==child2_mb):
                                continue
                            c_ += editDistance_branch(curr1,parent1,child2,curr2)
                        c = min(c,c_)
                    memT[(curr1,parent1,curr2,parent2)] = c
            return memT[(curr1,parent1,curr2,parent2)]
        #===============================================================================
        # If second tree empty, delete entire first subtree
        if(curr2<0):
            if((curr1,parent1,curr2,parent2) not in memT):
                #-----------------------------------------------------------------------
                # If first subtree has only one branch, return deletion cost of this branch
                if(len(topo1[curr1])==0):
                    memT[(curr1,parent1,curr2,parent2)] = editCost(nodes1[curr1],nodes1[parent1],None,None)
                #-----------------------------------------------------------------------
                # If first subtree has more than one branch, try all decompositions
                else:
                    c = float("inf")
                    for child1_mb in topo1[curr1]:
                        c_ = editDistance_branch(child1_mb,parent1,curr2,parent2)
                        for child1 in topo1[curr1]:
                            if(child1==child1_mb):
                                continue
                            c_ += editDistance_branch(child1,curr1,curr2,parent2)
                        c = min(c,c_)
                    memT[(curr1,parent1,curr2,parent2)] = c
            return memT[(curr1,parent1,curr2,parent2)]
        #===============================================================================
        # If both trees not empty, find optimal edit operation
        if((curr1,parent1,curr2,parent2) not in memT):
            #---------------------------------------------------------------------------
            # If both trees only have one branch, return edit cost between the two branches
            if(len(topo1[curr1])==0 and len(topo2[curr2])==0):
                memT[(curr1,parent1,curr2,parent2)] = editCost(nodes1[curr1],nodes1[parent1],nodes2[curr2],nodes2[parent2])
            #---------------------------------------------------------------------------
            # If first tree only has one branch, try all decompositions of second tree
            elif(len(topo1[curr1])==0):
                d = float("inf")
                for child2_mb in topo2[curr2]:
                    d_ = editDistance_branch(curr1,parent1,child2_mb,parent2)
                    for child2 in topo2[curr2]:
                        if(child2==child2_mb):
                            continue
                        d_ += editDistance_branch(-1,-1,child2,curr2)
                    d = min(d,d_)
                memT[(curr1,parent1,curr2,parent2)] = d
            #---------------------------------------------------------------------------
            # If second tree only has one branch, try all decompositions of first tree
            elif(len(topo2[curr2])==0):
                d = float("inf")
                for child1_mb in topo1[curr1]:
                    d_ = editDistance_branch(child1_mb,parent1,curr2,parent2)
                    for child1 in topo1[curr1]:
                        if(child1==child1_mb):
                            continue
                        d_ += editDistance_branch(child1,curr1,-1,-1)
                    d = min(d,d_)
                memT[(curr1,parent1,curr2,parent2)] = d
            #---------------------------------------------------------------------------
            # If both trees have more than one branch, try all decompositions of both trees
            else:
                d = float("inf")
                #-----------------------------------------------------------------------
                # Try all possible main branches of first tree (child1_mb) and all possible main branches of second tree (child2_mb)
                # Then try all possible matchings of subtrees
                # Special case of binary trees is treated differently for performance
                if(len(topo1[curr1])==2 and len(topo2[curr2])==2):
                    child11 = topo1[curr1][0]
                    child12 = topo1[curr1][1]
                    child21 = topo2[curr2][0]
                    child22 = topo2[curr2][1]
                    d = min(d,editDistance_branch(child11,parent1,child21,parent2) + editDistance_branch(child12,curr1,child22,curr2))
                    d = min(d,editDistance_branch(child12,parent1,child22,parent2) + editDistance_branch(child11,curr1,child21,curr2))
                    d = min(d,editDistance_branch(child11,parent1,child22,parent2) + editDistance_branch(child12,curr1,child21,curr2))
                    d = min(d,editDistance_branch(child12,parent1,child21,parent2) + editDistance_branch(child11,curr1,child22,curr2))
                # For non-binary trees use compute distance through maximum matching
                else:
                    for child1_mb in topo1[curr1]:
                        topo1_ = topo1[curr1].copy()
                        topo1_.remove(child1_mb)
                        for child2_mb in topo2[curr2]:
                            d_ = editDistance_branch(child1_mb,parent1,child2_mb,parent2)
                            topo2_ = topo2[curr2].copy()
                            topo2_.remove(child2_mb)                            
                            deg = max(len(topo1_),len(topo2_))
                            matchMatrix = np.zeros((deg,deg))
                            for i in range(deg):
                                child1 = topo1_[i] if i<len(topo1_) else -1
                                for j in range(deg):
                                    child2 = topo2_[j] if j<len(topo2_) else -1
                                    matchMatrix[i,j] = editDistance_branch(child1,curr1,child2,curr2)
                            row_ind, col_ind = linear_sum_assignment(matchMatrix)
                            d_ += matchMatrix[row_ind, col_ind].sum()
                            d = min(d,d_)
                #-----------------------------------------------------------------------
                # Try to continue main branch on one child of first tree and delete all other subtrees
                # Then match continued branch to current branch in second tree
                for child1_mb in topo1[curr1]:
                    d_ = editDistance_branch(child1_mb,parent1,curr2,parent2)
                    for child1 in topo1[curr1]:
                        if(child1 == child1_mb):
                            continue
                        d_ += editDistance_branch(child1,curr1,-1,-1)
                    d = min(d,d_)
                #-----------------------------------------------------------------------
                # Try to continue main branch on one child of second tree and delete all other subtrees
                # Then match continued branch to current branch in first tree
                for child2_mb in topo2[curr2]:
                    d_ = editDistance_branch(curr1,parent1,child2_mb,parent2)
                    for child2 in topo2[curr2]:
                        if(child2 == child2_mb):
                            continue
                        d_ += editDistance_branch(-1,-1,child2,curr2)
                    d = min(d,d_)
                memT[(curr1,parent1,curr2,parent2)] = d
        return memT[(curr1,parent1,curr2,parent2)]

    #===================================================================================
    # Recursive helper function that computes the optimal edit mapping between two subtrees rooted in (parent1,curr1),(parent2,curr2) given the memoization table from distance computation
    def editDistance_branch_traceback(curr1,parent1,curr2,parent2):
        #===============================================================================
        # base case
        if(curr1<0 and curr2<0):
            return []
        #===============================================================================
        # base case (first tree null)
        if(curr1<0):
            if(len(topo2[curr2])==0):
                return [((-1,-1),(curr2,parent2))]
            else:
                c = memT[(curr1,parent1,curr2,parent2)]
                for child2_mb in topo2[curr2]:
                    c_ = editDistance_branch(curr1,parent1,child2_mb,parent2)
                    for child2 in topo2[curr2]:
                        if(child2==child2_mb):
                            continue
                        c_ += editDistance_branch(curr1,parent1,child2,curr2)
                    if(c==c_):
                        match = editDistance_branch_traceback(curr1,parent1,child2_mb,parent2)
                        for child2 in topo2[curr2]:
                            if(child2==child2_mb):
                                continue
                            match += editDistance_branch_traceback(curr1,parent1,child2,curr2)
                        return match
        #===============================================================================
        # base case (second tree null)
        if(curr2<0):
            if(len(topo1[curr1])==0):
                return [((curr1,parent1),(-1,-1))]
            else:
                c = memT[(curr1,parent1,curr2,parent2)]
                for child1_mb in topo1[curr1]:
                    c_ = editDistance_branch(child1_mb,parent1,curr2,parent2)
                    for child1 in topo1[curr1]:
                        if(child1==child1_mb):
                            continue
                        c_ += editDistance_branch(child1,curr1,curr2,parent2)
                    if(c==c_):
                        match = editDistance_branch_traceback(child1_mb,parent1,curr2,parent2)
                        for child1 in topo1[curr1]:
                            if(child1==child1_mb):
                                continue
                            match += editDistance_branch_traceback(child1,curr1,curr2,parent2)
                        return match
        #===============================================================================
        # both trees not null
        
        #------------------------------------------------
        # both trees leaves
        if(len(topo1[curr1])==0 and len(topo2[curr2])==0):
            #print((curr1,parent1)," ",(curr2,parent2))
            return [((curr1,parent1),(curr2,parent2))]
        #------------------------------------------------
        # first tree leave
        elif(len(topo1[curr1])==0):
            d = memT[(curr1,parent1,curr2,parent2)]
            for child2_mb in topo2[curr2]:
                d_ = editDistance_branch(curr1,parent1,child2_mb,parent2)
                for child2 in topo2[curr2]:
                    if(child2==child2_mb):
                        continue
                    d_ += editDistance_branch(-1,-1,child2,curr2)
                if(d==d_):
                    match = editDistance_branch_traceback(curr1,parent1,child2_mb,parent2)
                    for child2 in topo2[curr2]:
                        if(child2==child2_mb):
                            continue
                        match += editDistance_branch_traceback(-1,-1,child2,curr2)
                    return match
        #------------------------------------------------
        # second tree leave
        elif(len(topo2[curr2])==0):
            d = memT[(curr1,parent1,curr2,parent2)]
            for child1_mb in topo1[curr1]:
                d_ = editDistance_branch(child1_mb,parent1,curr2,parent2)
                for child1 in topo1[curr1]:
                    if(child1==child1_mb):
                        continue
                    d_ += editDistance_branch(child1,curr1,-1,-1)
                if(d==d_):
                    match = editDistance_branch_traceback(child1_mb,parent1,curr2,parent2)
                    for child1 in topo1[curr1]:
                        if(child1==child1_mb):
                            continue
                        match += editDistance_branch_traceback(child1,curr1,-1,-1)
                    return match
        #------------------------------------------------
        # both trees inner nodes
        else:
            d = memT[(curr1,parent1,curr2,parent2)]
            if(len(topo1[curr1])==2 and len(topo2[curr2])==2):
                child11 = topo1[curr1][0]
                child12 = topo1[curr1][1]
                child21 = topo2[curr2][0]
                child22 = topo2[curr2][1]
                if(d == editDistance_branch(child11,parent1,child21,parent2) + editDistance_branch(child12,curr1,child22,curr2)):
                    return editDistance_branch_traceback(child11,parent1,child21,parent2) + editDistance_branch_traceback(child12,curr1,child22,curr2)
                if(d == editDistance_branch(child12,parent1,child22,parent2) + editDistance_branch(child11,curr1,child21,curr2)):
                    return editDistance_branch_traceback(child12,parent1,child22,parent2) + editDistance_branch_traceback(child11,curr1,child21,curr2)
                if(d == editDistance_branch(child11,parent1,child22,parent2) + editDistance_branch(child12,curr1,child21,curr2)):
                    return editDistance_branch_traceback(child11,parent1,child22,parent2) + editDistance_branch_traceback(child12,curr1,child21,curr2)
                if(d == editDistance_branch(child12,parent1,child21,parent2) + editDistance_branch(child11,curr1,child22,curr2)):
                    return editDistance_branch_traceback(child12,parent1,child21,parent2) + editDistance_branch_traceback(child11,curr1,child22,curr2)
            else:
                for child1_mb in topo1[curr1]:
                    topo1_ = topo1[curr1].copy()
                    topo1_.remove(child1_mb)
                    for child2_mb in topo2[curr2]:
                        d_ = editDistance_branch(child1_mb,parent1,child2_mb,parent2)
                        topo2_ = topo2[curr2].copy()
                        topo2_.remove(child2_mb)                            
                        deg = max(len(topo1_),len(topo2_))
                        matchMatrix = np.zeros((deg,deg))
                        for i in range(deg):
                            child1 = topo1_[i] if i<len(topo1_) else -1
                            for j in range(deg):
                                child2 = topo2_[j] if j<len(topo2_) else -1
                                matchMatrix[i,j] = editDistance_branch(child1,curr1,child2,curr2)
                        row_ind, col_ind = linear_sum_assignment(matchMatrix)
                        d_ += matchMatrix[row_ind, col_ind].sum()
                        if(d == d_):
                            match = editDistance_branch_traceback(child1_mb,parent1,child2_mb,parent2)
                            for i in range(len(row_ind)):
                                child1 = topo1_[row_ind[i]] if row_ind[i]<len(topo1_) else -1
                                child2 = topo2_[col_ind[i]] if col_ind[i]<len(topo2_) else -1
                                match += editDistance_branch_traceback(child1,curr1,child2,curr2)
                            return match
            for child1_mb in topo1[curr1]:
                d_ = editDistance_branch(child1_mb,parent1,curr2,parent2)
                for child1 in topo1[curr1]:
                    if(child1 == child1_mb):
                        continue
                    d_ += editDistance_branch(child1,curr1,-1,-1)
                if(d==d_):
                    match_ = editDistance_branch_traceback(child1_mb,parent1,curr2,parent2)
                    for child1 in topo1[curr1]:
                        if(child1 == child1_mb):
                            continue
                        match_ += editDistance_branch_traceback(child1,curr1,-1,-1)
                    return match_
            for child2_mb in topo2[curr2]:
                d_ = editDistance_branch(curr1,parent1,child2_mb,parent2)
                for child2 in topo2[curr2]:
                    if(child2 == child2_mb):
                        continue
                    d_ += editDistance_branch(-1,-1,child2,curr2)
                if(d==d_):
                    match_ = editDistance_branch_traceback(curr1,parent1,child2_mb,parent2)
                    for child2 in topo2[curr2]:
                        if(child2 == child2_mb):
                            continue
                        match_ += editDistance_branch_traceback(-1,-1,child2,curr2)
                    return match_
    
    #===============================================================================
    # if traceback flag set, return distance and mapping, otherwise only distance
    if(traceback):
        if sqrt:
            return math.sqrt(editDistance_branch(topo1[rootID1][0],rootID1,topo2[rootID2][0],rootID2)), editDistance_branch_traceback(topo1[rootID1][0],rootID1,topo2[rootID2][0],rootID2)
        else:
            return editDistance_branch(topo1[rootID1][0],rootID1,topo2[rootID2][0],rootID2), editDistance_branch_traceback(topo1[rootID1][0],rootID1,topo2[rootID2][0],rootID2)
    else:
        if sqrt:
            return math.sqrt(editDistance_branch(topo1[rootID1][0],rootID1,topo2[rootID2][0],rootID2))
        else:
            return editDistance_branch(topo1[rootID1][0],rootID1,topo2[rootID2][0],rootID2)
