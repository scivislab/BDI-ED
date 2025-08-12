import networkx as nx
import time
import pulp as pl
import time

def rooted_unordered_deformation_distance_with_solution(tree1, root1, tree2, root2, isSubProblem=False,forced_pair=None, hard_time_limit=None,initial_time_limit=None,known_upper_bound=None):

    down_min_cost_fun_with_IP_buf = {}
    up_min_cost_fun_with_IP_buf   = {}
    
    down_min_cost_fun_with_IP_buf_to_cache = {}
    up_min_cost_fun_with_IP_buf_to_cache   = {}
    iteration = 0

    def rooted_unordered_deformation_distance_with_limits(tree1, root1, tree2, root2, upper_bound, time_limit, solution=None, use_lower_bound=False, isSubProblem=False, forced_pair=None, hard_time_limit=None):
        labels1    = [0 for x in range(tree1.number_of_nodes())]
        for i in tree1.nodes:
            labels1[i] = tree1.nodes(data=True)[i]["scalar"]

        labels2 = [0 for x in range(tree2.number_of_nodes())]
        for i in tree2.nodes:
            labels2[i] = tree2.nodes(data=True)[i]["scalar"]
            
        oriented_tree1 = nx.dfs_tree(tree1, root1)
        oriented_tree2 = nx.dfs_tree(tree2, root2)

        # print(labels1)
        # print(tree1.number_of_nodes(),root1,oriented_tree1.number_of_nodes())
        # print("")
        # print(labels2)
        # print(tree2.number_of_nodes(),root2,oriented_tree2.number_of_nodes())
        # print("")
        
        def is_ancestor1(n1, n2):
            d = nx.descendants(oriented_tree1, n1)
            return (n2 in d)
        def is_ancestor_or_equal1(n1, n2):
            return (is_ancestor1(n1,n2) or n1==n2)
        # print(nx.descendants(oriented_tree1, root1))
            
        def is_ancestor2(n1, n2):
            d = nx.descendants(oriented_tree2, n1)
            return (n2 in d)
        def is_ancestor_or_equal2(n1, n2):
            return (is_ancestor2(n1,n2) or n1==n2)
        # print(nx.descendants(oriented_tree2, root2))
        
        # print("")

        # predecessors[i][j] contains all vertices between i and j
        # ToDo: should be optimized!
        predecessors1 = [[[] for x in range(tree1.number_of_nodes())] for y in range(tree1.number_of_nodes())]
        predecessors2 = [[[] for x in range(tree2.number_of_nodes())] for y in range(tree2.number_of_nodes())]
        depth1 = 0
        depth2 = 0
        for i in tree1.nodes:
            for j in nx.descendants(oriented_tree1,i):
                curr = list(oriented_tree1.in_edges(j))[0][0]
                while curr != i:
                    predecessors1[i][j].append(curr)
                    curr = list(oriented_tree1.in_edges(curr))[0][0]
                depth1 = max(depth1,len(predecessors1[i][j]))
        for i in tree2.nodes:
            for j in nx.descendants(oriented_tree2,i):
                curr = list(oriented_tree2.in_edges(j))[0][0]
                while curr != i:
                    predecessors2[i][j].append(curr)
                    curr = list(oriented_tree2.in_edges(curr))[0][0]
                depth2 = max(depth2,len(predecessors2[i][j]))
        # print(predecessors1)
        # print(predecessors2)

        subtree_del_costs1 = [0] * tree1.number_of_nodes()
        subtree_del_costs_wo_parent1 = [0] * tree1.number_of_nodes()
        def fill_subtree_del_costs1(i):
            subtree_size = 0
            for c in oriented_tree1.neighbors(i):
                fill_subtree_del_costs1(c)
                subtree_del_costs1[c] += abs(labels1[c] - labels1[i])
                subtree_size += subtree_del_costs1[c]
            subtree_del_costs1[i] = subtree_size
            subtree_del_costs_wo_parent1[i] = subtree_size

        subtree_del_costs2 = [0] * tree2.number_of_nodes()
        subtree_del_costs_wo_parent2 = [0] * tree2.number_of_nodes()
        def fill_subtree_del_costs2(i):
            subtree_size = 0
            for c in oriented_tree2.neighbors(i):
                fill_subtree_del_costs2(c)
                subtree_del_costs2[c] += abs(labels2[c] - labels2[i])
                subtree_size += subtree_del_costs2[c]
            subtree_del_costs2[i] = subtree_size
            subtree_del_costs_wo_parent2[i] = subtree_size

        fill_subtree_del_costs1(root1)
        fill_subtree_del_costs2(root2)

        def path_cost1(n1,p1):
            last = n1
            curr = list(oriented_tree1.in_edges(n1))[0][0]
            size = 0
            # print("---------------")
            # print(predecessors1[p1][n1])
            for i in range(len(predecessors1[p1][n1])):
                curr = predecessors1[p1][n1][i]
                for c in oriented_tree1.neighbors(curr):
                    if c!=last:
                        assert(subtree_del_costs1[c] > 0)
                        size += subtree_del_costs1[c]
                last = curr
            # print(size)
            # print("---------------")
            return size

        def path_cost2(n2,p2):
            last = n2
            curr = list(oriented_tree2.in_edges(n2))[0][0]
            size = 0
            for i in range(len(predecessors2[p2][n2])):
                curr = predecessors2[p2][n2][i]
                for c in oriented_tree2.neighbors(curr):
                    if c!=last:
                        assert(subtree_del_costs2[c] > 0)
                        size += subtree_del_costs2[c]
                last = curr
            return size

        furthest_descendant1 = [0] * tree1.number_of_nodes()
        furthest_descendant2 = [0] * tree2.number_of_nodes()

        for i in tree1.nodes:
            for j in nx.descendants(oriented_tree1,i):
                furthest_descendant1[i] = max(furthest_descendant1[i],abs(labels1[i]-labels1[j]))
        for i in tree2.nodes:
            for j in nx.descendants(oriented_tree2,i):
                furthest_descendant2[i] = max(furthest_descendant2[i],abs(labels2[i]-labels2[j]))
            
        def d_relabel(n1,p1, n2,p2):
            ln1 = labels1[n1]
            lp1 = labels1[p1]
            ln2 = labels2[n2]
            lp2 = labels2[p2]
            return abs(abs(ln1-lp1) - abs(ln2-lp2))
        
        def d_delete(n2,p2):
            # print(n2,p2)
            return abs(labels2[n2]-labels2[p2])
            
        def d_insert(n1,p1):
            # print(n1,p1)
            return abs(labels1[n1]-labels1[p1])
            
        # def d_relabel(n1, n2):
        #     l1 = labels1[n1]
        #     l2 = labels2[n2]
        #     if(tree1_type[n1] == tree2_type[n2]):
        #         return max(l1 - l2, l2 - l1)
        #     else:
        #         return (l1 + l2)
        
        # def d_delete(n2):
        #     return labels2[n2]
            
        # def d_insert(n1):
        #     return labels1[n1]
        
        #mapping variables
        m_vars = [[0 for x in range(tree2.number_of_nodes())] for y in range(tree1.number_of_nodes())] # is i mapped to j
        pm_vars  = [[[[0 for a in range(tree2.number_of_nodes())] # is path (a-b) mapped to path (c-d)
                      for b in range(tree2.number_of_nodes())]
                      for c in range(tree1.number_of_nodes())]
                      for d in range(tree1.number_of_nodes())]
        all_tuples = []
        all_path_tuples = []
        all_path_tuples1 = []
        all_path_tuples2 = []

        # structural information about paths in T1
        d_vars1 = [0 for x in range(tree1.number_of_nodes())] # is i deleted
        dd_vars1 = [0 for x in range(tree1.number_of_nodes())] # are all descendants of i and i itself deleted
        two_children_not_dd_vars1 = [0 for x in range(tree1.number_of_nodes())] # does i have more than one child remaining
        one_child_not_dd_vars1 = [0 for x in range(tree1.number_of_nodes())] # does i have one child remaining 
        pruned_vars1 = [0 for x in range(tree1.number_of_nodes())] # is i pruned
        parent_vars1 = [[0 for x in range(tree1.number_of_nodes())] for y in range(tree1.number_of_nodes())] # is i the imaginary parent of j

        # structural information about paths in T2
        d_vars2 = [0 for x in range(tree2.number_of_nodes())] # is j deleted
        dd_vars2 = [0 for x in range(tree2.number_of_nodes())] # are all descendants of j and j itself deleted
        two_children_not_dd_vars2 = [0 for x in range(tree2.number_of_nodes())] # does i have more than one child remaining
        one_child_not_dd_vars2 = [0 for x in range(tree2.number_of_nodes())] # does i have one child remaining 
        pruned_vars2 = [0 for x in range(tree2.number_of_nodes())] # is j pruned
        parent_vars2 = [[0 for x in range(tree2.number_of_nodes())] for y in range(tree2.number_of_nodes())] # is i the imaginary parent of j
        
        for i in tree1.nodes:
            for j in tree2.nodes:
                all_tuples += [[i, j]]
                m_vars[i][j] = pl.LpVariable("m_"+str(i)+"_"+str(j), 0, 1, pl.LpInteger)
        for i in tree1.nodes:
            d_vars1[i] = pl.LpVariable("d_"+str(i)+"_t1", 0, 1, pl.LpInteger)
            dd_vars1[i] = pl.LpVariable("dd_"+str(i)+"_t1", 0, 1, pl.LpInteger)
            two_children_not_dd_vars1[i] = pl.LpVariable("two_children_not_dd_"+str(i)+"_t1", 0, 1, pl.LpInteger)
            one_child_not_dd_vars1[i] = pl.LpVariable("one_child_not_dd_"+str(i)+"_t1", 0, 1, pl.LpInteger)
            pruned_vars1[i] = pl.LpVariable("pruned_"+str(i)+"_t1", 0, 1, pl.LpInteger)
        for j in tree2.nodes:
            d_vars2[j] = pl.LpVariable("d_"+str(j)+"_t2", 0, 1, pl.LpInteger)
            dd_vars2[j] = pl.LpVariable("dd_"+str(j)+"_t2", 0, 1, pl.LpInteger)
            two_children_not_dd_vars2[j] = pl.LpVariable("two_children_not_dd_"+str(j)+"_t2", 0, 1, pl.LpInteger)
            one_child_not_dd_vars2[j] = pl.LpVariable("one_child_not_dd_"+str(j)+"_t2", 0, 1, pl.LpInteger)
            pruned_vars2[j] = pl.LpVariable("pruned_"+str(j)+"_t2", 0, 1, pl.LpInteger)
        for i in tree1.nodes:
            for j in tree1.nodes:
                if is_ancestor1(i,j):
                    parent_vars1[i][j] = pl.LpVariable("parent_"+str(i)+"_"+str(j)+"_t1", 0, 1, pl.LpInteger)
                    all_path_tuples1.append([j,i])
        for i in tree2.nodes:
            for j in tree2.nodes:
                if is_ancestor2(i,j):
                    parent_vars2[i][j] = pl.LpVariable("parent_"+str(i)+"_"+str(j)+"_t2", 0, 1, pl.LpInteger)
                    all_path_tuples2.append([j,i])
        # for i in tree1.nodes:
        #     for j in tree1.nodes:
        #         for k in tree2.nodes:
        #             for l in tree2.nodes:
        #                 all_path_tuples += [[i,j,k,l]]
        
        down_min_cost_fun_with_IP_buf_to = {}
        up_min_cost_fun_with_IP_buf_to = {}
        
        ip_time_limit = max(time_limit, 80)

        def up_min_cost_fun(n1,p1,n2,p2):
            up_min_cost = abs(abs(labels1[root1]-labels1[p1]) - abs(labels2[root2]-labels2[p2]))
            if len(predecessors1[p1][n1]) == 0:
                first_on_path1 = n1
            else:
                first_on_path1 = predecessors1[p1][n1][-1]
            if len(predecessors2[p2][n2]) == 0:
                first_on_path2 = n2
            else:
                first_on_path2 = predecessors2[p2][n2][-1]
            up_total1 = subtree_del_costs1[root1] - subtree_del_costs1[first_on_path1]
            up_total2 = subtree_del_costs2[root2] - subtree_del_costs2[first_on_path2]
            up_min_cost = max(up_min_cost,abs(up_total1-up_total2))
            return up_min_cost
        
        def up_min_cost_fun_with_IP(n1,p1,n2,p2):
            # return up_min_cost_fun(n1,p1,n2,p2)
            nonlocal ip_time_limit
            
            if len(predecessors1[p1][n1]) == 0:
                first_on_path1 = n1
            else:
                first_on_path1 = predecessors1[p1][n1][-1]
            if len(predecessors2[p2][n2]) == 0:
                first_on_path2 = n2
            else:
                first_on_path2 = predecessors2[p2][n2][-1]
                
            if (first_on_path1, first_on_path2) in up_min_cost_fun_with_IP_buf.keys() and (first_on_path1, first_on_path2) not in up_min_cost_fun_with_IP_buf_to_cache.keys():
                return up_min_cost_fun_with_IP_buf[(first_on_path1,first_on_path2)]
            if (first_on_path1, first_on_path2) in up_min_cost_fun_with_IP_buf.keys() and (first_on_path1, first_on_path2) in up_min_cost_fun_with_IP_buf_to_cache.keys() and iteration < 4:
                return up_min_cost_fun_with_IP_buf[(first_on_path1,first_on_path2)]
            
            if p1 == root1 or p2 == root2:
                return up_min_cost_fun(n1,p1,n2,p2)
                
            subtree1_non_nodes = list(nx.descendants(oriented_tree1,first_on_path1)) + [first_on_path1]
            subtree2_non_nodes = list(nx.descendants(oriented_tree2,first_on_path2)) + [first_on_path2]
            
            subtree1_nodes = list(set(list(tree1.nodes)) - set(subtree1_non_nodes))
            subtree2_nodes = list(set(list(tree2.nodes)) - set(subtree2_non_nodes))
            
            assert root1 in subtree1_nodes
            assert root2 in subtree2_nodes
            
            assert p1 in subtree1_nodes
            assert p2 in subtree2_nodes
            
            assert n1 not in subtree1_nodes
            assert n2 not in subtree2_nodes
            
            sz1 = len(subtree1_nodes)
            sz2 = len(subtree2_nodes)
            up_min_cost = up_min_cost_fun(n1,p1,n2,p2)
            if(sz1 <= 22 and sz2 <= 22 and sz1 > 1 and sz2 > 1  and (iteration == 2 or iteration >= 3) and ip_time_limit > 0):
                if not isSubProblem:
                    print("IP_up_min_cost", first_on_path1, "-", first_on_path2, " ", ip_time_limit)
                subtree_cutout1  = tree1.subgraph(subtree1_nodes).copy()
                subtree_cutout2  = tree2.subgraph(subtree2_nodes).copy()
                
                subtree1_map = {subtree1_nodes[i]: i for i in range(len(subtree1_nodes))}
                subtree2_map = {subtree2_nodes[i]: i for i in range(len(subtree2_nodes))}
                
                subtree_cutout1 = nx.relabel_nodes(subtree_cutout1, subtree1_map)
                subtree_cutout2 = nx.relabel_nodes(subtree_cutout2, subtree2_map)
                
                assert tree1.nodes(data=True)[root1]["scalar"] == subtree_cutout1.nodes(data=True)[subtree1_map[root1]]["scalar"]  
                assert tree2.nodes(data=True)[root2]["scalar"] == subtree_cutout2.nodes(data=True)[subtree2_map[root2]]["scalar"]  
                
                assert tree1.nodes(data=True)[p1]["scalar"] == subtree_cutout1.nodes(data=True)[subtree1_map[p1]]["scalar"]  
                assert tree2.nodes(data=True)[p2]["scalar"] == subtree_cutout2.nodes(data=True)[subtree2_map[p2]]["scalar"]  
                
                assert subtree1_map[root1] in list(subtree_cutout1.nodes)
                assert subtree2_map[root2] in list(subtree_cutout2.nodes)
                
                assert subtree1_map[p1] in list(subtree_cutout1.nodes)
                assert subtree2_map[p2] in list(subtree_cutout2.nodes)
                
                #lower_bound = rooted_unordered_deformation_distance(subtree_cutout1, subtree1_map[root1], subtree_cutout2, subtree2_map[root2], isSubProblem=True, forced_pair=(subtree1_map[p1], subtree2_map[p2]), hard_time_limit = min(time_limit/4,22.1),initial_time_limit=2.55)
                #if iteration == 2:
                htl = 5.2
                if (first_on_path1, first_on_path2) in up_min_cost_fun_with_IP_buf_to_cache.keys():
                    htl = 120
                    up_min_cost_fun_with_IP_buf_to_cache.pop((first_on_path1, first_on_path2), None)
                #else:
                #    htl = 10.5
                start = time.time()
                lower_bound = rooted_unordered_deformation_distance(subtree_cutout1, subtree1_map[root1], subtree_cutout2, subtree2_map[root2], isSubProblem=True, forced_pair=(subtree1_map[p1], subtree2_map[p2]), hard_time_limit = htl,initial_time_limit=2.55)
                end = time.time()
                ip_time_limit -= (end-start)
                #if lower_bound == None  and iteration <= 2:
                #    print("time out ", ip_time_limit)
                #    up_min_cost_fun_with_IP_buf_to[(first_on_path1,first_on_path2)] = up_min_cost
                #    return up_min_cost
                if lower_bound == None:
                    print("time out (buffered) ", ip_time_limit)
                    lower_bound = up_min_cost
                    if htl < 20:
                        up_min_cost_fun_with_IP_buf_to_cache[(first_on_path1,first_on_path2)] = htl
                #print(lower_bound, "/", up_min_cost)
                up_min_cost = max(lower_bound, up_min_cost)
                up_min_cost_fun_with_IP_buf[(first_on_path1,first_on_path2)] = up_min_cost
            return up_min_cost
        
        def down_min_cost_fun(n1,p1,n2,p2):
            if (n1, n2) in down_min_cost_fun_with_IP_buf.keys():
                return down_min_cost_fun_with_IP_buf[(n1,n2)]
            down_min_cost = abs(furthest_descendant1[n1]-furthest_descendant2[n2])
            full_subtree_size1 = 0
            for c in oriented_tree1.neighbors(n1):
                full_subtree_size1 += subtree_del_costs1[c]
            full_subtree_size2 = 0
            for c in oriented_tree2.neighbors(n2):
                full_subtree_size2 += subtree_del_costs2[c]
            down_min_cost = max(down_min_cost,abs(full_subtree_size1-full_subtree_size2))
            return down_min_cost
            
            
        
        def down_min_cost_fun_with_IP(n1,p1,n2,p2):
            # return down_min_cost_fun(n1,p1,n2,p2)
            nonlocal ip_time_limit
            if (n1, n2) in down_min_cost_fun_with_IP_buf.keys():
                return down_min_cost_fun_with_IP_buf[(n1,n2)]
            if (n1, n2) in down_min_cost_fun_with_IP_buf_to.keys():
                return down_min_cost_fun_with_IP_buf_to[(n1,n2)]
            sz1 = len(nx.descendants(oriented_tree1,n1))
            sz2 = len(nx.descendants(oriented_tree2,n2))
            down_min_cost = down_min_cost_fun(n1,p1,n2,p2)
            if(sz1 <= 30 and sz2 <= 30 and sz1 > 1 and sz2 > 1 and (iteration >= 1 or iteration >= 3) and ip_time_limit > 0):
                subtree1_nodes = list(nx.descendants(oriented_tree1,n1)) + [n1]
                subtree2_nodes = list(nx.descendants(oriented_tree2,n2)) + [n2]
                
                if not isSubProblem:
                    print("IP_down_min_cost", n1, "-", n2, " ", ip_time_limit)
                subtree_cutout1  = tree1.subgraph(subtree1_nodes).copy()
                subtree_cutout2  = tree2.subgraph(subtree2_nodes).copy()
                
                subtree1_map = {subtree1_nodes[i]: i for i in range(len(subtree1_nodes))}
                subtree2_map = {subtree2_nodes[i]: i for i in range(len(subtree2_nodes))}
                
                subtree_cutout1 = nx.relabel_nodes(subtree_cutout1, subtree1_map, copy=True)
                subtree_cutout2 = nx.relabel_nodes(subtree_cutout2, subtree2_map, copy=True)           
                
                #lower_bound = rooted_unordered_deformation_distance(subtree_cutout1, subtree1_map[n1], subtree_cutout2, subtree2_map[n2], isSubProblem=True, hard_time_limit = min(time_limit/4,22.1),initial_time_limit=2.55)
                #if iteration == 2:
                htl = 5.2
                #else:
                #    htl = 10.5
                start = time.time()
                lower_bound = rooted_unordered_deformation_distance(subtree_cutout1, subtree1_map[n1], subtree_cutout2, subtree2_map[n2], isSubProblem=True, hard_time_limit = htl,initial_time_limit=2.55)
                end = time.time()
                ip_time_limit -= (end-start)
                #if lower_bound == None and iteration <= 2:
                #    print("time out ", ip_time_limit)
                #    down_min_cost_fun_with_IP_buf_to[(n1,n2)] = down_min_cost
                #    return down_min_cost
                if lower_bound == None:
                    print("time out (buffered) ", ip_time_limit)
                    lower_bound = down_min_cost
                down_min_cost = lower_bound
                down_min_cost_fun_with_IP_buf[(n1,n2)] = down_min_cost
            return down_min_cost
        
        def down_min_cost_fun_with_mapping_cases(n1,p1,n2,p2):
            return down_min_cost_fun(n1,p1,n2,p2)
            if (n1, n2) in down_min_cost_fun_with_IP_buf.keys():
                return down_min_cost_fun_with_IP_buf[(n1,n2)]
            children1 = list(oriented_tree1.neighbors(n1))
            children2 = list(oriented_tree2.neighbors(n2))
            if not (len(children1)==2 and len(children2)==2):
                return down_min_cost_fun(n1,p1,n2,p2)
            c11 = children1[0]
            c12 = children1[1]
            c21 = children2[0]
            c22 = children2[1]
            min_map_cost_c11toc21 = abs(subtree_del_costs1[c11] - subtree_del_costs2[c21])
            min_map_cost_c12toc22 = abs(subtree_del_costs1[c12] - subtree_del_costs2[c22])
            min_map_cost_c11toc22 = abs(subtree_del_costs1[c11] - subtree_del_costs2[c22])
            min_map_cost_c12toc21 = abs(subtree_del_costs1[c12] - subtree_del_costs2[c21])
            min_map_cost_c11tof2 = abs(subtree_del_costs1[c11] - (subtree_del_costs2[c21]+subtree_del_costs2[c22]))
            min_map_cost_c12tof2 = abs(subtree_del_costs1[c12] - (subtree_del_costs2[c21]+subtree_del_costs2[c22]))
            min_map_cost_f1toc21 = abs((subtree_del_costs1[c11]+subtree_del_costs1[c12]) - subtree_del_costs2[c21])
            min_map_cost_f1toc22 = abs((subtree_del_costs1[c11]+subtree_del_costs1[c12]) - subtree_del_costs2[c22])

            down_min_map_cost_c11toc21_c12toc22 = min_map_cost_c11toc21 + min_map_cost_c12toc22
            down_min_map_cost_c11toc22_c12toc21 = min_map_cost_c11toc22 + min_map_cost_c12toc21
            down_min_map_cost_c11tof2 = min_map_cost_c11tof2 + subtree_del_costs1[c12]
            down_min_map_cost_c12tof2 = min_map_cost_c12tof2 + subtree_del_costs1[c11]
            down_min_map_cost_f1toc21 = min_map_cost_f1toc21 + subtree_del_costs2[c22]
            down_min_map_cost_f1toc22 = min_map_cost_f1toc22 + subtree_del_costs2[c21]

            down_min_map_cost = min([down_min_map_cost_c11toc21_c12toc22,
                                     down_min_map_cost_c11toc22_c12toc21,
                                     down_min_map_cost_c11tof2,
                                     down_min_map_cost_c12tof2,
                                     down_min_map_cost_f1toc21,
                                     down_min_map_cost_f1toc22])

            down_min_furthestpath_cost_simple = abs(furthest_descendant1[n1]-furthest_descendant2[n2])
            down_min_cost_simple = down_min_cost_fun(n1,p1,n2,p2)
            #print('')
            #print(down_min_map_cost_c11toc21_c12toc22,
            #      down_min_map_cost_c11toc22_c12toc21,
            #      down_min_map_cost_c11tof2,
            #      down_min_map_cost_c12tof2,
            #      down_min_map_cost_f1toc21,
            #      down_min_map_cost_f1toc22)
            #print(min_map_cost_c11toc21,
            #    min_map_cost_c12toc22,
            #    min_map_cost_c11toc22,
            #    min_map_cost_c12toc21)
            #print(down_min_cost_simple,down_min_map_cost,down_min_furthestpath_cost_simple)
            #print('')
            down_min_cost = max(down_min_map_cost,down_min_furthestpath_cost_simple)
            return down_min_cost
        
        all_root_path_tuples = []
        for [n1,p1] in all_path_tuples1:
            for [n2,p2] in all_path_tuples2:
                if p1!=root1 or p2!=root2:
                    continue

                down_min_cost = down_min_cost_fun(n1,p1,n2,p2)

                if d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + down_min_cost <= upper_bound:
                    all_root_path_tuples.append([n1,p1,n2,p2])
        
        longest_root_path1 = []
        longest_root_path2 = []
        root_min_relabel = float('inf')
        for [n1,p1,n2,p2] in all_root_path_tuples:
            if len(predecessors1[p1][n1])>len(longest_root_path1):
                longest_root_path1 = predecessors1[p1][n1]
            if len(predecessors2[p2][n2])>len(longest_root_path2):
                longest_root_path2 = predecessors2[p2][n2]
            root_min_relabel = min(root_min_relabel,d_relabel(n1,root1,n2,root2))
        
        shortest_root_path1 = []
        shortest_root_path2 = []
        for v1 in longest_root_path1:
            v1_in_all = True
            for [n1,p1,n2,p2] in all_root_path_tuples:
                if not v1 in (predecessors1[root1][n1]+[n1]):
                    v1_in_all = False
            if v1_in_all:
                shortest_root_path1.append(v1)
        for v2 in longest_root_path2:
            v2_in_all = True
            for [n1,p1,n2,p2] in all_root_path_tuples:
                if not v2 in (predecessors2[root2][n2]+[n2]):
                    v2_in_all = False
            if v2_in_all:
                shortest_root_path2.append(v2)

        root_min_cost = root_min_relabel
        if not isSubProblem:
            print(root_min_cost, len(shortest_root_path1), len(shortest_root_path2))
        if len(shortest_root_path1)>0 and len(shortest_root_path2)>0:
            n1 = shortest_root_path1[0]
            n2 = shortest_root_path2[0]
            if not isSubProblem:
                print(path_cost1(n1,root1), len(shortest_root_path1))
                print(path_cost2(n2,root2), len(shortest_root_path2))
            root_min_cost += path_cost1(n1,root1) + path_cost2(n2,root2)
        if not isSubProblem:
            print(root_min_cost)

        full_size = 0
        simple_opt = 0
        simple_opt2 = 0
        simple_opt3 = 0
        simple_opt4 = 0
        simple_opt5 = 0
        simple_opt6 = 0
        opt = 0
        tuple_to_cost = dict()
        # restrict path mapping variables to those pairs of paths that are possible under the given upper bound
        for [n1,p1] in all_path_tuples1:
            for [n2,p2] in all_path_tuples2:
                full_size += 1
                # if the cost of relabeling the paths is larger than upper bound, mapping them is impossible
                if d_relabel(n1,p1,n2,p2) <= upper_bound:
                    simple_opt += 1

                down_min_cost = down_min_cost_fun(n1,p1,n2,p2)
                down_min_cost_cases = down_min_cost_fun_with_mapping_cases(n1,p1,n2,p2)
                # print(down_min_cost,down_min_cost_cases)
                # down_min_cost_cases = down_min_cost_fun(n1,p1,n2,p2)
                up_min_cost = up_min_cost_fun(n1,p1,n2,p2)
                # up_min_cost = max(up_min_cost,root_min_cost)
                # print(up_min_cost,root_min_cost)

                # if n1-p1 is mapped to n2-p2, then all subtrees branching from nodes between n1/n2 and p1/p1 have to be deleted
                if d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) <= upper_bound:
                    simple_opt2 += 1

                # n1-p1 mapped to n2-p2 implies mapping of subtrees rooted in n1/n2
                if d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + down_min_cost <= upper_bound:
                    simple_opt3 += 1
                
                # n1-p1 mapped to n2-p2 implies mapping of whole trees "above" p1/p2
                if d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + down_min_cost + up_min_cost <= upper_bound:
                    simple_opt4 += 1
                cost_lb = d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + down_min_cost_cases + up_min_cost
                if cost_lb <= upper_bound:
                    simple_opt5 += 1
                
                if cost_lb < upper_bound and time_limit > 20.1: # + upper_bound/4 # cost_lb > (upper_bound/5) and 
                    cost_lb = max(cost_lb, d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2)  + up_min_cost + down_min_cost_fun_with_IP(n1,p1,n2,p2))
                    if cost_lb < upper_bound:
                        cost_lb = max(cost_lb, d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2)  + up_min_cost_fun_with_IP(n1,p1,n2,p2) + down_min_cost_fun_with_IP(n1,p1,n2,p2))
                    
                if cost_lb <= upper_bound:
                    simple_opt6 += 1
                    if isSubProblem or (not ((p1==root1 and p2!=root2) or (p1!=root1 and p2==root2))):
                        opt += 1
                        isLeaf_n1 = len(list(oriented_tree1.neighbors(n1)))==0
                        isLeaf_n2 = len(list(oriented_tree2.neighbors(n2)))==0
                        if (isLeaf_n1==isLeaf_n2):
                            all_path_tuples.append([n1,p1,n2,p2])
                    

                tuple_to_cost[(n1,p1,n2,p2)] = d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + down_min_cost_cases + up_min_cost
                # print(d_relabel(n1,p1,n2,p2), path_cost1(n1,p1), path_cost2(n2,p2))
                
        

        if not isSubProblem:
            print("")
            print("Total persistences:",subtree_del_costs1[root1],subtree_del_costs2[root2])
            print("Upper bound:", upper_bound)
            
            print("")
            print("Tuple optimizations:")
            print(full_size,simple_opt,simple_opt2,simple_opt3,simple_opt4,simple_opt5,simple_opt6,opt,len(all_path_tuples))

        # only do look-ahead optimization if not in a subproblem IP
        if not isSubProblem:
            all_path_tuples1_filtered = []
            all_path_tuples2_filtered = []
            for [n1,p1,n2,p2] in all_path_tuples:
                all_path_tuples1_filtered.append([n1,p1])
                all_path_tuples2_filtered.append([n2,p2])
            all_path_tuples1_filtered = list(set(tuple(i) for i in all_path_tuples1_filtered))
            all_path_tuples2_filtered = list(set(tuple(i) for i in all_path_tuples2_filtered))

            all_path_tuples_filtered = []

            counter = 0
            limit = 2500
            if not isSubProblem:
                print(len(all_path_tuples1_filtered),len(all_path_tuples2_filtered))

            # next we exclude mappings that are impossible for any other mapped pair
            for [n1,p1,n2,p2] in all_path_tuples:

                # insdel_cost = d_delete(p2,root2) + d_insert(p1,root1)

                # only do this optimization for a few tuples
                if counter>limit:
                    all_path_tuples_filtered.append([n1,p1,n2,p2])
                    continue
                counter += 1

                if p1==root1 or p2==root2:
                    # check if n1-p1 and m1-q1 are compatible, i.e. if they can appear in the same mapping
                    candidates1 = []
                    for (m1,q1) in all_path_tuples1_filtered:
                        if not is_ancestor_or_equal1(n1,q1): 
                            continue
                        candidates1.append([m1,q1])

                    # check if n2-p2 and m2-q2 are compatible, i.e. if they can appear in the same mapping
                    candidates2 = []
                    for (m2,q2) in all_path_tuples2_filtered:
                        if not is_ancestor_or_equal2(n2,q2):
                            continue
                        candidates2.append([m2,q2])
                    
                    if len(predecessors1[p1][n1]) == 0:
                        first_on_path1 = n1
                    else:
                        first_on_path1 = predecessors1[p1][n1][-1]
                    if len(predecessors2[p2][n2]) == 0:
                        first_on_path2 = n2
                    else:
                        first_on_path2 = predecessors2[p2][n2][-1]
                    assert(first_on_path1 in oriented_tree1.neighbors(p1))
                    assert(first_on_path2 in oriented_tree2.neighbors(p2))

                    full_subtree_size1 = 0
                    for c in oriented_tree1.neighbors(n1):
                        full_subtree_size1 += subtree_del_costs1[c]
                    full_subtree_size2 = 0
                    for c in oriented_tree2.neighbors(n2):
                        full_subtree_size2 += subtree_del_costs2[c]
                    up_total1 = subtree_del_costs1[root1] - subtree_del_costs1[first_on_path1]
                    up_total2 = subtree_del_costs2[root2] - subtree_del_costs2[first_on_path2]
                    lower_bound_cost = d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + up_total1 + up_total2 + full_subtree_size1 + full_subtree_size2
                    
                    for [m1,q1] in candidates1:
                        if lower_bound_cost <= upper_bound:
                            break
                        for [m2,q2] in candidates2:
                            if lower_bound_cost <= upper_bound:
                                break

                            assert(is_ancestor1(p1,n1))
                            assert(is_ancestor2(p2,n2))
                            assert(is_ancestor1(q1,m1))
                            assert(is_ancestor2(q2,m2))

                            assert(is_ancestor_or_equal1(n1,q1))
                            assert(is_ancestor_or_equal2(n2,q2))

                            down_min_cost = down_min_cost_fun(n1,p1,n2,p2)

                            up_min_cost = 0

                            if len(predecessors1[q1][m1]) == 0:
                                first_on_qm_path1 = m1
                            else:
                                first_on_qm_path1 = predecessors1[q1][m1][-1]
                            if len(predecessors2[q2][m2]) == 0:
                                first_on_qm_path2 = m2
                            else:
                                first_on_qm_path2 = predecessors2[q2][m2][-1]
                            assert(first_on_qm_path1 in oriented_tree1.neighbors(q1))
                            assert(first_on_qm_path2 in oriented_tree2.neighbors(q2))

                            between_size1_alt = subtree_del_costs_wo_parent1[n1] - subtree_del_costs1[first_on_qm_path1]
                            between_size1_alt -= d_insert(q1,n1)

                            between_size2_alt = subtree_del_costs_wo_parent2[n2] - subtree_del_costs2[first_on_qm_path2]
                            between_size2_alt -= d_delete(q2,n2)

                            # between_min_cost = max(abs(between_size1_alt - between_size2_alt),abs(d_delete(p2,m2)-d_insert(p1,m1)))
                            between_min_cost = abs(between_size1_alt - between_size2_alt)+abs(d_delete(q2,n2)-d_insert(q1,n1))

                            lower_bound_cost = min(lower_bound_cost,d_relabel(n1,p1,n2,p2)+ d_relabel(m1,q1,m2,q2)
                                                                    + path_cost1(n1,p1) + path_cost2(n2,p2) + path_cost1(m1,q1) + path_cost2(m2,q2)
                                                                    + down_min_cost + up_min_cost + between_min_cost)

                    if lower_bound_cost <= upper_bound:
                        all_path_tuples_filtered.append([n1,p1,n2,p2])
                    
                    # all_path_tuples_filtered.append([n1,p1,n2,p2])
                    # continue

                else:
                    # check if n1-p1 and m1-q1 are compatible, i.e. if they can appear in the same mapping
                    candidates1 = []
                    for (m1,q1) in all_path_tuples1_filtered:
                        if not is_ancestor_or_equal1(m1,p1): 
                            continue
                        candidates1.append([m1,q1])

                    # check if n2-p2 and m2-q2 are compatible, i.e. if they can appear in the same mapping
                    candidates2 = []
                    for (m2,q2) in all_path_tuples2_filtered:
                        if not is_ancestor_or_equal2(m2,p2):
                            continue
                        candidates2.append([m2,q2])

                    if len(candidates1)==0 or len(candidates2)==0:
                        print("no candidates")
                        continue

                    if len(predecessors1[p1][n1]) == 0:
                        first_on_path1 = n1
                    else:
                        first_on_path1 = predecessors1[p1][n1][-1]
                    if len(predecessors2[p2][n2]) == 0:
                        first_on_path2 = n2
                    else:
                        first_on_path2 = predecessors2[p2][n2][-1]
                    assert(first_on_path1 in oriented_tree1.neighbors(p1))
                    assert(first_on_path2 in oriented_tree2.neighbors(p2))

                    full_subtree_size1 = 0
                    for c in oriented_tree1.neighbors(n1):
                        full_subtree_size1 += subtree_del_costs1[c]
                    full_subtree_size2 = 0
                    for c in oriented_tree2.neighbors(n2):
                        full_subtree_size2 += subtree_del_costs2[c]
                    up_total1 = subtree_del_costs1[root1] - subtree_del_costs1[first_on_path1]
                    up_total2 = subtree_del_costs2[root2] - subtree_del_costs2[first_on_path2]
                    lower_bound_cost = d_relabel(n1,p1,n2,p2) + path_cost1(n1,p1) + path_cost2(n2,p2) + up_total1 + up_total2 + full_subtree_size1 + full_subtree_size2

                    # assert(lower_bound_cost <= upper_bound)
                    for [m1,q1] in candidates1:
                        if lower_bound_cost <= upper_bound:
                            break
                        for [m2,q2] in candidates2:
                            if lower_bound_cost <= upper_bound:
                                break

                            # if n1-n2 and m1-m2 break ancestor-preservation, final mapping cannot contain both tuples
                            if is_ancestor1(n1, m1) != is_ancestor2(n2, m2):
                                assert(False)
                            if is_ancestor1(m1, n1) != is_ancestor2(m2, n2):
                                assert(False)
                            if is_ancestor1(p1,q1) and is_ancestor2(q2,p2):
                                assert(False)
                            if is_ancestor1(q1,p1) and is_ancestor2(p2,q2):
                                assert(False)
                            assert(is_ancestor1(p1,n1))
                            assert(is_ancestor2(p2,n2))
                            assert(is_ancestor1(q1,m1))
                            assert(is_ancestor2(q2,m2))
                            assert(not (is_ancestor1(n1,m1) and is_ancestor2(m2,n2)))
                            assert(not (is_ancestor1(m1,n1) and is_ancestor2(n2,m2)))
                            assert(not (is_ancestor1(p1,q1) and is_ancestor2(q2,p2)))
                            assert(not (is_ancestor1(q1,p1) and is_ancestor2(p2,q2)))
                            assert(m1!=n1)
                            assert(is_ancestor1(p1,n1))
                            assert(is_ancestor_or_equal1(m1,p1))
                            assert(is_ancestor1(q1,m1))
                            assert(m2!=n2)
                            assert(is_ancestor2(p2,n2))
                            assert(is_ancestor_or_equal2(m2,p2))
                            assert(is_ancestor2(q2,m2))

                            # compute implied costs of subtrees rooted in n1,m1,n2,m2
                            down_min_cost = down_min_cost_fun(n1,p1,n2,p2)

                            up_min_cost = 0
                            up_size1 = subtree_del_costs1[root1]
                            up_size2 = subtree_del_costs2[root2]
                            if is_ancestor1(n1,m1):
                                assert(False)
                            elif is_ancestor1(m1,n1):
                                assert(is_ancestor_or_equal2(m2,n2))
                                assert(is_ancestor_or_equal1(m1,p1))
                                assert(is_ancestor_or_equal2(m2,p2))
                                for c in oriented_tree1.neighbors(q1):
                                    if c==m1 or (c in predecessors1[q1][m1]):
                                        up_size1 -= subtree_del_costs1[c]
                                for c in oriented_tree2.neighbors(q2):
                                    if c==m2 or (c in predecessors2[q2][m2]):
                                        up_size2 -= subtree_del_costs2[c]
                            else:
                                assert(False)

                            up_min_cost = abs(up_size1-up_size2)

                            between_min_cost = 0
                            
                            if is_ancestor1(n1,m1):
                                assert(False)
                            
                            assert(is_ancestor1(m1,n1))
                            assert(is_ancestor_or_equal1(m1,p1))
                            assert(is_ancestor_or_equal2(m2,p2))

                            between_size1_alt = subtree_del_costs_wo_parent1[m1] - subtree_del_costs1[first_on_path1]
                            between_size1_alt -= d_insert(p1,m1)

                            between_size2_alt = subtree_del_costs_wo_parent2[m2] - subtree_del_costs2[first_on_path2]
                            between_size2_alt -= d_delete(p2,m2)

                            # between_min_cost = max(abs(between_size1_alt - between_size2_alt),abs(d_delete(p2,m2)-d_insert(p1,m1)))
                            between_min_cost = abs(between_size1_alt - between_size2_alt)+abs(d_delete(p2,m2)-d_insert(p1,m1))

                            # assert(abs(between_size1-between_size1_alt)<0.0001),str(between_size1)+" "+str(between_size1_alt)
                            # assert(abs(between_size2-between_size2_alt)<0.0001),str(between_size2)+" "+str(between_size2_alt)

                            lower_bound_cost = min(lower_bound_cost,d_relabel(n1,p1,n2,p2)+ d_relabel(m1,q1,m2,q2)
                                                                    + path_cost1(n1,p1) + path_cost2(n2,p2) + path_cost1(m1,q1) + path_cost2(m2,q2)
                                                                    + down_min_cost + up_min_cost + between_min_cost)

                    if lower_bound_cost <= upper_bound:
                        all_path_tuples_filtered.append([n1,p1,n2,p2])

            n_root_paths = 0
            for [n1,p1,n2,p2] in all_path_tuples:
                if p1==root1 or p2==root2:
                    n_root_paths += 1
            if not isSubProblem:
                print(full_size,simple_opt,simple_opt2,simple_opt3,simple_opt4,simple_opt5,opt,len(all_path_tuples),len(all_path_tuples_filtered))
            all_path_tuples = all_path_tuples_filtered
        
        #print(counter)
        #print(n_root_paths)
        #print("")

        # all_path_tuples1_filtered = []
        # all_path_tuples2_filtered = []
        # for [n1,p1,n2,p2] in all_path_tuples:
        #     all_path_tuples1_filtered.append([n1,p1])
        #     all_path_tuples2_filtered.append([n2,p2])
        # all_path_tuples1_filtered = list(set(tuple(i) for i in all_path_tuples1_filtered))
        # all_path_tuples2_filtered = list(set(tuple(i) for i in all_path_tuples2_filtered))

        # for [j,i] in all_path_tuples1_filtered:
        #     assert(is_ancestor1(i,j))
        #     parent_vars1[i][j] = pl.LpVariable("parent_"+str(i)+"_"+str(j)+"_t1", 0, 1, pl.LpInteger)
        #     # all_path_tuples1.append([j,i])
        # for [j,i] in all_path_tuples2_filtered:
        #     assert(is_ancestor2(i,j))
        #     parent_vars2[i][j] = pl.LpVariable("parent_"+str(i)+"_"+str(j)+"_t2", 0, 1, pl.LpInteger)
        #     # all_path_tuples2.append([j,i])

        # print(all_path_tuples)
        # print("")
        
        for [n1,p1,n2,p2] in all_path_tuples:
            pm_vars [n1][p1][n2][p2] = pl.LpVariable("pm_"+str(n1)+"_"+str(p1)+"_"+str(n2)+"_"+str(p2), 0, 1, pl.LpInteger)
        prob = pl.LpProblem("treeEditDistance", pl.LpMinimize)
        
        constant1 = pl.LpVariable("one", 1, 1, pl.LpInteger)       
        insert_sum = 0
        ddescendants1 = [[] for x in range(tree1.number_of_nodes())] # descendants of vertex i in T1 including i
        ddescendants2 = [[] for x in range(tree2.number_of_nodes())] # descendants of vertex j in T2 including j
        for i in tree1.nodes:
            ddescendants1[i] = list(nx.descendants(oriented_tree1,i))+[i]
            for c in oriented_tree1.neighbors(i):
                insert_sum += d_insert(c,i)
            # insert_sum += d_insert(i)
        delete_sum = 0
        for i in tree2.nodes:
            ddescendants2[i] = list(nx.descendants(oriented_tree2,i))+[i]
            for c in oriented_tree2.neighbors(i):
                delete_sum += d_delete(c,i)
            # delete_sum += d_delete(i)
        # print(ddescendants1)
        # print(ddescendants2)
            
        # check that this is actually minimized
        # prob += pl.lpSum([(d_relabel(t[0], t[1]) - d_insert(t[0]) - d_delete(t[1])) * m_vars[t[0]][t[1]] for t in all_tuples] + [insert_sum * constant1, delete_sum * constant1])
        prob += pl.lpSum([(d_relabel(t[0],t[1],t[2],t[3]) - d_insert(t[0],t[1]) - d_delete(t[2],t[3])) * (pm_vars [t[0]][t[1]][t[2]][t[3]]) for t in all_path_tuples] + [insert_sum * constant1, delete_sum * constant1])
        
        if use_lower_bound:
            prob += pl.lpSum([(d_relabel(t[0],t[1],t[2],t[3]) - d_insert(t[0],t[1]) - d_delete(t[2],t[3])) * (pm_vars [t[0]][t[1]][t[2]][t[3]]) for t in all_path_tuples] + [insert_sum * constant1, delete_sum * constant1]) >= abs(insert_sum-delete_sum)


        # roots always need to be mapped
        if not isSubProblem:
            prob += m_vars[root1][root2] == 1
        
        
        if forced_pair != None:
            (p1, p2) = forced_pair
            prob += d_vars1[p1] + d_vars2[p2] + 2*m_vars[p1][p2] == 2

        for t in tree2.nodes:
            # every vertex may only be mapped once
            prob += pl.lpSum([m_vars[s][t] for s in tree1.nodes]) <= 1

            # deletion of node in T1 is one big NOR of all mappings to nodes in T2
            prob += pl.lpSum([m_vars[s][t] for s in tree1.nodes]) == 1 - d_vars2[t]

            # deletion of subtree is one big AND of deletions for all nodes in the subtree
            prob += 1 - dd_vars2[t] <= pl.lpSum([(1 - d_vars2[ct]) for ct in ddescendants2[t]])
            prob += pl.lpSum([(1 - d_vars2[ct]) for ct in ddescendants2[t]]) <= len(ddescendants2[t])*(1 - dd_vars2[t])

            if subtree_del_costs2[t] > upper_bound:
                prob += dd_vars2[t] == 0
            # if t!=root2 and d_delete(t,list(oriented_tree2.in_edges(t))[0][0]) > upper_bound:
            #         prob += d_vars2[t] == 0

            # at least two children without full subtree deletion
            prob += 2*two_children_not_dd_vars2[t] <= pl.lpSum([1-dd_vars2[ct] for ct in list(oriented_tree2.neighbors(t))])
            prob += pl.lpSum([1-dd_vars2[ct] for ct in list(oriented_tree2.neighbors(t))]+[-1]) <= len(list(oriented_tree2.neighbors(t)))*two_children_not_dd_vars2[t]
            
            # at least one child without full subtree deletion
            prob += one_child_not_dd_vars2[t] <= pl.lpSum([1-dd_vars2[ct] for ct in list(oriented_tree2.neighbors(t))])
            prob += pl.lpSum([1-dd_vars2[ct] for ct in list(oriented_tree2.neighbors(t))]) <= len(list(oriented_tree2.neighbors(t)))*one_child_not_dd_vars2[t]
            
            # pruned if exactly one remaining child and not deleted
            if len(list(oriented_tree2.neighbors(t))) > 1:
                prob += 1-pruned_vars2[t] <= pl.lpSum([1-one_child_not_dd_vars2[t],two_children_not_dd_vars2[t]])
                prob += pl.lpSum([1-one_child_not_dd_vars2[t],two_children_not_dd_vars2[t]]) <= 2*(1-pruned_vars2[t])
            else:
                prob += pruned_vars2[t] == 0
        
        for s in tree1.nodes:
            # every vertex may only be mapped once
            prob += pl.lpSum([m_vars[s][t] for t in tree2.nodes]) <= 1
            
            # deletion of node in T2 is one big NOR of all mappings to nodes in T1
            prob += pl.lpSum([m_vars[s][t] for t in tree2.nodes]) == 1 - d_vars1[s]
            
            # deletion of subtree is one big AND of deletions for all nodes in the subtree
            prob += 1 - dd_vars1[s] <= pl.lpSum([(1 - d_vars1[cs]) for cs in ddescendants1[s]])
            prob += pl.lpSum([(1 - d_vars1[cs]) for cs in ddescendants1[s]]) <= len(ddescendants1[s])*(1 - dd_vars1[s])

            if subtree_del_costs1[s] > upper_bound:
                prob += dd_vars1[s] == 0
            # if s!=root1 and d_insert(s,list(oriented_tree1.in_edges(s))[0][0]) > upper_bound:
            #         prob += d_vars1[s] == 0
            
            # at least two children without full subtree deletion
            prob += 2*two_children_not_dd_vars1[s] <= pl.lpSum([1-dd_vars1[cs] for cs in list(oriented_tree1.neighbors(s))])
            prob += pl.lpSum([1-dd_vars1[cs] for cs in list(oriented_tree1.neighbors(s))]+[-1]) <= len(list(oriented_tree1.neighbors(s)))*two_children_not_dd_vars1[s]
            
            # at least one child without full subtree deletion
            prob += one_child_not_dd_vars1[s] <= pl.lpSum([1-dd_vars1[cs] for cs in list(oriented_tree1.neighbors(s))])
            prob += pl.lpSum([1-dd_vars1[cs] for cs in list(oriented_tree1.neighbors(s))]) <= len(list(oriented_tree1.neighbors(s)))*one_child_not_dd_vars1[s]
            
            # pruned if exactly one remaining child and not deleted
            if len(list(oriented_tree1.neighbors(s))) > 1:
                prob += 1-pruned_vars1[s] <= pl.lpSum([1-one_child_not_dd_vars1[s],two_children_not_dd_vars1[s]])
                prob += pl.lpSum([1-one_child_not_dd_vars1[s],two_children_not_dd_vars1[s]]) <= 2*(1-pruned_vars1[s])
            else:
                prob += pruned_vars1[s] == 0
        
        # for [j,i] in all_path_tuples1_filtered:
        #     assert(is_ancestor1(i,j))
        #     prob += 1-parent_vars1[i][j] <= pl.lpSum([(1-pruned_vars1[s]) for s in predecessors1[i][j]]+[d_vars1[j],pruned_vars1[i],pruned_vars1[j]])
        #     prob += pl.lpSum([(1-pruned_vars1[s]) for s in predecessors1[i][j]]+[d_vars1[j],pruned_vars1[i],pruned_vars1[j]]) <= (len(predecessors1[i][j])+3)*(1-parent_vars1[i][j])

        # for [j,i] in all_path_tuples2_filtered:
        #     assert(is_ancestor2(i,j))
        #     prob += 1-parent_vars2[i][j] <= pl.lpSum([(1-pruned_vars2[s]) for s in predecessors2[i][j]]+[d_vars2[j],pruned_vars2[i],pruned_vars2[j]])
        #     prob += pl.lpSum([(1-pruned_vars2[s]) for s in predecessors2[i][j]]+[d_vars2[j],pruned_vars2[i],pruned_vars2[j]]) <= (len(predecessors2[i][j])+3)*(1-parent_vars2[i][j])
        
        for i in tree1.nodes:
            for j in tree1.nodes:
                
                # if i is not an ancestor of j or they are equal, i cannot be imaginary parent
                if (not is_ancestor1(i,j) or i==j): 
                    prob += parent_vars1[i][j] == 0
                
                # otherwise, i is imaginary parent of j if
                #   (a) all nodes in between (i.e. predecessors[i][j]) are pruned
                #    AND
                #   (b) j is neither pruned nor deleted and i is not pruned
                else: 
                    prob += 1-parent_vars1[i][j] <= pl.lpSum([(1-pruned_vars1[s]) for s in predecessors1[i][j]]+[d_vars1[j],pruned_vars1[i],pruned_vars1[j]])
                    prob += pl.lpSum([(1-pruned_vars1[s]) for s in predecessors1[i][j]]+[d_vars1[j],pruned_vars1[i],pruned_vars1[j]]) <= (len(predecessors1[i][j])+3)*(1-parent_vars1[i][j])

        for i in tree2.nodes:
            for j in tree2.nodes:
                
                # if i is not an ancestor of j or they are equal, i cannot be imaginary parent
                if (not is_ancestor2(i,j) or i==j): 
                    prob += parent_vars2[i][j] == 0
                
                # otherwise, i is imaginary parent of j if
                #   (a) all nodes in between (i.e. predecessors[i][j]) are pruned
                #    AND
                #   (b) j is neither pruned nor deleted and i is not pruned
                else: 
                    prob += 1-parent_vars2[i][j] <= pl.lpSum([(1-pruned_vars2[s]) for s in predecessors2[i][j]]+[d_vars2[j],pruned_vars2[i],pruned_vars2[j]])
                    prob += pl.lpSum([(1-pruned_vars2[s]) for s in predecessors2[i][j]]+[d_vars2[j],pruned_vars2[i],pruned_vars2[j]]) <= (len(predecessors2[i][j])+3)*(1-parent_vars2[i][j])
        
        # for i in tree1.nodes:
        #     # every vertex may only have one imaginary parent
        #     prob += pl.lpSum([parent_vars1[p][i] for p in predecessors1[i][root1]+[root1]]) <= 1
        # for j in tree2.nodes:
        #     # every vertex may only have one imaginary parent
        #     prob += pl.lpSum([parent_vars2[p][j] for p in predecessors2[j][root2]+[root2]]) <= 1
        
        for tuple1 in all_tuples:
            for tuple2 in all_tuples:
                # preserving ancestors
                if is_ancestor1(tuple1[0], tuple2[0]) != (is_ancestor2(tuple1[1], tuple2[1])):
                    prob += (m_vars[tuple1[0]][tuple1[1]] + m_vars[tuple2[0]][tuple2[1]]) <= 1
        
        # for [n1,p1] in all_path_tuples1:
        #     for [n2,p2] in all_path_tuples2:
        for [n1,p1,n2,p2] in all_path_tuples:
            # derive path mappings from node mappings and imaginary parents
            prob += 1-pm_vars [n1][p1][n2][p2] <= pl.lpSum([ 1-parent_vars1[p1][n1], 1-parent_vars2[p2][n2], 1-m_vars[n1][n2] ])
            prob += pl.lpSum([ 1-parent_vars1[p1][n1], 1-parent_vars2[p2][n2], 1-m_vars[n1][n2] ]) <= 3*(1-pm_vars [n1][p1][n2][p2])
        
        #if solution:
        #    for v in prob.variables():
        #        v.setInitialValue(solution[v.name])

        # prob.writeLP("test.lp")
        start = time.time()
        # solver = pl.getSolver('CPLEX_PY', timeLimit=time_limit, warmStart=True)
        if not isSubProblem:
            solver = pl.getSolver('GUROBI', timeLimit=time_limit, warmStart=False, msg=1, LogToConsole=1, OutputFlag=1, Threads=1)
        else:
            solver = pl.getSolver('GUROBI', timeLimit=time_limit, warmStart=False, msg=0, LogToConsole=0, OutputFlag=0, Threads=1)
        
        #solver.solver_parameter("Threads", 1) 
        # solver = pulp.CPLEX_PY()
        #solver.buildSolverModel(prob)
        # solver.solverModel.parameters.timelimit.set(time_limit)
        #solver.solverModel.parameters.mip.tolerances.absmipgap = 0.1
        #solver.solverModel.parameters.threads.set(1)
        #solver.callSolver(prob)
        #solver.SetSolverSpecificParametersAsString('Threads 1')
        #status = solver.findSolutionValues(prob)
        status = prob.solve(solver)
        #status = prob.solve()
        end = time.time()

        # root1_path_found = 0
        # root2_path_found = 0
        # for [n1,p1,n2,p2] in all_path_tuples:
        #     print(pm_vars [n1][p1][n2][p2].value())
        #     if pm_vars [n1][p1][n2][p2].value() == 1 and p1==root1:
        #         root1_path_found += 1
        #     if pm_vars [n1][p1][n2][p2].value() == 1 and p2==root2:
        #         root2_path_found += 1
        # assert(root1_path_found==1)
        # assert(root2_path_found==1)


        # print(prob.variables())
        
        # for i in tree1.nodes:
        #     for j in tree2.nodes:
        #         print(m_vars[i][j],m_vars[i][j].value())
        # print("")
        # for i in tree1.nodes:
        #     print(d_vars1[i],d_vars1[i].value())
        #     print(dd_vars1[i],dd_vars1[i].value())
        #     print(one_child_not_dd_vars1[i],one_child_not_dd_vars1[i].value())
        #     print(two_children_not_dd_vars1[i],two_children_not_dd_vars1[i].value())
        #     print(pruned_vars1[i],pruned_vars1[i].value())
        # print("")
        # for j in tree2.nodes:
        #     print(d_vars2[j],d_vars2[j].value())
        #     print(dd_vars2[j],dd_vars2[j].value())
        #     print(one_child_not_dd_vars2[j],one_child_not_dd_vars2[j].value())
        #     print(two_children_not_dd_vars2[j],two_children_not_dd_vars2[j].value())
        #     print(pruned_vars2[j],pruned_vars2[j].value())
        # print("")
        # print("")
        # # for i in tree1.nodes:
        # #     for j in tree1.nodes:
        # for [j,i] in all_path_tuples1:
        #     print(parent_vars1[i][j],parent_vars1[i][j].value())
        # print("")
        # # for i in tree2.nodes:
        # #     for j in tree2.nodes:
        # for [j,i] in all_path_tuples2:
        #     print(parent_vars2[i][j],parent_vars2[i][j].value())
        # print("")
        # # for [n1,p1] in all_path_tuples1:
        # #     for [n2,p2] in all_path_tuples2:
        # for [n1,p1,n2,p2] in all_path_tuples:
        #         if pm_vars [n1][p1][n2][p2].value() == 1.0:
        #             print("("+str(n1)+","+str(p1)+") - ("+str(n2)+","+str(p2)+") : ", d_relabel(n1,p1,n2,p2))
        # print("")
        # # for t in all_tuples:
        # #     if m_vars [t[0]][t[1]].value() == 1.0:
        # #         print(str(t[0])+"-"+str(t[1])+": ", d_relabel(t[0],t[1]))
        # print("")
        # print(end-start)

        # print(pl.LpStatus[status])
        #print(prob.objective)
        return (pl.LpStatus[status],prob.solverModel.ObjVal,prob.variables())

    labels1    = [0 for x in range(tree1.number_of_nodes())]
    tree1_type = [0 for x in range(tree1.number_of_nodes())]
    for i in tree1.nodes:
        labels1[i] = tree1.nodes(data=True)[i]["scalar"]
    for i in tree1.nodes:
        if tree1.degree[i] == 1 and (labels1[i] <= labels1[list(tree1.neighbors(i))[0]]):
            tree1_type[i] = 2
        if tree1.degree[i] == 1 and (labels1[i] >= labels1[list(tree1.neighbors(i))[0]]):
            tree1_type[i] = 1

    labels2 = [0 for x in range(tree2.number_of_nodes())]
    tree2_type = [0 for x in range(tree2.number_of_nodes())]
    for i in tree2.nodes:
        labels2[i] = tree2.nodes(data=True)[i]["scalar"]
    for i in tree2.nodes:
        if tree2.degree[i] == 1 and (labels2[i] <= labels2[list(tree2.neighbors(i))[0]]):
            tree2_type[i] = 2
        if tree2.degree[i] == 1 and (labels2[i] >= labels2[list(tree2.neighbors(i))[0]]):
            tree2_type[i] = 1
        
    oriented_tree1 = nx.dfs_tree(tree1, root1)
    oriented_tree2 = nx.dfs_tree(tree2, root2)
    
    def d_delete(n2,p2):
        return abs(labels2[n2]-labels2[p2])
        
    def d_insert(n1,p1):
        return abs(labels1[n1]-labels1[p1])
    
    insert_sum = 0
    for i in tree1.nodes:
        for c in oriented_tree1.neighbors(i):
            insert_sum += d_insert(c,i)
        # insert_sum += d_insert(i)
    delete_sum = 0
    for i in tree2.nodes:
        for c in oriented_tree2.neighbors(i):
            delete_sum += d_delete(c,i)

    time_limit = 10.0
    if initial_time_limit != None:
        time_limit = initial_time_limit
    upper_bound = delete_sum+insert_sum
    if known_upper_bound != None:
        upper_bound = min(upper_bound, known_upper_bound)
    best_solution = None

    while True:
        if not isSubProblem:
            print("")
            print("Computing IP for",tree1.number_of_nodes(),tree2.number_of_nodes(),"nodes")
            print("Time limit:",time_limit)
            print("Upper bound:",upper_bound)
        (status,res,sol) = rooted_unordered_deformation_distance_with_limits(tree1,root1,tree2,root2,upper_bound+0.01,time_limit,None,isSubProblem=isSubProblem,forced_pair=forced_pair, hard_time_limit=hard_time_limit)
        if not isSubProblem:
            print("Problem status:",status)
            print("Found solution value:",res)
        if status == "Optimal":
            return res,sol
        elif status == 103:
            return known_upper_bound, None
        elif status == "Not Solved":
            iteration += 1
            time_limit *= 2
            if hard_time_limit != None:
                if time_limit > hard_time_limit:
                    return None, None
            if res<upper_bound:
                upper_bound = res
                #best_solution = dict()
                #for v in sol:
                #    best_solution[v.name] = v.value()
        else:
            iteration += 1
            time_limit *= 2
            if hard_time_limit != None:
                if time_limit > hard_time_limit:
                    return None, None

def rooted_unordered_deformation_distance(tree1, root1, tree2, root2, isSubProblem=False,forced_pair=None,hard_time_limit=None,initial_time_limit=None,known_upper_bound=None):
    res,sol = rooted_unordered_deformation_distance_with_solution(tree1, root1, tree2, root2, isSubProblem=isSubProblem,forced_pair=forced_pair, hard_time_limit=hard_time_limit, initial_time_limit=initial_time_limit,known_upper_bound=known_upper_bound)
    return res
