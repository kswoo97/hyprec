import numpy as np
from tqdm import tqdm
from itertools import combinations

'''
How total_tuple looks like?

{1 : {(0) : 1, # size : Node : degree
      (1) : 3, },

 2 : {(0, 1) : 3, # size : Node : degree
    ...}
 }
'''


def naive_sampling(head_index, tail_index, seed):  # Naive Baseline Method!
    np.random.seed(seed)
    new_head = []
    new_tail = []

    head_size, tail_size = [], []

    v_max = 0

    total_head_nodes = 0
    total_tail_nodes = 0

    for i in range(head_index.shape[0]):

        cur_head_size = len(head_index[i])
        cur_tail_size = len(tail_index[i])

        total_head_nodes += cur_head_size
        total_tail_nodes += cur_tail_size

        head_size.append(cur_head_size)
        tail_size.append(cur_tail_size)
        union_max = max(head_index[i].union(tail_index[i]))
        if v_max < union_max:
            v_max = union_max

    in_degree = np.zeros(int(v_max + 1))
    out_degree = np.zeros(int(v_max + 1))

    for i in range(head_index.shape[0]):
        cur_head = list(head_index[i])
        cur_tail = list(tail_index[i])

        in_degree[cur_head] += 1
        out_degree[cur_tail] += 1

    total_degree = in_degree + out_degree
    required_index = np.where(total_degree != 0)[0]

    entire_nodes = np.arange(required_index.shape[0])

    entire_head_random = np.random.choice(entire_nodes, total_head_nodes)
    entire_tail_random = np.random.choice(entire_nodes, total_tail_nodes)

    head_indpt = 0
    tail_indpt = 0

    for i in (range(head_index.shape[0])):

        cond = True

        h_size = head_size[i]
        t_size = tail_size[i]

        next_head_indpt = head_indpt + h_size
        next_tail_indpt = tail_indpt + t_size

        cur_H = set(entire_head_random[head_indpt:next_head_indpt])
        cur_T = set(entire_tail_random[tail_indpt:next_tail_indpt])

        cond1 = len(cur_H) == h_size
        cond2 = len(cur_T) == t_size
        cond3 = len(cur_H.intersection(cur_T)) == 0

        if cond1 * cond2 * cond3 == 1:

            None

        else:
            while cond:

                cur_H = set(np.random.choice(entire_nodes, h_size, replace=False))
                cur_T = set(np.random.choice(entire_nodes, t_size, replace=False))

                if len(cur_H.intersection(cur_T)) == 0:
                    cond = False

        new_head.append((cur_H))
        new_tail.append((cur_T))
        head_indpt = next_head_indpt
        tail_indpt = next_tail_indpt

    return np.array(new_head), np.array(new_tail)

def find_interaction(total_tuple, valid_nodes, required_size, degree_list, delta):
    try:
        deg_lists = np.array(list(total_tuple[required_size].values()))
        probs = (deg_lists) / np.sum(deg_lists)
        c_inter = int(np.random.choice(a=np.arange(probs.shape[0]),
                                       size=1,
                                       p=probs,
                                       replace=False))
        add_list = np.array(list(total_tuple[required_size].keys()))[c_inter]
    except:
        valid_nodes = list(set(valid_nodes))
        deg_list = degree_list[:len(valid_nodes)].copy() + delta
        deg_p = deg_list / np.sum(deg_list)
        add_list = np.random.choice(a=valid_nodes,
                                    p=deg_p,
                                    size=required_size, replace=False)

    return set(add_list)


def update_interaction(total_tuple, new_edge):
    new_E = list(new_edge)

    for i in range(1, int(len(new_E) + 1)):
        possible_comb = combinations(new_E, i)
        for t1 in possible_comb:
            try:
                total_tuple[i][t1] += 1
            except:
                total_tuple[i][t1] = 1

    return total_tuple


def get_recip_set(degree, candid_list, NPick, delta):
    deg_ = degree[list(candid_list)].copy()
    deg_ += delta
    prob = deg_ / np.sum(deg_)
    chosen_list = set(np.random.choice(a=list(candid_list),
                                       size=NPick,
                                       p=prob,
                                       replace=False))

    return chosen_list


def ReDi(head_index, tail_index, beta1, beta2, mode
         , seed = 0):
    v_max = 0
    head_size = []
    tail_size = []
    np.random.seed(seed)
    delta = 1

    for i in range(head_index.shape[0]):

        union_max = max(head_index[i].union(tail_index[i]))
        if v_max < union_max:
            v_max = union_max
        if ((len(head_index[i])) <= 10) & ((len(tail_index[i])) <= 10):
            head_size.extend([len(head_index[i])])
            tail_size.extend([len(tail_index[i])])

    node_filter1 = np.zeros(int(v_max + 1))
    node_filter2 = np.zeros(int(v_max + 1))

    for heads, tails in zip(head_index, tail_index):
        node_filter1[list(heads)] += 1
        node_filter2[list(tails)] += 1

    condition_indicator = node_filter1 + node_filter2
    total_N = np.where(condition_indicator != 0)[0].shape[0]  # Number of nodes

    add_N = max(int(len(head_index) - total_N), 0)
    NP = np.ones(total_N)

    tmp_NP = np.ones(v_max + 1)

    if mode == 0:
        add_pos = np.random.choice(a=np.arange(int(tmp_NP.shape[0])),
                                   size=add_N)

    else:
        for h_, t_ in zip(head_index, tail_index):
            this_node = max(max(h_), max(t_))
            tmp_NP[this_node] += 1

        tmpP = tmp_NP / np.sum(tmp_NP)
        np.random.shuffle(tmpP)

        add_pos = np.random.choice(a=np.arange(int(tmpP.shape[0])),
                                   p=tmpP,
                                   size=add_N)

    for i in add_pos:
        NP[i] += 1  # Add base degree

    each_add_N = NP

    total_heads = []
    total_tails = []

    ## Initialization
    max_S = max(np.max(head_size), np.max(tail_size))

    head_bucket = {i: {} for i in range(int(max_S + 1))}  # Where we save interactions
    tail_bucket = {i: {} for i in range(int(max_S + 1))}  # Where we save interactions

    current_nodes = []

    total_indegree = np.zeros(total_N)
    total_outdegree = np.zeros(total_N)

    is_recip_vector = np.random.choice(a=[True, False],
                                       size=int(np.sum(each_add_N)),
                                       p=[beta1, 1 - beta1])

    size_vectors = np.random.choice(np.arange(len(head_size)), size=int(np.sum(each_add_N) * 3))

    if mode > 1:
        init_n = (max_S * mode)

    else:
        init_n = int(max_S)

    init_nodes = np.arange(int(init_n * 2))

    for i in range(int(init_n)):
        tmp_head_list = {init_nodes[int(i)]}
        tmp_tail_list = {init_nodes[int(init_n + i)]}
        total_heads.append(tmp_head_list)
        total_tails.append(tmp_tail_list)
        head_bucket = update_interaction(total_tuple=head_bucket,
                                         new_edge=tmp_head_list)
        tail_bucket = update_interaction(total_tuple=tail_bucket,
                                         new_edge=tmp_tail_list)
        current_nodes.append(int(i))
        current_nodes.append(int(init_n + i))
        total_indegree[list(tmp_head_list)] += 1
        total_outdegree[list(tmp_tail_list)] += 1

    idx = -1
    for i in tqdm(range(total_N)):

        k = each_add_N[i]  # Decide number of arcs to be created
        indptr = len(total_heads)
        for j in range(int(k)):

            cond = True
            idx += 1
            idx2 = idx - 1
            recip = is_recip_vector[idx]

            while cond:
                idx2 += 1
                cur_idx = size_vectors[idx2]
                N_head = int(head_size[cur_idx])
                N_tail = int(tail_size[cur_idx])

                is_head = np.random.choice([True, False], 1)

                if recip:  # Let arc reciprocal for this time
                    # If first iteration ; Nothing but choice one edge and make as reciprocal

                    if j == 0:

                        to_be_reciprocal = np.arange(len(total_heads))
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    else:

                        to_be_reciprocal = np.arange(indptr, len(total_heads))  # Reciprocal candidates
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))  # Choose one of them
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    max_tail_recips = min(N_tail, len(reciprocal_head)) - 1
                    max_head_recips = min(N_head, len(reciprocal_tail)) - 1

                    if max_head_recips > 1:
                        Nhead_reciprocal = np.random.binomial(n=max_head_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    else:
                        Nhead_reciprocal = 1

                    if max_tail_recips > 1:
                        Ntail_reciprocal = np.random.binomial(n=max_tail_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    else:
                        Ntail_reciprocal = 1

                    head_idx = set(get_recip_set(total_indegree,
                                                 reciprocal_tail,
                                                 Nhead_reciprocal,
                                                 delta))

                    tail_idx = set(get_recip_set(total_outdegree,
                                                 reciprocal_head,
                                                 Ntail_reciprocal,
                                                 delta))

                    if (len(head_idx) == N_head) & (len(tail_idx) == N_tail):

                        None  # Nothing to do

                    elif len(head_idx) == N_head:
                        tail_idx = tail_idx.union({i})

                    else:
                        head_idx = head_idx.union({i})

                    N_add_head = N_head - len(head_idx)
                    N_add_tail = N_tail - len(tail_idx)

                    if N_add_head > 0:
                        part_head = find_interaction(head_bucket,
                                                     current_nodes, N_add_head,
                                                     total_indegree, delta)
                        head_idx = head_idx.union(part_head)

                    if N_add_tail > 0:
                        part_tail = find_interaction(tail_bucket,
                                                     current_nodes, N_add_tail,
                                                     total_outdegree, delta)
                        tail_idx = tail_idx.union(part_tail)

                else:  # Naively modify
                    if is_head:
                        part_head = find_interaction(head_bucket,
                                                     current_nodes, int(N_head - 1),
                                                     total_indegree, delta)
                        head_idx = {i}.union(part_head)
                        tail_idx = find_interaction(tail_bucket, current_nodes, int(N_tail),
                                                    total_outdegree, delta)
                    else:
                        head_idx = find_interaction(head_bucket, current_nodes, int(N_head),
                                                    total_indegree, delta)
                        part_tail = find_interaction(tail_bucket, current_nodes, int(N_tail - 1),
                                                     total_outdegree, delta)
                        tail_idx = {i}.union(part_tail)

                if len(head_idx.intersection(tail_idx)) < 1:
                    cond = False
                    head_bucket = update_interaction(head_bucket, head_idx)  # In-interaction
                    tail_bucket = update_interaction(tail_bucket, tail_idx)  # Out-Interaction
                    total_indegree[list(head_idx)] += 1  # In-degree
                    total_outdegree[list(tail_idx)] += 1  # Out-degree
                    total_heads.append(head_idx)
                    total_tails.append(tail_idx)

        current_nodes.append(i)

    return np.array(total_heads), np.array(total_tails)

def ReDi_degreewise(head_index, tail_index, beta1, beta2, mode,
                    seed = 0) :

    v_max = 0
    head_size = np.zeros(head_index.shape[0])
    tail_size = np.zeros(tail_index.shape[0])
    delta = 1

    np.random.seed(seed)

    for i in range(head_index.shape[0]):

        union_max = max(head_index[i].union(tail_index[i]))
        if v_max < union_max:
            v_max = union_max
        if len(head_index[i]) <= 10:
            head_size[i] = len(head_index[i])
        else :
            head_size[i] = 10
        if len(head_index[i]) <= 10:
            tail_size[i] = len(tail_index[i])
        else :
            tail_size[i] = 10

    node_filter1 = np.zeros(int(v_max + 1))
    node_filter2 = np.zeros(int(v_max + 1))

    for heads, tails in zip(head_index, tail_index):
        node_filter1[list(heads)] += 1
        node_filter2[list(tails)] += 1

    condition_indicator = node_filter1 + node_filter2
    total_N = np.where(condition_indicator != 0)[0].shape[0]  # Number of nodes

    add_N = max(int(len(head_index) - total_N), 0)
    NP = np.ones(int(v_max) + 1)

    if mode == 0:  # mode 0 indicates naive addition of edges
        add_pos = np.random.choice(a=np.arange(int(total_N - 1)),
                                   size=add_N,
                                   replace=True)

        for i in add_pos:
            NP[i] += 1  # Add degree one by one

        each_add_N = np.random.choice(NP, total_N, replace=True)

    else:
        tmp_NP = np.ones(NP.shape[0])
        for h_, t_ in zip(head_index, tail_index):
            this_node = max(max(h_), max(t_))
            tmp_NP[this_node] += 1

        tmpP = tmp_NP / np.sum(tmp_NP)
        np.random.shuffle(tmpP)

        add_pos = np.random.choice(a=np.arange(int(tmpP.shape[0])),
                                   p=tmpP,
                                   size=add_N)

        for i in add_pos:
            NP[i] += 1  # Add degree one by one

        each_add_N = np.random.choice(NP, total_N, replace=True)

    total_heads = []
    total_tails = []

    in_norm = 0
    out_norm = 0

    ## Initialization
    max_S = max(np.max(head_size), np.max(tail_size))
    if mode > 1:
        init_n = int(max_S) * mode

    else:
        init_n = int(max_S)
    init_nodes = np.arange(int(init_n * 2))

    current_nodes = []

    total_indegree = np.zeros(total_N)
    total_outdegree = np.zeros(total_N)
    IND = int(init_n * 2)

    for i in range(int(init_n)):
        tmp_head_list = {init_nodes[int(i)]}
        tmp_tail_list = {init_nodes[int(init_n + i)]}

        total_heads.append(tmp_head_list)
        total_tails.append(tmp_tail_list)

        in_norm += len(tmp_head_list)
        out_norm += len(tmp_head_list)

        total_indegree[list(tmp_head_list)] += 1
        total_outdegree[list(tmp_tail_list)] += 1

    is_recip_vector = np.random.choice(a=[True, False],
                                            size=int(np.sum(each_add_N)),
                                            p=[beta1, 1 - beta1])

    size_vectors = np.random.choice(np.arange(head_index.shape[0]),
                                         size=int(np.sum(each_add_N) * 3))
    idx = -1
    idx2 = -1

    for i in tqdm(each_add_N):

        if i > IND:
            IND = i

        k = each_add_N[i]  # Decide number of arcs to be created
        indptr = len(total_heads)
        for j in range(int(k)):
            idx += 1
            cond = True
            recip = is_recip_vector[idx] # Check whether an arc is reciprocal
            idx2 = idx - 1
            while cond:
                idx2 += 1
                cur_idx = size_vectors[idx2]
                N_head = int(head_size[cur_idx])
                N_tail = int(tail_size[cur_idx])

                is_head = np.random.choice([True, False], 1)

                if recip:  # Let this arc reciprocal
                    # If first iteration ; Nothing but choice one edge and make as reciprocal

                    if j == 0:
                        to_be_reciprocal = np.arange(len(total_heads))
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    else:
                        to_be_reciprocal = np.arange(indptr, len(total_heads))  # Reciprocal Candiate
                        recip_idx = int(np.random.choice(to_be_reciprocal, 1))  # Choose one among them
                        reciprocal_head = total_heads[recip_idx]
                        reciprocal_tail = total_tails[recip_idx]

                    max_tail_recips = min(N_tail, len(reciprocal_head)) - 1
                    max_head_recips = min(N_head, len(reciprocal_tail)) - 1

                    if max_head_recips > 1:
                        Nhead_reciprocal = np.random.binomial(n=max_head_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    else:
                        Nhead_reciprocal = 1

                    if max_tail_recips > 1:
                        Ntail_reciprocal = np.random.binomial(n=max_tail_recips,
                                                              p=beta2,
                                                              size=1) + 1
                    else:
                        Ntail_reciprocal = 1

                    head_idx = set(get_recip_set(total_indegree,
                                                 reciprocal_tail,
                                                 Nhead_reciprocal,
                                                 delta))

                    tail_idx = set(get_recip_set(total_outdegree,
                                                 reciprocal_head,
                                                 Ntail_reciprocal,
                                                 delta))

                    if (len(head_idx) == N_head) & (len(tail_idx) == N_tail):

                        None  # Nothing to do

                    elif len(head_idx) == N_head:
                        tail_idx = tail_idx.union({i})

                    else:
                        head_idx = head_idx.union({i})

                    N_add_head = N_head - len(head_idx)
                    N_add_tail = N_tail - len(tail_idx)

                    if N_add_head > 0:
                        deg = total_indegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_head = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=N_add_head,
                                                              replace=False)))
                        head_idx = head_idx.union(part_head)

                    if N_add_tail > 0:
                        deg = total_outdegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_tail = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=N_add_tail,
                                                              replace=False)))
                        tail_idx = tail_idx.union(part_tail)

                else:  # Let this arc be naively created
                    if is_head:
                        deg = total_indegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_head = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=int(N_head - 1),
                                                              replace=False)))
                        head_idx = {i}.union(part_head)

                        deg = total_outdegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        tail_idx = set(list(np.random.choice(a=np.arange(IND),
                                                             p=prop,
                                                             size=int(N_tail),
                                                             replace=False)))
                    else:
                        deg = total_indegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        head_idx = set(list(np.random.choice(a=np.arange(IND),
                                                             p=prop,
                                                             size=int(N_head),
                                                             replace=False)))

                        deg = total_outdegree[:IND] + delta
                        prop = deg / np.sum(deg)
                        part_tail = set(list(np.random.choice(a=np.arange(IND),
                                                              p=prop,
                                                              size=int(N_tail - 1),
                                                              replace=False)))
                        tail_idx = {i}.union(part_tail)

                if len(head_idx.intersection(tail_idx)) < 1:
                    cond = False

                    total_indegree[list(head_idx)] += 1  # In-degree
                    total_outdegree[list(tail_idx)] += 1  # Out-degree
                    total_heads.extend([head_idx])
                    total_tails.extend([tail_idx])

        current_nodes.append(i)

    return np.array(total_heads), np.array(total_tails)