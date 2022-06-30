## Tools for HyperRc Measurement

import numpy as np

import itertools
from tqdm import tqdm
from scipy.spatial import distance
from scipy import sparse


def returning_head_tail_matrix(head_index, tail_index, as_sparse = True) :

    if as_sparse :
        head_rows, head_cols = [], []
        tail_rows, tail_cols = [], []

        for i, s1 in enumerate(head_index) :

            for s2 in s1 :

                head_rows.append(s2)
                head_cols.append(i)

        for j, t1 in enumerate(tail_index) :

            for t2 in t1 :

                tail_rows.append(t2)
                tail_cols.append(j)

        s_values = np.repeat(1, len(head_rows))
        t_values = np.repeat(1, len(tail_rows))

        n_dim = max(head_rows + tail_rows) + 1

        shead = sparse.csr_matrix((s_values, (head_rows, head_cols)),
                                  shape = (n_dim, len(head_index)))
        stail = sparse.csr_matrix((t_values, (tail_rows, tail_cols)),
                                  shape = (n_dim, len(tail_index)))

    else :

        head_rows, head_cols = [], []
        tail_rows, tail_cols = [], []

        for i, s1 in enumerate(head_index):

            for s2 in s1:
                head_rows.append(s2)
                head_cols.append(i)

        for j, t1 in enumerate(tail_index):

            for t2 in t1:
                tail_rows.append(t2)
                tail_cols.append(j)

        n_dim = max(head_rows + tail_rows) + 1

        shead = np.zeros((n_dim, head_index.shape[0]))
        stail = np.zeros((n_dim, tail_index.shape[0]))

        shead[(head_rows, head_cols)] = 1
        stail[(tail_rows, tail_cols)] = 1

    return shead, stail

def efficient_distance_matching(true_S, true_T, e, candidates, head_index, tail_index,
                                is_sparse=True, alpha=1e-6,
                               special_case = False) :

    MAX = np.log(2)
    max_recip = 0

    if len(candidates) == 1 : # 1 vs 1 simple matching

        if is_sparse : # Sparse matrix
            true_probs = np.array(true_T[:, e].todense()).reshape(1, -1)[0] / np.sum(true_T[:, e])
            candid = candidates[0]
            our_dist = np.array(true_S[:, candid].todense()).reshape(1, -1)[0] / np.sum(true_S[:, candid])

        else :  # Dense matrix
            true_probs = np.array(true_T[:, e]) / np.sum(true_T[:, e])
            candid = candidates[0]
            our_dist = np.array(true_S[:, candid]) / np.sum(true_S[:, candid])

        penalty = (distance.jensenshannon(true_probs, our_dist)) ** 2
        S1 = len(head_index[e].intersection(tail_index[candid]))
        S2 = len(head_index[e]) - S1

        max_recip = 1 - (S1 * penalty + S2 * MAX) / (len(head_index[e]) * MAX)

    else : # Sparse Matrix
        # Step 1. Reduce V -> V ; only nodes who exist in arcs
        if is_sparse :
            true_probs = np.array(true_T[:, e].todense()).reshape(1, -1)[0] / np.sum(true_T[:, e])
            partialS = set(np.where(true_S[:, candidates].todense() > 0)[0]) # Target Arc Nodes are inside here
            partialT = set(np.where(true_probs > 0)[0]) # Candidate Nodes are inside here
            idx = list(partialS.union(partialT))
            S = np.array(true_S[idx, :].todense())  # V' x |E|

        else :
            true_probs = np.array(true_T[:, e]) / np.sum(true_T[:, e])
            partialS = set(np.where(true_S[:, candidates] > 0)[0])
            partialT = set(np.where(true_probs > 0)[0])
            idx = list(partialS.union(partialT))
            S = true_S[idx, :]

        if S.ndim == 1 :
            S = S.reshape(-1, 1)
        true_probs = true_probs[idx]

        head_tail_intersect = np.zeros((len(head_index[e]),
                                        len(candidates)))  ## |S| * |C| Matrix

        for i, head_ind in enumerate(head_index[e]) :
            for j, tail_ind in enumerate(candidates) :
                if head_ind in tail_index[tail_ind] :
                    head_tail_intersect[i, j] = 1 # Intersection among head & tail set nodes


        full_probs = np.array(S[:, candidates])
        full_probs = full_probs / np.sum(full_probs, axis=0)
        full_probs = full_probs.transpose()  # C x V'

        entire_candids = np.arange(len(candidates))
        
        if special_case == "FDHG" : # Iteratively add one by one

            each_head_size = [len(head_index[x]) for x in candidates]
            candid_e = []

            for cardin_ in range(1, (len(candidates) + 1)):

                this_time_add = np.argmin(each_head_size)
                candid_e.append(this_time_add)
                each_head_size[this_time_add] = 1000000

                if cardin_ == 1 :

                    new_S = head_tail_intersect[:, candid_e].reshape(-1)  # |S| * |C'|
                    new_P = full_probs[candid_e, :]  # C' x V'

                    new_P = new_P.reshape(-1)
                    try:
                        penalty = (distance.jensenshannon(true_probs, new_P)) ** 2
                    except:
                        penalty = MAX

                    cur_recip = (1 - (penalty * np.sum(new_S) + MAX * (new_S.shape[0] - 1)) / (  # This has shape of 1
                                MAX * new_S.shape[0]))*((1/cardin_)**alpha)

                    if cur_recip > max_recip :
                        max_recip = cur_recip
                    else :
                        break

                else:  # more than 1

                    new_S = head_tail_intersect[:, candid_e]  # |S| * |C'|
                    new_P = full_probs[candid_e, :]  # C' x V'

                    new_S = new_S.transpose()
                    new_S /= np.sum(new_S, axis=0)
                    new_S = new_S.transpose()  # |S| * |C'|
                    if np.sum(np.isnan(new_S)) > 0:
                        new_S[np.isnan(new_S)] = 0

                    new_P = np.matmul(new_S, new_P)  # |S| * |V'|
                    penalty = 0

                    for tmp_e in range(new_P.shape[0]):

                        if np.sum(new_P[tmp_e, :]) == 0:
                            penalty += MAX

                        else:
                            try:
                                penalty += (distance.jensenshannon(true_probs, new_P[tmp_e, :])) ** 2
                            except:
                                penalty += MAX

                    cur_recip = (1 - (penalty) / (new_P.shape[0] * MAX))*(((1 / len(candid_e)) ** alpha))

                    if cur_recip > max_recip :
                        max_recip = cur_recip
                    else :
                        break
            
        else : # General Case
            
            for cardin_ in range(1, (len(candidates) + 1)):

                cur_candid = list(itertools.combinations(entire_candids, cardin_))

                if cardin_ == 1 :
                    for combin in cur_candid:
                        candid_e = list(combin)

                        new_S = head_tail_intersect[:, candid_e].reshape(-1)  # |S| * |C'|
                        new_P = full_probs[candid_e, :]  # C' x V'

                        new_P = new_P.reshape(-1)
                        try:
                            penalty = (distance.jensenshannon(true_probs, new_P)) ** 2
                        except:
                            penalty = MAX

                        cur_recip = 1 - (penalty * np.sum(new_S) + MAX * (new_S.shape[0] - np.sum(new_S))) / (MAX * new_S.shape[0])

                        if cur_recip > max_recip:
                            max_recip = cur_recip

                else:  # more than 1
                    
                    for combin in cur_candid:
                        candid_e = list(combin)
                        new_S = head_tail_intersect[:, candid_e]  # |S| * |C'|
                        new_P = full_probs[candid_e, :]  # C' x V'

                        new_S = new_S.transpose()
                        new_S /= np.sum(new_S, axis=0)
                        new_S = new_S.transpose()  # |S| * |C'|
                        if np.sum(np.isnan(new_S)) > 0 : 
                            new_S[np.isnan(new_S)] = 0

                        new_P = np.matmul(new_S, new_P)  # |S| * |V'|
                        penalty = 0

                        for tmp_e in range(new_P.shape[0]):

                            if np.sum(new_P[tmp_e, :]) == 0:
                                penalty += MAX

                            else:
                                try:
                                    penalty += (distance.jensenshannon(true_probs, new_P[tmp_e, :])) ** 2
                                except:
                                    penalty += MAX
                                    

                        cur_recip = 1 - (penalty) / (new_P.shape[0] * MAX)
                        cur_recip *= (1/len(candid_e))**alpha

                        if cur_recip > max_recip:
                            max_recip = cur_recip

    return max_recip

def measure_reciprocity(head_index, tail_index, data_type = 'normal_sparse',
                       omega_dist = {}, given_interval = np.arange(10), alpha = 1e-6, special_case = "normal") :

    new_omega = np.zeros(head_index.shape[0])
    reciprocity = np.zeros(len(head_index))
    S, T = returning_head_tail_matrix(head_index, tail_index, as_sparse=True)

    if data_type == 'normal_sparse' :

        source = S.transpose().dot(T)  # E x E of Checking start
        entire = source.multiply(T.transpose().dot(S))

        omega_dist = np.array(np.sum(entire, axis=1)).reshape(1, -1)[0]

        iterator_ = tqdm(range(omega_dist.shape[0]))
            
        for e in iterator_ :
            
            is_pass = False

            if omega_dist[e] == 0:  # No Overlapping
                continue
            
            else:
                true_headset = head_index[e]  # Headset |Hi|
                true_tailset = tail_index[e]  # Tailset |Ti|

                head_idx = {j: i for i, j in enumerate(true_headset)}
                tail_idx = {j: i for i, j in enumerate(true_tailset)}
                
                candids = np.where(np.array(entire[e, :].todense())[0] > 0)[0]
                H_ = np.zeros((candids.shape[0], (len(true_headset)) + len(true_tailset)))
                size_order = np.zeros((candids.shape[0], 2))  # Head and tail size

                for cur_e, c_ in enumerate(candids):  # Iteratively moving candid

                    head_over = true_headset.intersection(tail_index[c_])
                    tail_over = true_tailset.intersection(head_index[c_])

                    for i_ in head_over:
                        H_[cur_e, head_idx[i_]] = 1

                    for j_ in tail_over:
                        H_[cur_e, len(head_idx) + tail_idx[j_]] = 1

                    size_order[cur_e, 0] = len(head_index[c_])
                    size_order[cur_e, 1] = len(tail_index[c_])

                tmp_entire = np.unique(H_, axis=0, return_inverse=True)[1]
                candids_list = []

                for cur_edge in np.unique(tmp_entire):
                    inds = np.where(tmp_entire == cur_edge)[0]  # They are group of edges which overlap on the same place
                    head_min = np.argmin(size_order[inds, 0])
                    candids_list.append(candids[inds[head_min]])

                candids_list = list(set(candids_list))

                for c_e in candids_list :

                    tmp_head = head_index[c_e]
                    tmp_tail = tail_index[c_e]

                    if (tmp_head == true_tailset) & (tmp_tail == true_headset) :

                        reciprocity[e] = 1
                        is_pass = True

                if is_pass :
                    continue

                max_recip = efficient_distance_matching(e = e,
                                                candidates = candids_list,
                                                head_index = head_index,
                                                tail_index = tail_index,
                                                true_S = S,
                                                true_T = T,
                                                is_sparse=True,
                                                alpha = alpha,
                                                special_case = special_case)

                reciprocity[e] = max_recip
                
    elif data_type == "bitcoin" : 
        
        reciprocity = np.zeros(given_interval.shape[0])
        new_omega = np.zeros(given_interval.shape[0])
        starting_point = min(given_interval)
        
        for e in tqdm(given_interval) :

            is_pass = False

            if omega_dist[e] == [] :  # No Overlapping
                continue
            
            else:
                
                true_headset = head_index[e]  # Headset |Hi|
                true_tailset = tail_index[e]  # Tailset |Ti|

                head_idx = {j: i for i, j in enumerate(true_headset)}
                tail_idx = {j: i for i, j in enumerate(true_tailset)}
                
                candids = omega_dist[e]
                H_ = np.zeros((len(candids), (len(true_headset)) + len(true_tailset)))
                size_order = np.zeros((len(candids), 2))  # Head and tail size

                for cur_e, c_ in enumerate(candids):  # Iteratively moving candid

                    head_over = true_headset.intersection(tail_index[c_])
                    tail_over = true_tailset.intersection(head_index[c_])

                    for i_ in head_over:
                        H_[cur_e, head_idx[i_]] = 1

                    for j_ in tail_over:
                        H_[cur_e, len(head_idx) + tail_idx[j_]] = 1

                    size_order[cur_e, 0] = len(head_index[c_])
                    size_order[cur_e, 1] = len(tail_index[c_])

                tmp_entire = np.unique(H_, axis=0, return_inverse=True)[1]
                candids_list = []

                for cur_edge in np.unique(tmp_entire):
                    inds = np.where(tmp_entire == cur_edge)[0]  # They are group of edges which overlap on the same place
                    head_min = np.argmin(size_order[inds, 0])
                    candids_list.append(candids[inds[head_min]])

                candids_list = list(set(candids_list))

                for c_e in candids_list :

                    tmp_head = head_index[c_e]
                    tmp_tail = tail_index[c_e]

                    if (tmp_head == true_tailset) & (tmp_tail == true_headset) :

                        reciprocity[e] = 1
                        is_pass = True

                if is_pass :
                    continue

                max_recip = efficient_distance_matching(e = e,
                                            candidates = candids_list,
                                            head_index = head_index,
                                            tail_index = tail_index,
                                            true_S = S,
                                            true_T = T,
                                            is_sparse=True,
                                            alpha = alpha,
                                            special_case = special_case)
                idx = int(e - starting_point)
                reciprocity[idx] = max_recip
    
    else : # Normal dense matrix
        
        source = np.matmul(S.transpose(), T)  # E x E of Checking start
        entire = source * np.matmul(T.transpose(), S)
        omega_dist = np.sum(entire, axis=1)
        

        iterator_ = tqdm(range(omega_dist.shape[0]))

        for e in iterator_ :

            is_pass = False
            
            max_recip = 0

            if omega_dist[e] == 0:  # No Overlapping
                continue
                
            else:
                true_headset = head_index[e]  # Headset |Hi|
                true_tailset = tail_index[e]  # Tailset |Ti|

                head_idx = {j: i for i, j in enumerate(true_headset)}
                tail_idx = {j: i for i, j in enumerate(true_tailset)}
                candids = np.where(entire[e, :] > 0)[0]  # Overlapped Edges
                H_ = np.zeros((candids.shape[0], (len(true_headset)) + len(true_tailset)))
                size_order = np.zeros((candids.shape[0], 2))  # Head and tail size

                for cur_e, c_ in enumerate(candids):  # Iteratively moving candid

                    head_over = true_headset.intersection(tail_index[c_])
                    tail_over = true_tailset.intersection(head_index[c_])

                    for i_ in head_over:
                        H_[cur_e, head_idx[i_]] = 1

                    for j_ in tail_over:
                        H_[cur_e, len(head_idx) + tail_idx[j_]] = 1

                    size_order[cur_e, 0] = len(head_index[c_])
                    size_order[cur_e, 1] = len(tail_index[c_])

                tmp_entire = np.unique(H_, axis=0, return_inverse=True)[1]
                candids_list = []

                for cur_edge in np.unique(tmp_entire):
                    inds = np.where(tmp_entire == cur_edge)[0]  # They are group of edges which overlap on the same place
                    head_min = np.argmin(size_order[inds, 0])
                    candids_list.append(candids[inds[head_min]])

                candids_list = list(set(candids_list))

                for c_e in candids_list :

                    tmp_head = head_index[c_e]
                    tmp_tail = tail_index[c_e]

                    if (tmp_head == true_tailset) & (tmp_tail == true_headset) :

                        reciprocity[e] = 1
                        is_pass = True

                if is_pass :
                    continue
                
                max_recip = efficient_distance_matching(e = e,
                                            candidates = candids_list, 
                                            head_index = head_index, 
                                            tail_index = tail_index,
                                            true_S = S, 
                                            true_T = T, 
                                            is_sparse=False, 
                                            alpha = alpha,
                                            special_case = special_case)
                

                reciprocity[e] = max_recip

    return reciprocity


def return_search_space(S, T, head_index, tail_index, data_type='normal_sparse',
                        omega_dist={}, given_interval=np.arange(10)):

    original_omega = np.zeros(head_index.shape[0])
    new_omega = np.zeros(head_index.shape[0])
    return_candids = []

    if data_type == 'normal_sparse':

        source = S.transpose().dot(T)  # E x E of Checking start
        entire = source.multiply(T.transpose().dot(S))

        omega_dist = np.array(np.sum(entire, axis=1)).reshape(1, -1)[0]

        for e in tqdm(range(omega_dist.shape[0])):

            max_recip = 0
            if omega_dist[e] == 0:  # No Overlapping
                original_omega[e] = 0
                new_omega[e] = 0
                return_candids.extend([0])

            else:

                true_headset = head_index[e]  # Headset |Hi|
                true_tailset = tail_index[e]  # Tailset |Ti|

                head_idx = {j: i for i, j in enumerate(true_headset)}
                tail_idx = {j: i for i, j in enumerate(true_tailset)}

                candids = np.where(np.array(entire[e, :].todense())[0] > 0)[0]
                H_ = np.zeros((candids.shape[0], (len(true_headset)) + len(true_tailset)))
                size_order = np.zeros((candids.shape[0], 2))  # Head and tail size

                for cur_e, c_ in enumerate(candids):  # Iteratively moving candid

                    head_over = true_headset.intersection(tail_index[c_])
                    tail_over = true_tailset.intersection(head_index[c_])

                    for i_ in head_over:
                        H_[cur_e, head_idx[i_]] = 1

                    for j_ in tail_over:
                        H_[cur_e, len(head_idx) + tail_idx[j_]] = 1

                    size_order[cur_e, 0] = len(head_index[c_])
                    size_order[cur_e, 1] = len(tail_index[c_])

                tmp_entire = np.unique(H_, axis=0, return_inverse=True)[1]
                candids_list = []

                for cur_edge in np.unique(tmp_entire):
                    inds = np.where(tmp_entire == cur_edge)[
                        0]  # They are group of edges which overlap on the same place
                    head_min = np.argmin(size_order[inds, 0])
                    candids_list.append(candids[inds[head_min]])

                candids_list = list(set(candids_list))

                original_omega[e] = candids.shape[0]
                new_omega[e] = len(candids_list)
                return_candids.append(candids_list)

    elif data_type == "bitcoin":

        reciprocity = np.zeros(given_interval.shape[0])
        new_omega = np.zeros(given_interval.shape[0])
        starting_point = min(given_interval)

        for e in tqdm(given_interval):

            max_recip = 0
            if omega_dist[e] == []:  # No Overlapping
                new_omega[e] = 0
                original_omega[e] = 0
                return_candids.append([0])
                continue

            else:

                true_headset = head_index[e]  # Headset |Hi|
                true_tailset = tail_index[e]  # Tailset |Ti|

                head_idx = {j: i for i, j in enumerate(true_headset)}
                tail_idx = {j: i for i, j in enumerate(true_tailset)}

                candids = omega_dist[e]
                H_ = np.zeros((len(candids), (len(true_headset)) + len(true_tailset)))
                size_order = np.zeros((len(candids), 2))  # Head and tail size

                for cur_e, c_ in enumerate(candids):  # Iteratively moving candid

                    head_over = true_headset.intersection(tail_index[c_])
                    tail_over = true_tailset.intersection(head_index[c_])

                    for i_ in head_over:
                        H_[cur_e, head_idx[i_]] = 1

                    for j_ in tail_over:
                        H_[cur_e, len(head_idx) + tail_idx[j_]] = 1

                    size_order[cur_e, 0] = len(head_index[c_])
                    size_order[cur_e, 1] = len(tail_index[c_])

                tmp_entire = np.unique(H_, axis=0, return_inverse=True)[1]
                candids_list = []

                for cur_edge in np.unique(tmp_entire):
                    inds = np.where(tmp_entire == cur_edge)[
                        0]  # They are group of edges which overlap on the same place
                    head_min = np.argmin(size_order[inds, 0])
                    candids_list.append(candids[inds[head_min]])

                candids_list = list(set(candids_list))

                original_omega[e] = len(omega_dist)
                new_omega[e] = len(candids_list)
                return_candids.append(candids_list)

    else:  # Normal dense matrix

        source = np.matmul(S.transpose(), T)  # E x E of Checking start
        entire = source * np.matmul(T.transpose(), S)
        omega_dist = np.sum(entire, axis=1)

        for e in tqdm(range(omega_dist.shape[0])):

            max_recip = 0

            if omega_dist[e] == 0:  # No Overlapping
                new_omega[e] = 0
                original_omega[e] = 0
                return_candids.append([0])

            else:
                true_headset = head_index[e]  # Headset |Hi|
                true_tailset = tail_index[e]  # Tailset |Ti|

                head_idx = {j: i for i, j in enumerate(true_headset)}
                tail_idx = {j: i for i, j in enumerate(true_tailset)}
                candids = np.where(entire[e, :] > 0)[0]  # Overlapped Edges
                H_ = np.zeros((candids.shape[0], (len(true_headset)) + len(true_tailset)))
                size_order = np.zeros((candids.shape[0], 2))  # Head and tail size

                for cur_e, c_ in enumerate(candids):  # Iteratively moving candid

                    head_over = true_headset.intersection(tail_index[c_])
                    tail_over = true_tailset.intersection(head_index[c_])

                    for i_ in head_over:
                        H_[cur_e, head_idx[i_]] = 1

                    for j_ in tail_over:
                        H_[cur_e, len(head_idx) + tail_idx[j_]] = 1

                    size_order[cur_e, 0] = len(head_index[c_])
                    size_order[cur_e, 1] = len(tail_index[c_])

                tmp_entire = np.unique(H_, axis=0, return_inverse=True)[1]
                candids_list = []

                for cur_edge in np.unique(tmp_entire):
                    inds = np.where(tmp_entire == cur_edge)[
                        0]  # They are group of edges which overlap on the same place
                    head_min = np.argmin(size_order[inds, 0])
                    candids_list.append(candids[inds[head_min]])

                candids_list = list(set(candids_list))

                new_omega[e] = len(candids_list)
                original_omega[e] = candids.shape[0]
                return_candids.append(candids_list)

    return original_omega, new_omega, return_candids

def bitcoin_dataset_reciprocity_computation(head_index, tail_index, alpha) : #

    ## As bitcoin dataset is extremely big, do efficient computation
    true_head_idx = head_index
    true_tail_idx = tail_index

    head_rows = []
    head_cols = []

    tail_rows = []
    tail_cols = []

    for i in (range((head_index.shape[0]))) :

        for head_ in head_index[i]:
            head_rows.append(head_)
            head_cols.append(i)

        for tail_ in tail_index[i]:
            tail_rows.append(tail_)
            tail_cols.append(i)

    row_vals = np.ones(len(head_rows)).astype(int)
    col_vals = np.ones(len(tail_rows)).astype(int)

    E = head_index.shape[0]
    V = int(max(max(head_rows), max(tail_rows))) + 1

    S = sparse.csr_matrix((row_vals,
                    (head_rows, head_cols)), shape=(V, E))
    T = sparse.csr_matrix((col_vals,
                    (tail_rows, tail_cols)), shape=(V, E))

    interval = np.arange(0, E + 200000, 200000)
    omega = {i: [] for i in range(E)}

    for i in (range(interval.shape[0] - 1)):

        print("Preprocessing {0}/{1} is on process".format(i + 1, interval.shape[0] - 1))
        E1 = S.T[int(interval[i]): int(interval[i + 1])].dot(T)
        E2 = T.T[int(interval[i]): int(interval[i + 1])].dot(S)
        E1, E2 = E1.multiply(E2).nonzero()

        for e1, e2 in (zip(E1, E2)):
            if e1 != e2:
                omega[int(e1 + i * 200000)].extend([e2])

    prev_sum = 0
    interval_ranges = {}
    entire_reciprocity = np.zeros(head_index.shape[0])

    for k in range(40):
        if prev_sum + 50000 > head_index.shape[0]:
            end_point = head_index.shape[0]
            cur_intervals = np.arange(prev_sum, end_point)
            interval_ranges[k] = cur_intervals
            prev_sum += 50000
            break
        else:
            end_point = prev_sum + 50000
            cur_intervals = np.arange(prev_sum, end_point)
            interval_ranges[k] = cur_intervals
            prev_sum += 50000

    print("Preprocessing Done. Reciprocity Computation...")

    for i in range(len(interval_ranges)):
        try:
            cur_interval = interval_ranges[i]
            print("{0}/{1} process is on".format(i, len(interval_ranges)))
            true_omega1, true_omega2, true_recip, candid_ = measure_reciprocity(true_head_idx, true_tail_idx,
                                                                                data_type='bitcoin',
                                                                                omega_dist=omega,
                                                                                given_interval=cur_interval,
                                                                                alpha=alpha)
            entire_reciprocity[cur_interval] = true_recip
        except:
            continue

    return entire_reciprocity
