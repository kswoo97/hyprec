
import numpy as np
import matplotlib.pyplot as plt
from HyperRec import *
import os
import pandas as pd
from scipy.signal import savgol_filter


def compute_nodewise_reciprocity(head_index, tail_index, reciprocity):
    n_node = 0

    for e1, e2 in zip(head_index, tail_index) :

        tmp_max = max(max(e1), max(e2))

        if tmp_max > n_node :
            n_node = tmp_max

    n_node += 1

    my_edge_neighbors = {i: [] for i in range(n_node)}
    out_degree_dist = np.zeros(n_node)
    in_degree_dist = np.zeros(n_node)
    avg_recip = np.zeros(n_node)

    for s in range(len(head_index)):

        tmp_head = list(head_index[s])
        tmp_tail = list(tail_index[s])

        in_degree_dist[tmp_head] += 1
        out_degree_dist[tmp_tail] += 1

        for e1 in tmp_head:
            my_edge_neighbors[e1].extend([s])
        for e2 in tmp_tail:
            my_edge_neighbors[e2].extend([s])

    for v in range(len(my_edge_neighbors)):

        my_edges = my_edge_neighbors[v]
        if len(my_edges) > 0:
            avg_ = np.mean(reciprocity[my_edges])
            avg_recip[v] = avg_

    total_recip = in_degree_dist + out_degree_dist
    real_nodes = np.where(total_recip != 0)[0]
    return in_degree_dist[real_nodes], out_degree_dist[real_nodes], avg_recip[real_nodes]


def nodewise_reciprocity(in_d, out_d, n_r, param):
    diff = np.log(in_d + 1) - np.log(out_d + 1)
    N1 = pd.DataFrame({'diff': diff,
                       'recip': n_r})
    N2 = N1.groupby(by=['diff'], as_index=False).mean()
    X = N2['diff'].values
    Y = savgol_filter(N2.recip.values, param, 3)

    return X, Y


def return_XY(head_index, tail_index, reciprocity, param) :
    inD, outD, recip = compute_nodewise_reciprocity(head_index, tail_index, reciprocity)
    X, Y = nodewise_reciprocity(inD, outD, recip, param)

    return X, Y


def mean_gap(real_x, real_y, part_x, part_y) :
    overlapped_x = list(set(real_x).intersection(part_x))
    true_x_ind = list(np.where(np.isin(real_x, overlapped_x))[0])
    part_x_ind = list(np.where(np.isin(part_x, overlapped_x))[0])
    diffs = np.mean(np.abs(real_y[true_x_ind] - part_y[part_x_ind]))

    return diffs

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

def observation2(true_head, true_tail, true_recip, null_head, null_tail, null_recip,
                 title, saving_dir) :

    ticks = ["Real", "Null"]

    r_not = np.where(true_recip == 0)[0]
    r_is = np.where(true_recip != 0)[0]

    r_not2 = np.where(null_recip == 0)[0]
    r_is2 = np.where(null_recip != 0)[0]

    v_max = 0
    v_max2 = 0

    for i, j in zip(true_head, true_tail):

        tmp_max = max(max(i), max(j))

        if tmp_max > v_max:
            v_max = tmp_max

    for i2, j2 in zip(null_head, null_tail):

        tmp_max = max(max(i2), max(j2))

        if tmp_max > v_max2:
            v_max2 = tmp_max

    in_degrees = np.zeros(int(v_max + 1))
    out_degrees = np.zeros(int(v_max + 1))

    null_in_degrees = np.zeros(int(v_max2 + 1))
    null_out_degrees = np.zeros(int(v_max2 + 1))

    for i, j in zip(true_head, true_tail):
        in_degrees[list(i)] += 1
        out_degrees[list(j)] += 1

    for i2, j2 in zip(null_head, null_tail):
        null_in_degrees[list(i2)] += 1
        null_out_degrees[list(j2)] += 1

    edge_degree_vector_head = np.zeros(true_head.shape[0])
    null_edge_degree_vector_head = np.zeros(null_head.shape[0])

    edge_degree_vector_tail = np.zeros(true_head.shape[0])
    null_edge_degree_vector_tail = np.zeros(null_head.shape[0])

    idx = 0
    idx2 = 0

    for head_, tail_ in zip(true_head, true_tail):

        this_headD = np.mean(out_degrees[list(head_)])
        this_tailD = np.mean(in_degrees[list(tail_)])

        edge_degree_vector_head[idx] = this_headD
        edge_degree_vector_tail[idx] = this_tailD

        idx += 1

    for head_, tail_ in zip(null_head, null_tail) :

        this_headD = np.mean(null_out_degrees[list(head_)])
        this_tailD = np.mean(null_in_degrees[list(tail_)])

        null_edge_degree_vector_head[idx2] = this_headD
        null_edge_degree_vector_tail[idx2] = this_tailD

        idx2 += 1

    # We compute degree values of edge_degree_vector

    true_NNZ_head_degree = (edge_degree_vector_head[r_is])
    true_Z_head_degree = (edge_degree_vector_head[r_not])

    true_NNZ_tail_degree = (edge_degree_vector_tail[r_is])
    true_Z_tail_degree = (edge_degree_vector_tail[r_not])

    null_NNZ_head_degree = (null_edge_degree_vector_head[r_is2])
    null_Z_head_degree = (null_edge_degree_vector_head[r_not2])

    null_NNZ_tail_degree = (null_edge_degree_vector_tail[r_is2])
    null_Z_tail_degree = (null_edge_degree_vector_tail[r_not2])


    if title in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"] :
        zero_recips_head = [list(true_Z_head_degree + 1), list(null_Z_head_degree + 1)]
        non_zero_recips_head = [list(true_NNZ_head_degree + 1), list(null_NNZ_head_degree + 1)]

        zero_recips_tail = [list(true_Z_tail_degree + 1), list(null_Z_tail_degree + 1)]
        non_zero_recips_tail = [list(true_NNZ_tail_degree + 1), list(null_NNZ_tail_degree + 1)]

    else :
        zero_recips_head = [list(true_Z_head_degree), list(null_Z_head_degree)]
        non_zero_recips_head = [list(true_NNZ_head_degree), list(null_NNZ_head_degree)]

        zero_recips_tail = [list(true_Z_tail_degree), list(null_Z_tail_degree)]
        non_zero_recips_tail = [list(true_NNZ_tail_degree), list(null_NNZ_tail_degree)]

        plt.figure(figsize = (9, 3))

        plt.subplot(121)
        bpl = plt.boxplot(non_zero_recips_head , positions=[-1, -0.1], sym='', widths=0.3, )
        bpr = plt.boxplot(zero_recips_head , positions=[-0.6, 0.3], sym='', widths=0.3,  )
        set_box_color(bpl, '#D7191C')
        set_box_color(bpr, '#2C7BB6')

        plt.plot([], c='#D7191C', label='Nonzero')
        plt.plot([], c='#2C7BB6', label='Zero')

        plt.xticks([-0.8, 0.1], ticks, size = 25)
        if title in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"] :
            plt.yscale('log')
        plt.xlim(-1.2, 0.5)
        plt.tight_layout()
        plt.tick_params(axis='x', direction='in', size = 2)
        plt.tick_params(axis='y', direction='in', size = 2)
        plt.title(title + " Head Set Outdegree")
        plt.legend()

        plt.subplot(122)
        bpl = plt.boxplot(non_zero_recips_tail, positions=[-1, -0.1], sym='', widths=0.3, )
        bpr = plt.boxplot(zero_recips_tail, positions=[-0.6, 0.3], sym='', widths=0.3, )
        set_box_color(bpl, '#D7191C')
        set_box_color(bpr, '#2C7BB6')

        plt.plot([], c='#D7191C', label='Nonzero')
        plt.plot([], c='#2C7BB6', label='Zero')

        plt.xticks([-0.8, 0.1], ticks, size=25)
        if title in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]:
            plt.yscale('log')
        plt.xlim(-1.2, 0.5)
        plt.tight_layout()
        plt.tick_params(axis='x', direction='in', size=2)
        plt.tick_params(axis='y', direction='in', size=2)
        plt.title(title + " Tail Set Indegree")
        plt.legend()

        saving_path = os.path.join(saving_dir, title + "_Obs2.svg")
        plt.savefig(saving_path, format='svg', dpi=1200)
        plt.show()

def observation2_gen(gen_head, gen_tail, gen_recip, saving_dir, title) :

    r_not = np.where(gen_recip == 0)[0]
    r_is = np.where(gen_recip != 0)[0]

    v_max = 0

    ticks = ["Head Out", "Tail In"]

    for i, j in zip(gen_head, gen_tail):

        tmp_max = max(max(i), max(j))

        if tmp_max > v_max:
            v_max = tmp_max

    in_degrees = np.zeros(int(v_max + 1))
    out_degrees = np.zeros(int(v_max + 1))

    for i2, j2 in zip(gen_head, gen_tail):
        in_degrees[list(i2)] += 1
        out_degrees[list(j2)] += 1

    edge_head_vec = np.zeros(gen_head.shape[0])
    edge_tail_vec = np.zeros(gen_tail.shape[0])

    idx = 0
    idx2 = 0

    for head_, tail_ in zip(gen_head, gen_tail):

        this_inD = np.mean(out_degrees[list(head_)])
        this_outD = np.mean(in_degrees[list(tail_)])

        edge_head_vec[idx] = this_inD
        edge_tail_vec[idx] = this_outD

        idx += 1

    head_NNZ = edge_head_vec[r_is]
    head_Z = edge_head_vec[r_not]

    tail_NNZ = edge_tail_vec[r_is]
    tail_Z = edge_tail_vec[r_not]

    zero_recips = [list(head_Z), list(tail_Z)]
    non_zero_recips = [list(head_NNZ), list(tail_NNZ)]

    plt.figure(figsize = (5, 3))

    bpl = plt.boxplot(non_zero_recips , positions=[-1, -0.1], sym='', widths=0.3, )
    bpr = plt.boxplot(zero_recips , positions=[-0.6, 0.3], sym='', widths=0.3,  )

    set_box_color(bpl, '#D7191C')
    set_box_color(bpr, '#2C7BB6')

    plt.plot([], c='#D7191C', label='Nonzero')
    plt.plot([], c='#2C7BB6', label='Zero')

    plt.xticks([-0.8, 0.1], ticks, size = 15)
    plt.xlim(-1.2, 0.5)
    plt.tight_layout()
    plt.legend()

    plt.tick_params(axis='x', direction='in', size = 2)
    plt.tick_params(axis='y', direction='in', size = 2)
    saving_path = os.path.join(saving_dir, title + "_Obs2.svg")
    plt.savefig(saving_path, format='svg', dpi=1200)
    plt.show()


def observation3(true_head, true_tail, true_recip, gen_head, gen_tail, gen_recip,
                 null_head, null_tail, null_recip, base_head, base_tail, base_recip, saving_dir,
                 title) :

    param = 101 # Cab be changed if required

    real_X, real_Y = return_XY(true_head, true_tail, true_recip, param)
    gen_X, gen_Y = return_XY(gen_head, gen_tail, gen_recip, param)
    null_X, null_Y = return_XY(null_head, null_tail, null_recip, 5)
    base_X, base_Y = return_XY(base_head, base_tail, base_recip, 5)

    plt.figure(figsize=(5, 5))
    ## Seed 20
    plt.plot(real_X, real_Y, color="dodgerblue", linewidth=3, label = "Real World")
    plt.plot(gen_X, gen_Y, color="forestgreen", linewidth=3, label = "ReDi")
    plt.plot(base_X, base_Y, color="orange", linewidth=3, label = "Baseline")
    plt.plot(null_X, null_Y, color="red", linewidth=3, label = "Null")
    plt.legend()
    plt.xlim(-5, 5)
    plt.xticks([-3, 0, 3], size=15)
    plt.tight_layout()
    plt.tick_params(axis='x', direction='in', size=4, width=2)
    saving_directory = os.path.join(saving_dir, title + "_Obs3.svg")
    plt.savefig(saving_directory, format='svg', dpi=1200)
    plt.show()

    print("Mean Gap between Real World vs ReDi: ", mean_gap(real_X, real_Y, gen_X, gen_Y))
    print("Mean Gap between Real World vs Null: ", mean_gap(real_X, real_Y, null_X, null_Y))
    print("Mean Gap between Real World vs Baseline: ", mean_gap(real_X, real_Y, base_X, base_Y))