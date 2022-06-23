## This file is to measure reciprocity for a given dataset

import sys
import argparse
from figure import *
from ReDi import * # Hypergraph Generator
import pickle

if __name__ == "__main__" :

    '''
    plot type supports three ; 
    1. type = 'obs2' ; this gives box-plot
    2. type = 'obs3' ; this gives line plot of three different cases
    3. type = 'obs2_gen' ; this gives box-plot of generated graphs
    '''

    parser = argparse.ArgumentParser("Official Impelmentation for Directed Hypergraph Reciprocity Figure")
    parser.add_argument("-d", "--dataset", type = str, help = "Select dataset for observing result")
    parser.add_argument("-t", "--fig_type", type = str, help = "Select which type of plot you wish to create")

    args = parser.parse_args()
    data_lists = ["metabolic_iaf1260b", "metabolic_ijo1366", "email_enron", "email_eu", "citation_dm", "citation_software",
                  "qna_math", "qna_server", "bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]

    if not os.path.exists(os.path.join(os.getcwd(), "figures")) :
        os.mkdir(os.path.join(os.getcwd(), "figures"))

    figure_path = os.path.join(os.getcwd(), "figures")

    if args.dataset in data_lists :
        data_path = "data/{0}/{0}".format(args.dataset)
        recip_path = "reciprocity/{0}/{0}_reciprocity_00.npy".format(args.dataset)
        gen_path = "generated/{0}/{0}".format(args.dataset)

    else :
        print("Unable to interpret dataset type. Please check dataset name once again")
        sys.exit(0)

######## Load four types of reciproicity

    with open(data_path + "_head.pickle", 'rb') as handle:
        true_head = np.array(pickle.load(handle))
    with open(data_path + "_tail.pickle", 'rb') as handle:
        true_tail = np.array(pickle.load(handle))
    true_recip = np.load(recip_path, allow_pickle=True)

    with open(gen_path + "_redi_head.pickle", 'rb') as handle:
        redi_head = np.array(pickle.load(handle))
    with open(gen_path + "_redi_tail.pickle", 'rb') as handle:
        redi_tail = np.array(pickle.load(handle))
    redi_recip = np.load(gen_path + "_redi_recip.npy", allow_pickle=True)

    with open(gen_path + "_base_head.pickle", 'rb') as handle:
        base_head = np.array(pickle.load(handle))
    with open(gen_path + "_base_tail.pickle", 'rb') as handle:
        base_tail = np.array(pickle.load(handle))
    base_recip = np.load(gen_path + "_base_recip.npy", allow_pickle=True)

    with open(gen_path + "_null_head.pickle", 'rb') as handle:
        null_head = np.array(pickle.load(handle))
    with open(gen_path + "_null_tail.pickle", 'rb') as handle:
        null_tail = np.array(pickle.load(handle))
    null_recip = np.load(gen_path + "_null_recip.npy", allow_pickle=True)

    if args.fig_type == "obs2" :
        observation2(true_head, true_tail, true_recip, null_head, null_tail, null_recip,
                     args.dataset, figure_path)
        sys.exit(0)
    elif args.fig_type == "obs3" :
        observation3(true_head, true_tail, true_recip, redi_head, redi_tail, redi_recip,
                     null_head, null_tail, null_recip, base_head, base_tail, base_recip,
                     figure_path, args.dataset)
        sys.exit(0)
    elif args.fig_type == "obs2_gen":
        observation2_gen(redi_head, redi_tail, redi_recip, figure_path, args.dataset)
        sys.exit(0)

    else :
        print("Task should be given one of obs2, obs3, obs2_gen")
        sys.exit(0)