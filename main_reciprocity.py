## This file is to measure reciprocity for a given dataset

import os
import sys
import argparse

from HyperRec import * # Reciprocity Measuring
from ReDi import * # Hypergraph Generator

if __name__ == "__main__" :

    parser = argparse.ArgumentParser("Official Impelmentation for Directed Hypergraph Reciprocity")
    parser.add_argument("-d", "--dataset", type = str, help = "Select dataset for reciprocity measuring or generation")
    parser.add_argument("-t", "--task", type = str, help = "Select task you want to process")
    parser.add_argument("-a",
                        "--alpha",
                        type = float,
                        default=1e-6,
                        help = "Choose alpha parameter for reciprocity calculation of HyperRec")

    parser.add_argument("-b1",
                        "--beta1",
                        type = float,
                        default = 0.1,
                        help = "Choose beta 1 parameter for generation of ReDi")

    parser.add_argument("-b2",
                        "--beta2",
                        type=float,
                        default=0.1,
                        help="Choose beta 2 parameter for generation of ReDi")

    parser.add_argument("-in_n",
                        "--inN",
                        type=float,
                        default=1,
                        help="Choose number of initial nodes parameter x, which is x times 2*max(|H_{i}|, |T_{i}|)")

    parser.add_argument("-gen_recip",
                        "--gen_recip",
                        type=bool,
                        default=False,
                        help="Measuring generated hypergraphs reciprocity right after the generation")

    args = parser.parse_args()
    data_lists = ["metabolic_iaf1260b", "metabolic_ijo1366", "email_enron", "email_eu", "citation_dm", "citation_software",
                  "qna_math", "qna_server", "bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]

    if not os.path.exists(os.path.join(os.getcwd(), "results")) :
        os.mkdir(os.path.join(os.getcwd(), "results"))

    if args.dataset in data_lists :
        data_path = "data/{0}/{0}".format(args.dataset)
    else :
        print("Unable to interpret dataset type. Please check dataset name once again")
        sys.exit(0)

######## Real-World Hypergraph Reciprocity Computation

    with open(data_path + "_head.pickle", 'rb') as handle:
        head_index = np.array(pickle.load(handle))
    with open(data_path + "_tail.pickle", 'rb') as handle:
        tail_index = np.array(pickle.load(handle))

    if args.task == "reciprocity" :

        if args.dataset in ["email_enron", "email_eu"] :
            result = measure_reciprocity(head_index = head_index,
                                         tail_index = tail_index,
                                         alpha=args.alpha,
                                         special_case="FDHG")

        elif args.dataset in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"] : # Requires scalable addition steps
            result = bitcoin_dataset_reciprocity_computation(head_index = head_index,
                                                             tail_index = tail_index,
                                                             alpha = args.alpha)
        else :
            result = measure_reciprocity(head_index=head_index,
                                         tail_index=tail_index,
                                         alpha=args.alpha,
                                         special_case="normal")

        print("Dataset : {0}".format(args.dataset)) ; print()
        print("Computed Reciprocity : {0}".format(np.mean(result)))
        saving_reciprocity = os.path.join(os.getcwd(), "results" , args.dataset + "_reciprocity_{0}.npy".format(args.alpha))
        np.save(saving_reciprocity, result)

    # Generating ReDi Hypergraph

    elif args.task == "redi_generation" :

        if args.dataset in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"] : # Degreewise generation of graph
            gen_head, gen_tail = ReDi_degreewise(head_index=head_index,
                                                 tail_index=tail_index,
                                                 beta1 = args.beta1,
                                                 beta2 = args.beta2,
                                                 mode = args.inN)

        else :
            gen_head, gen_tail = ReDi(head_index=head_index,
                                                 tail_index=tail_index,
                                                 beta1=args.beta1,
                                                 beta2=args.beta2,
                                                mode=args.inN)

        head_path = os.path.join(os.getcwd(), "results" ,
                                 args.dataset + "_redi_head.pickle")
        tail_path = os.path.join(os.getcwd(), "results" ,
                                 args.dataset + "_redi_tail.pickle")

        with open(head_path, 'wb') as handle:
            pickle.dump(list(gen_head), handle)

        with open(tail_path, 'wb') as handle:
            pickle.dump(list(gen_tail), handle)

        if args.gen_recip :

            if args.dataset in ["email_enron", "email_eu"]:
                result = measure_reciprocity(head_index=gen_head,
                                             tail_index=gen_tail,
                                             alpha=args.alpha,
                                             special_case="FDHG")

            elif args.dataset in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]:  # Requires scalable addition steps
                result = bitcoin_dataset_reciprocity_computation(head_index=gen_head,
                                                                tail_index=gen_tail,
                                                                 alpha=args.alpha)
            else:
                result = measure_reciprocity(head_index=gen_head,
                                             tail_index=gen_tail,
                                             alpha=args.alpha,
                                             special_case="normal")

            print("ReDi Generated Dataset : {0}".format(args.dataset));
            print()
            print("Computed Reciprocity : {0}".format(np.mean(result)))
            saving_reciprocity = os.path.join(os.getcwd(), "results",
                                              args.dataset + "_redi_reciprocity_{0}.npy".format(args.alpha))
            np.save(saving_reciprocity, result)

    # Generating Null Hypergraph

    elif args.task == "null_generation" :

        gen_head, gen_tail = naive_sampling(head_index = head_index,
                       tail_index = tail_index,
                       seed = 0)

        head_path = os.path.join(os.getcwd(), "results" ,
                                 args.dataset + "_null_head.pickle")
        tail_path = os.path.join(os.getcwd(), "results" ,
                                 args.dataset + "_null_tail.pickle")

        with open(head_path, 'wb') as handle:
            pickle.dump(list(gen_head), handle)

        with open(tail_path, 'wb') as handle:
            pickle.dump(list(gen_tail), handle)

        if args.gen_recip :

            if args.dataset in ["email_enron", "email_eu"]:
                result = measure_reciprocity(head_index=gen_head,
                                             tail_index=gen_tail,
                                             alpha=args.alpha,
                                             special_case="FDHG")

            elif args.dataset in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]:  # Requires scalable addition steps
                result = bitcoin_dataset_reciprocity_computation(head_index=gen_head,
                                                                tail_index=gen_tail,
                                                                 alpha=args.alpha)
            else:
                result = measure_reciprocity(head_index=gen_head,
                                             tail_index=gen_tail,
                                             alpha=args.alpha,
                                             special_case="normal")

            print("Null Generated Dataset : {0}".format(args.dataset));
            print()
            print("Computed Reciprocity : {0}".format(np.mean(result)))
            saving_reciprocity = os.path.join(os.getcwd(), "results",
                                              args.dataset + "_null_reciprocity_{0}.npy".format(args.alpha))
            np.save(saving_reciprocity, result)

    # Generating Baseline Hypergraph

    elif args.task == "baseline_generation" :

        if args.dataset in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]:  # Degreewise generation of graph
            gen_head, gen_tail = ReDi_degreewise(head_index=head_index,
                                                 tail_index=tail_index,
                                                 beta1=0.0,
                                                 beta2=0.0,
                                                mode = args.inN)

        else:
            gen_head, gen_tail = ReDi(head_index=head_index,
                                      tail_index=tail_index,
                                      beta1=0.0,
                                      beta2=0.0,
                                      mode = args.inN)

        head_path = os.path.join(os.getcwd(), "results" ,
                                 args.dataset + "_base_head.pickle")
        tail_path = os.path.join(os.getcwd(), "results" ,
                                 args.dataset + "_base_tail.pickle")

        with open(head_path, 'wb') as handle:
            pickle.dump(list(gen_head), handle)

        with open(tail_path, 'wb') as handle:
            pickle.dump(list(gen_tail), handle)

        if args.gen_recip :

            if args.dataset in ["email_enron", "email_eu"]:
                result = measure_reciprocity(head_index=gen_head,
                                             tail_index=gen_tail,
                                             alpha=args.alpha,
                                             special_case="FDHG")

            elif args.dataset in ["bitcoin_2014", "bitcoin_2015", "bitcoin_2016"]:  # Requires scalable addition steps
                result = bitcoin_dataset_reciprocity_computation(head_index=gen_head,
                                                                tail_index=gen_tail,
                                                                 alpha=args.alpha)
            else:
                result = measure_reciprocity(head_index=gen_head,
                                             tail_index=gen_tail,
                                             alpha=args.alpha,
                                             special_case="normal")

            print("Baseline Generated Dataset : {0}".format(args.dataset));
            print()
            print("Computed Reciprocity : {0}".format(np.mean(result)))
            saving_reciprocity = os.path.join(os.getcwd(), "results",
                                              args.dataset + "_base_reciprocity_{0}.npy".format(args.alpha))
            np.save(saving_reciprocity, result)

    else :
        print("Task should be given either reciprocity or generation")
        sys.exit(0)