# Reciprocity in Directed Hypergraphs: Measures, Findings, and Generators

We provide 
- **online appendix file**
- **source code**
for ***Reciprocity in Directed Hypergraphs: Measures, Findings, and Generators, which has been submitted to ICDM 22.***

In this work, we
- Propose a novel hypergraph reciprocity measure, ***HyperRec***.
- Propose an efficient searching algorithm for ***HyperRec***, a ***Ferret***.
- Discover real-world hypergraph's reciprocal pattern on 11 datasets.
- Propose an intuitive reciprocity-preserving hypergraph generator, ***ReDi***, which can captures the reciprocal patterns of real world hypergraphs.

***For online appendix, refer to the pdf file named [online_appendix_ICDM22_reciprocity_in_directed_hypergraphs.pdf]***

***For source codr, refer to the below descriptions.***

We provide source codes for following process. 

(1) Calculating proposed reciprocity ***HyperRec*** of the real-world hypergraphs.  
(2) Generating synthetic hypergraphs by using proposed hypergraph generator ***ReDi***, or null hypergraphs & baseline hypergraphs (see main paper for detail).  
(3) Reproducing all figures in the paper.    

## Requirements

To get requirements ready, run the following command on your terminal:
```
pip install -r requirements.txt
```

## Datasets

We provide 11 real-world hypergraph datasets
- metabolic_iaf1260b
- metabolic_ijo1366
- email_enron
- email_eu
- citation_dm
- citation_software
- qna_math
- qna_server
- bitcoin_2014
- bitcoin_2015
- bitcoin_2016

Please download dataset and put it inside the directory of src folder. For example
```
src
  |_data 
      |_iaf1260b
      |_ijo1366
      |_ ...
   main_figure.py
   ...
```

## Reproducing Reciprocity and Generation

To obtain results of reciprocity computation and generation, you may use main_reciprocity.py.  
Results will be saved in [/result] directory.

```
main_reciprocity.py [-d NAME_OF_DATASET] [-t TASK] [-a ALPHA_FOR_RECIPROCITY] [-b1 BETA1] [-b2 BETA2] [-in_n NUMBER_OF_INIT_NODES] [-gen_recip DIRECTLY_COMPUTING_RECIPROCITY_AFTER_GENERATION]

optional arguments:
-d NAME_OF_DATASET, Select name of dataset. String type with the exact name of the dataset (See Datasets subsection name).
-t TASK, Select task you want to proceed. String type with one of ["reciprocity", "redi_generation", "null_generation", "base_generation"]
  = "reciprocity" : Computes reciprocity
  = "redi_generation" : Generate ReDi's hypergraph (see main paper for ReDi)
  = "null_generation" : Generate Null hypergraph (see main paper for null hypergraph)
  = "base_generation" : Generate Baseline hypergraph (see main paper for baseline hypergraph)
-a ALPHA, Select alpha parameter. Float type with any alpha value greater than 0. (see main paper for alpha description)
-b1 BETA1, Select beta 1 parameter for ReDi. Float type between [0, 1]. (see main paper for beta1 description)
-b2 BETA2, Select beta 2 parameter for ReDi. Float type between [0, 1]. (see main paper for beta2 description)
-in_n NUMBER_OF_INIT_NODES, Select number of initial nodes parameter x, which implies (initial nodes = x times 2*max(H_{i}, T_{i}))
-gen_recip DIRECTLY_COMPUTING_RECIPROCITY_AFTER_GENERATION, Measuring generated hypergraphs right after the generation. Boolean type 
```

## Reproducing Figures

To obtain figures illustrated in the main paper and online appendix, you may use main_figure.py

As several dataset takes time in generation, we provide result of each generator in [/generated] directory.
Please download generated files and put them inside the directory of src folder. For example
```
src
  |_generated
      |_iaf1260b
      |_ijo1366
      |_ ...
   main_figure.py
   ...
```

Description regarding running is as follows.
Results sill be saved in [/figure] directory.

```
main_figure.py [-d NAME_OF_DATASET] [-t TASK] 
optional arguments:
-d NAME_OF_DATASET, Select name of dataset. String type with the exact name of the dataset (See Datasets subsection name).
-t TASK, Select task you want to proceed. String type with one of ["obs2", "obs3", "obs2_gen"]
  = "obs2" : Draw observation 2 figure of selected dataset (head-set out degree and tail-set in degree of zero & non-zero reciprocity hyperarcs)
  = "obs3" : Draw observation 3 figure of selected dataset (head-set out degree and tail-set in degree of zero & non-zero reciprocity hyperarcs)
  = "obs2_gen" : Draw observation 3 figure of selected dataset's generated hypergraph (ReDi).
```
