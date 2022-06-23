# Reciprocity in Directed Hypergraphs: Measures, Findings, and Generators

We provide source code for Reciprocity in Directed Hypergraphs: Measures, Findings, and Generators, which has been submitted to ICDM 22.  

In this work, we
- Propose a novel hypergraph reciprocity measure, ***HyperRec***.
- Propose an efficient searching algorithm for ***HyperRec***, a ***Ferret***.
- Discover real-world hypergraph's reciprocal pattern on 11 datasets.
- Propose an intuitive reciprocity-preserving hypergraph generator, ***ReDi***, which can captures the reciprocal patterns of real world hypergraphs.

Here, we provide source codes for following process. 

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

## Reproduction of Reciprocity and Generation

To obtain the results of reciprocity computation and generation, you may use main_reciprocity.py

```
main_reciprocity.py [-d NAME_OF_DATASET] [-t TASK] [-a ALPHA_FOR_RECIPROCITY] [-b1 BETA1] [-b2 BETA2] [-in_n NUMBER_OF_INIT_NODES] [-gen_recip DIRECTLY_COMPUTING_RECIPROCITY_AFTER_GENERATION]
```
