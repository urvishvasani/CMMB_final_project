# Multiple Sequence Alignment using multi-objective optimization approach

unzip the "MSA_G3.zip" file. 

## Content of MSA_G3

1. MSA - Julia package for multiple sequence alignment
2. MSA.ipynb - Jupyter notebook for multiple sequence alignment
3. input1.txt - sample input file.


## 1. Steps to generate alignments from the julia package : MSA

To import the package and generate the alignments, execute the following commands from the julia console:

1. cd("MSA_G3")
2. cd("MSA")
3. press ] to switch to package mode.
4. activate .
5. press ctrl+c to switch back to julia console.
6. cd("..")
7. import MSA
8. msa.generate_alignments("input1.txt")

Here, FASTA formatted input is present inside "input1.txt". Output alignments are stored in the "output1.txt" file.


## 2. Steps to execute the Jupyter Notebook: MSA.ipynb

* To run the Jupyter Notebook, Open the notebook "MSA.ipynb" and hit "cells -> run all" if you are running from jupyter notebook. 
* Code reads the input sequences from the "input1.txt" and writes the generated alignments back to the "output1.txt".
* To provide the custom file for input, change the filename in the function call of "generate_alignments".
* This notebook also contains the code to generate plots which demonstrates the performance of the proposed approach for various aspects. This code is not present inside the julia package.


### Note: Detailed description of each functions are present inside jupyter markdown. 