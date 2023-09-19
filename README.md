# ClassGraph
ClassGraph2 is a tool that classifies metagenomic reads starting from the output of other pre-existing binning tools

Current technologies allow the sequencing of microbial communities directly from the environment without prior culturing. One of the major problems when analyzing a microbial sample is to taxonomically annotate its reads to identify the species it contains. Most of the methods currently available focus on the classification of reads using a set of reference genomes and their k-mers. While in terms of precision these methods have reached percentages of correctness close to perfection, in terms of recall (the actual number of classified reads) the performances fall at around 50%. One of the reasons is the fact that the sequences in a sample can be very different from the corresponding reference genome. To address this issue, we propose ClassGraph, a new binning  method  that  makes  use  of  the  reads  overlap  graph  and  applies  a label propagation algorithm to refine the binning result of existing tools.

<p align="center">
  <img src="Images/CompletePipeline.png" width="600" title="ClassGraph workflow" alt="ClassGraph workflow">
</p>

## ClassGraph2 Download
It's possible to download ClassGraph2 by cloning the repository in your machine.

```
git clone https://github.com/MargheritaCavattoni/ClassGraph2.git
```
## Installation
The istallation of ClassGraph2 requires python 3.7 or above. Besides the following dependency is needed:
* [python-igraph](https://igraph.org/python/)

### Setting the environment with Conda
It's possible to use [Conda](https://docs.conda.io/en/latest/) in order to set-up the proper environment of ClassGraph2. [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) contain Conda.

In order to create the Conda environment and activate it run the following commands:
```
conda env create -f environment.yml
conda activate ClassGraph2
```

## Preprocessing
ClassGraph2 requires two input ﬁles: one representing a graph of reads and the other containing the result of the classiﬁcation process. The labelles assigned to the reads by the pre-existing binning tool will be propagated over the graph to the still unclassified reads. ClassGraph2 is thought to be used with both paired-end and single-end reads. It can be used for short and for long reads.

### Overlap Graph
The graph must be composed as follows: the nodes represent the reads and the weights of the edges their overlaps. The file must be presented in asqg format.
One possible tool that satisfies these requirements is [**SGA**](https://github.com/jts/sga) (String Graph Assembler), an assembler based on the overlap-layout-consensus.
If you are using an assembler that does not use exact overlaps to build the graph, you can modify the asqg format by adding a final column containing the values of the weights you want to assign to the edges (in this case set the option --assembler 2).

### Binning Tools
ClassGraph2 requires as input the result of the classification process presented as follows:
```
NODE.1  L1
NODE.2  L1
NODE.3  L1
NODE.4  L2
NODE.5  L2
...
```
Where NODE_i is the node ID and Li represents its label. Li must be equal to 0 if the read wasn't classified.
Hypothetically any metagenomic binning tool could be used for the classification. 

## Usage Example
In order to run ClassGraph2 it's necessary to type on the terminal a command structured as follows:
```
python3 ClassGraph2.py --graph $Assembler_DIR/Graph.asqg --output $ClassGraph2_Output_DIR/ --binned $Classifier_DIR/BinnedReads.out --preﬁx example1 --max_iteration 20 --read_type 1 --assembler 1
```
Where:
```
--graph           path to the graph file (asqg)
--output          path to the folder where we want to store the output
--binned          path to the file that stores the result of the initial classification with a preexisting binning tool
--prefix          name of the output file
--max_iteration   maximum number of iterations in the label propagation algorithm (default: 20)
--read_type       1 for paired-end (default), 2 for single-end reads
--assembler       1 for asqg format (default), 2 in the case of modified asqg
```
## Getting help
If you encounter bugs or have further questions or requests, you can raise an issue at the issue page. You can also contact Margherita Cavattoni at cavattonimargherita@gmail.com

## Citation
M. Cavattoni, M. Comin,
*"Improved Metagenomic Reads Classification with ClassGraph2"*,
under submission.
