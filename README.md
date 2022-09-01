**FAMSA -** Multiple Sequence Alignment | [LatchBio](https://console.latch.bio/)
---

[Workflow Page](https://console.latch.bio/explore/67306/info) | [Paper](https://www.nature.com/articles/srep33964) | [Source Documentation](https://github.com/refresh-bio/FAMSA/blob/master/README.md)

FAMSA (Fast and Accurate Multiple Sequence Alignment of huge protein families) is an algorithm for ultra-scale multiple sequence alignments (3M protein sequences in 5 minutes and 24 GB of RAM).

## **How to Use**

### **Parameters**

* `Input` - Upload an input FASTA file
* `Guide Tree` - Choose an existing guide tree method or upload your own file in Newick (`.dnd`) format.  The available guide tree options are:
    * Single Linkage Tree (*default*)
    * UPGMA
    * Neighbour Joining Tree
* `Medoid Tree` - Use MedoidTree heuristic for speeding up tree construction (*default: disabled*). If enabled, enter the threshold number of sequences (`n_seqs`) to use medoid trees only for sets with `n_seqs` or more
* `Output` - Provide a name for the output.  By default, the output is saved as a `.aln` file
* `Compress Output` - Enable gzipped output (default: disabled).  When enabled, compression level is set to `7`.  See the original documentation for this

### **Test Data**
Select `Use Test Data` `>` `Refresh Bio's Data` to make use of a sample input for alignment and a custom guide tree file

### **Original Documentation**
Check out the original [documentation](https://github.com/refresh-bio/FAMSA/blob/master/README.md) for FAMSA

## **Algorithms**
The major algorithmic features in FAMSA are:
- Pairwise distances based on the longest common subsequence (LCS). Thanks to the bit-level parallelism and utilization of SIMD extensions, LCS can be computed very fast.
- Single-linkage guide trees. While being very accurate, single-linkage trees can be established without storing entire distance matrix, which makes them suitable for large alignments. Although, the alternative guide tree algorithms like UPGMA and neighbour joining are also provided.
- The new heuristic based on K-Medoid clustering for generating fast guide trees. Medoid trees can be calculated in O(N logN) time and work with all types of subtrees (single linkage, UPGMA, NJ). The heuristic can be enabled with -medoidtree switch and allow aligning millions of sequences in minutes.

## **Results -** How does FAMSA compare to other algorithms?

The analysis was performed on our extHomFam 2 benchmark produced by combining Homstrad (March 2020) references with Pfam 33.1 families (NCBI variant). The data set was deposited at Zenodo: [https://zenodo.org/record/6524237](https://zenodo.org/record/6524237). The following algorithms were investigated:

| Name  | Version  | Command line  |
|---|---|---|
| Clustal&Omega;  | 1.2.4 |  `clustalo --threads=32 -i <input> -o <output>` |
| Clustal&Omega; iter2  | 1.2.4   | `clustalo --threads=32 --iter 2 -i <input> -o <output>` |
| MAFFT PartTree  |  7.453 | `mafft --thread 32 --anysymbol --quiet --parttree <input> -o <output>` |
| MAFFT DPPartTree  |  7.453 |  `mafft --thread 32 --anysymbol --quiet --dpparttree <input> -o <output>` |
| Kalign3 | 3.3.2 | `kalign -i <input> -o <output>` | 
| FAMSA  | 1.6.2  | `famsa -t 32 <input> <output>`  |
| FAMSA 2 | 2.0.1  | `famsa -t 32 -gz <input> <output>`  |
| FAMSA 2 Medoid | 2.0.1  | `famsa -t 32 -medoidtree -gt upgma -gz <input> <output>`  |


The tests were performed with 32 computing threads on a machine with AMD Ryzen Threadripper 3990X CPU and 256 GB of RAM. For each extHomFam 2 subset we measured a fraction of properly aligned columns (TC score) as well as a total running time and a maximum memory requirements. The results are presented in the figure below. Notches at boxplots indicate 95% confidence interval for median, triangle represent means. The missing series for some algorithm-set pairs indicate that the running times exceeded a week. Kalign3 failed to process 10 families (5 in second, 3 in fourth, and 2 in the largest subset). FAMSA 2 alignments were stored in gzip format (`-gz` switch). 

![extHomFam-v2-TC-comparison](https://user-images.githubusercontent.com/14868954/171652224-af88d980-5b49-4dcc-95e7-4de5dc152fb3.png)


The most important observations are as follows: 
* FAMSA 2 was superior in terms of accuracy to all the competitors. Only on the smallest families (*N* < 10k) Clustal&Omega; kept up with our algorithm.
* The advantage of FAMSA 2 increased with the number of sequences and reached 20-30 percent points for (100k, 250k] subset. 
* FAMSA 2 with medoid trees offered astonishing throughput (a familiy PF00005 of 3 million ABC transporters was aligned in 5 minutes) with accuracy only slightly inferior to that of the default single linkage trees.
* None of the competing algorithms was able to complete all the families in the largest [250k, 3M) subset.
* The memory requirements of FAMSA 2 allow ultra-scale analyzes at a desktop computer (24 GB for 3M sequences).

## **Citation**

[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. 
Scientific Reports, 6, 33964](https://www.nature.com/articles/srep33964)

---

**This documentation is derived from the [source documentation](https://github.com/refresh-bio/FAMSA/blob/master/README.md)**

**This workflow is authored by Shivaramakrishna Srinivasan.**
**Feel free to reach out to me via [email](mailto:shivaramakrishna.srinivasan@gmail.com) regarding any suggestions/feedback to improve this workflow.**