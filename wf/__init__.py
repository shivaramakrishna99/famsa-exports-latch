'''
Export guide trees in the Newick format and distance matrices as CSV using FAMSA
'''
import subprocess
from pathlib import Path

from latch import large_task, small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchAuthor, LatchDir, LatchFile, LatchMetadata, LatchParameter
from enum import Enum
from typing import Union, Optional, List
from latch.resources.launch_plan import LaunchPlan

class GuideTree(Enum):
    sl = 'Single Linkage Tree' # Single Linkage Tree
    upgma = 'UPGMA' # UPGMA
    nj = 'Neighbor Joining Tree' # Neighbor Joining Tree

class DistExport(Enum):
    dt = 'Distance Matrix' # Distance
    pid = 'Pairwise Identity Matrix' # Pairwise Identity

@large_task
def famsa_export_task(
    input: LatchFile,
    gtOutput: str = 'custom_guide_tree',
    gtExport: GuideTree = GuideTree.sl,
    ) -> LatchFile:

    input_path = Path(input).resolve()

    gtOutput+='.dnd'

    gtOptions = {
        'Single Linkage Tree':'sl',
        'UPGMA':'upgma',
        'Neighbor Joining Tree':'nj'
        }

    gtValue = str(gtExport.value)
    for gt in gtOptions.keys():
        if gt == gtValue:
            gtValue = gtOptions[gt]

    gt_export_cmd = [
        './FAMSA/famsa',
        '-gt',
        gtValue,
        '-gt_export',
        str(input_path),
        str(gtOutput)
    ]

    subprocess.run(gt_export_cmd)

    return LatchFile(str(gtOutput), f"latch:///FAMSA-Exports/{gtOutput}")

"""The metadata included here will be injected into your interface."""

metadata = LatchMetadata(
    display_name="FAMSA - Guide Tree Export",
    documentation="https://github.com/refresh-bio/FAMSA/blob/master/README.md",
    author=LatchAuthor(
        name="Shivaramakrishna Srinivasan",
        email="shivaramakrishna.srinivasan@gmail.com",
        github="github.com/shivaramakrishna99",
    ),
    repository="https://github.com/shivaramakrishna99/famsa-latch",
    license="GPL-3.0",
    parameters= {
        "input": LatchParameter(
            display_name="Input",
            description="Upload a FASTA file",
        ),
        "gtExport": LatchParameter(
            display_name="Tree Algorithm Type",
            description="Choose between three guide trees to export your file in the Newick format",
            section_title="Guide Tree Export"
        ),
        "gtOutput": LatchParameter(
            display_name="Guide Tree Output",
        )
    },
)

@workflow(metadata)
def famsa_export(
    input: LatchFile,
    gtOutput: str = 'custom_guide_tree',
    gtExport: GuideTree = GuideTree.sl
    ) -> LatchFile:
    """A progressive algorithm for large-scale multiple sequence alignments and guide tree generation
# **FAMSA -** Guide Tree Export
---

[GitHub Repository](https://github.com/shivaramakrishna99/famsa-exports-latch) | [Paper](https://www.nature.com/articles/srep33964) | [Source Documentation](https://github.com/refresh-bio/FAMSA/blob/master/README.md)

FAMSA (Fast and Accurate Multiple Sequence Alignment of huge protein families) is an algorithm for ultra-scale multiple sequence alignments (3M protein sequences in 5 minutes and 24 GB of RAM).

## **How to Use**

### **Parameters**

* `Input` - Upload an input FASTA file
* `Guide Tree` - Choose a guide tree method to generate your custom guide tree in the [Newick format](https://evolution.genetics.washington.edu/phylip/newicktree.html) (`.dnd`).  The available guide tree options are:
    * Single Linkage Tree (*default*)
    * UPGMA
    * Neighbour Joining Tree
* `Output` - Provide a name for your guide tree file.

You can import the generated guide tree in the [FAMSA Multiple Sequence Alignment workflow](https://console.latch.bio/workflows/67882/info) on Latch.

### **Test Data**
Select `Use Test Data` `>` `Refresh Bio's Data` to make use of a sample input to generate a guide tree

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
    """
    return famsa_export_task(
    input=input,
    gtOutput=gtOutput,
    gtExport=gtExport
)


"""
Add test data with a LaunchPlan. Provide default values in a dictionary with
the parameter names as the keys. These default values will be available under
the 'Test Data' dropdown at console.latch.bio.
"""
LaunchPlan(
    famsa_export,
    "Refresh Bio's Data",
    {
        "input": LatchFile("s3://latch-public/test-data/3701/FAMSA/test/adeno_fiber/adeno_fiber"),
    },
)
