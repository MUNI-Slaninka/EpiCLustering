# EpiCLustering
Differentiating between identical nucleotide sequences using epigenetic information 

## Description
This program clusters reads based on their epigenetic information(methylation status of CG sites).

## Instalation
Clone or download and extract this GitHub repository from https://github.com/MUNI-Slaninka/EpiCLustering.

Install dependencies pandas, numpy, scikit-learn, pysam.
This can beachieved in terminal with pip in following manner:
```bash
pip install pandas numpy scikit-learn pysam
```
This can be achieved in terminal with conda in following manner:
```bash
conda install -c anaconda numpy pandas scikit-learn
conda install -c bioconda pysam
```

## Usage
Run the main.py file inside the the downloaded repository from terminal with the use of python:
```bash
python3 main.py [arguments]
```
or
```bash
python main.py [arguments]
```
For help with arguments usage open readme.md file or run:
```bash
python main.py --help
```

### Input
Bam files that should be mapped and have Mm and M1 tags filled
### Output
Bam files with single reads stored in direstories based on their clusters
Optional: CSV from pandas dataframe of C and CG positions used for clustering
### Arguments
```bash
-i INPUT, --input INPUT
```
Specify path to input BAM file

```bash
-o OUT_DIR, --out-dir OUT_DIR
```
Specify path to directory in which output will be stored

```bash
-n REF_NAME, --ref-name REF_NAME
```
Specify name of reference

```bash
-s START, --start START
```
Specify starting position

```bash
-e END, --end END
```
Specify end position

```bash
-t CLUST_TYPES, --clust-types CLUST_TYPES
```
Specify types of clustering algorithms (default: ABDK)

```bash
-d IMPUTER_DIVISOR, --imputer-divisor IMPUTER_DIVISOR
```
Specify types of clustering algorithms (default: 4)

```bash
-c CLUSTERS, --clusters CLUSTERS
```
Specify number of expected clusters (for non-Density based clustering) (default: 2)

```bash
-m MIN_SAMPLES, --min-samples MIN_SAMPLES
```
Specify number of minimal samples in cluster (for Density based clustering) (default: None)

```bash
--min-divisor MIN_DIVISOR
```
Specify divisor for automatic minimum samples calculation (for Density based clustering) (default: 5)

```bash
--eps-divisor EPS_DIVISOR
```
Specify divisor for automatic epsilon value calculation (for Density based clustering) (default: 4)

```bash
--eps-step EPS_STEP
```
Specify step size for automatic epsilon value calculation (for Density based clustering) (default: 1)

```bash
--eps-range EPS_RANGE
```
Specify broadness of epsilon range for automatic epsilon value calculation (for Density based clustering) (default: 10)

```bash
--neg-weight NEG_WEIGHT
```
Specify weight of outliers for choice of the best Density based clustering (default: 3)

```bash
--out-bam OUT_BAM
```
Specify if clustering results will be output as bam file (default: True)

```bash
--out-positions OUT_POSITIONS
```
Specify if positions used for clustering will be output as csv file (default: False)
