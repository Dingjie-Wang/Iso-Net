# Manual for construction of co-expression networks using Iso-Net
This tutorial will help you infer co-expression networks using Iso-Net based on exon-level expression data. If you experience any problems following these steps, please don't hesitate to contact Dingjie.Wang@osumc.edu. <br>

# 1. Requirements
(1) System Requirements: Windows or Linux operating system, and R 3.6.1 or later <br>
(2) R package required: CompQuadForm and survey <br>

# 2. Download lastest version of Iso-Net
In order to run Iso-Net, you should first fetch the source code from Iso-Net git repository. <br>
$ git clone --recursive https://github.com/Dingjie-Wang/Iso-Net.git <br>
$ cd Iso-Net <br>

# 3. Run Iso-Net
Open R environment in Windows or Linux operating system, type the following commands: <br>
$ source("IsoNet.R"); <br>
$ PerformIsoNet("Input_file.txt","output_file.txt",sample_size); <br>
Here, "Input_file.txt" is the name of input file, "output_file.txt" is the name of result file and sample_size if the sample size of data.

# 4. The format of input file and the output file
(1) Input Format: An example of data file (Input_file.txt) is given in the package. In this file, the first column is the gene ID, the second column is the exon ID, and the expression of samples are from the third column to the last column. <br>
(2) Output Format: An example of result file (out_putfile.txt) is given in the package. In this file, each row represents an possible edge. The first and second column are the two genes' IDs. For RVNet method, it includes three columns and the third column is the RV cofficient. For MINet method, it includes four columns, the third column is the test statistics of the possible edge and the last column is the p-value. <br>

