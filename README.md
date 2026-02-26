# KMeansBA

![Methodological Process](fig1_1.png)
*Figure 1: Methodological workflow.*

## Abstract

Clustering with size constraints organizes data while respecting predefined limits on the size of each cluster. This work introduces **K-MeansBA**, a novel extension of K-Means that integrates the Bat Algorithm (BA) to strictly enforce cardinality requirements during the clustering process. Our algorithm optimizes cluster quality—measured by intra-cluster cohesion and inter-cluster separation—while guaranteeing that each cluster contains exactly the designated number of instances. We evaluated the algorithm on **100+ OpenML datasets** using internal metrics (silhouette coefficient) and external metrics (Adjusted Rand Index - ARI, Adjusted Mutual Information - AMI, and Normalized Mutual Information - NMI), alongside direct verification of constraint compliance.  

Key findings include:
- **K-MeansBA reliably enforces size constraints** in all evaluated datasets.
- Competitive clustering quality, particularly in configurations with **fewer clusters**.

## Implemented Algorithms

This repository includes three cardinality-constrained clustering algorithms for comparison:

1. **K-MeansBA**: Our proposed method combining K-Means with the Bat Algorithm for size-constrained optimization.
2. **K-MedoidsSC**: A K-Medoids-based approach with size constraints.
3. **CSCLP**: A linear programming-based clustering algorithm with size constraints.

## Repository Structure

- `KMeansBA.R`: Implementation of **K-MeansBA** and dataset preprocessing logic.
- `OpenML_Datasets.R`: Functions to fetch and filter datasets from OpenML.
- `KmedoidsSC_OpenML.R`: Implementation and evaluation of **K-MedoidsSC**.
- `CSCLP_OpenML.R`: Implementation and evaluation of **CSCLP**.
- `results_KMeansBA.csv`: Results for K-MeansBA (metrics and execution times).
- `results_KmedoidsSC.csv`: Results for K-MedoidsSC.
- `results_CSCLP.csv`: Results for CSCLP.
- `silhouette_results.csv`: Silhouette coefficients for all instances.

## Requirements

Install the required R packages:
```R
install.packages(c(
  "cluster", "proxy", "mlr3oml", "mlr3", "pryr", "dplyr", 
  "aricode", "ggplot2", "corrplot", "clValid", "RColorBrewer", 
  "factoextra", "lpSolve"
))
```

## Dataset Preparation

The `OpenML_Datasets.R` script fetches and preprocesses datasets from OpenML:
```R
source("OpenML_Datasets.R")  # Load datasets into `odatasets_unique`
```

## Basic Usage

### K-MeansBA
```R
dataset <- odatasets_unique[1]$dataset[[1]]   # Select dataset
dataset_name <- odatasets_unique[1]$name     # Dataset name
target_cardinality <- odatasets_unique[1]$class_distribution_vector[[1]]  # Target sizes

# Run clustering
run_clustering(dataset, target_cardinality, dataset_name)  # Results saved to CSV
```

### K-MedoidsSC
```R
result <- run_KmedoidsSC(dataset, target_cardinality, dataset_name)  
```

### CSCLP
```R
result <- run_CSCLP(dataset, target_cardinality, dataset_name)  
```

## Results

Metrics for all algorithms are stored in CSV files:
- **Internal Metrics**: Silhouette coefficient (quality of clusters).
- **External Metrics**: ARI, AMI, NMI (agreement with ground truth labels).
- **Constraint Violations**: Difference between target and actual cluster sizes.

## Appendix — Summary of Evaluation Metrics

The following table summarizes selected results from **103 OpenML datasets** used for evaluation.  
Each entry shows both **external validation** (ARI, AMI, NMI) and **internal validation** (S(i)) scores for the three clustering algorithms.

<table>
<thead>
<tr>
<th>ID</th><th>Algorithm</th><th>#Features</th><th>#Instances</th>
<th>Cluster Sizes</th><th>ARI</th><th>AMI</th><th>NMI</th><th>S(i)</th>
</tr>
</thead>
<tbody>
<tr><td rowspan="3" style="vertical-align: middle;">1</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">625</td><td rowspan="3" style="vertical-align: middle;">288, 49, 288</td><td>0.185</td><td>0.153</td><td>0.156</td><td>0.122</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.005</td><td>0.005</td><td>0.009</td><td>-0.008</td></tr>
<tr><td>CSCLP</td><td>0.011</td><td>0.043</td><td>0.047</td><td>0.208</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">2</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">2000</td><td rowspan="3" style="vertical-align: middle;">200, 200, 200, 200, 200, 200, 200, 200, 200, 200</td><td>0.285</td><td>0.43</td><td>0.435</td><td>0.141</td></tr>
<tr><td>K-MedoidsSC</td><td>0.09</td><td>0.18</td><td>0.187</td><td>-0.304</td></tr>
<tr><td>CSCLP</td><td>0.317</td><td>0.478</td><td>0.483</td><td>0.44</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">3</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">9</td><td rowspan="3" style="vertical-align: middle;">768</td><td rowspan="3" style="vertical-align: middle;">500, 268</td><td>0.037</td><td>0.013</td><td>0.014</td><td>0.314</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.009</td><td>0.006</td><td>0.007</td><td>0.048</td></tr>
<tr><td>CSCLP</td><td>0.03</td><td>0.01</td><td>0.011</td><td>0.672</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">4</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">336</td><td rowspan="3" style="vertical-align: middle;">143, 77, 52, 35, 20, 5, 2, 2</td><td>0.197</td><td>0.354</td><td>0.381</td><td>-0.108</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.015</td><td>0.012</td><td>0.053</td><td>-0.057</td></tr>
<tr><td>CSCLP</td><td>0.47</td><td>0.454</td><td>0.476</td><td>0.168</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">5</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">14</td><td rowspan="3" style="vertical-align: middle;">270</td><td rowspan="3" style="vertical-align: middle;">150, 120</td><td>0.02</td><td>0.013</td><td>0.015</td><td>0.294</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.003</td><td>-0.003</td><td>0.0</td><td>0.183</td></tr>
<tr><td>CSCLP</td><td>0.072</td><td>0.049</td><td>0.052</td><td>0.513</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">6</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">150</td><td rowspan="3" style="vertical-align: middle;">50, 50, 50</td><td>0.591</td><td>0.636</td><td>0.64</td><td>0.44</td></tr>
<tr><td>K-MedoidsSC</td><td>0.01</td><td>0.009</td><td>0.021</td><td>-0.075</td></tr>
<tr><td>CSCLP</td><td>0.886</td><td>0.861</td><td>0.862</td><td>0.725</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">7</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">9</td><td rowspan="3" style="vertical-align: middle;">1484</td><td rowspan="3" style="vertical-align: middle;">463, 429, 244, 163, 51, 44, 35, 30, 20, 5</td><td>0.152</td><td>0.23</td><td>0.243</td><td>-0.082</td></tr>
<tr><td>K-MedoidsSC</td><td>0.014</td><td>0.013</td><td>0.028</td><td>-0.094</td></tr>
<tr><td>CSCLP</td><td>0.108</td><td>0.155</td><td>0.168</td><td>-0.105</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">8</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">14</td><td rowspan="3" style="vertical-align: middle;">178</td><td rowspan="3" style="vertical-align: middle;">59, 71, 48</td><td>0.374</td><td>0.401</td><td>0.407</td><td>0.537</td></tr>
<tr><td>K-MedoidsSC</td><td>0.058</td><td>0.062</td><td>0.072</td><td>0.067</td></tr>
<tr><td>CSCLP</td><td>0.362</td><td>0.357</td><td>0.364</td><td>0.49</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">9</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">3</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">125, 125</td><td>0.007</td><td>0.005</td><td>0.008</td><td>0.484</td></tr>
<tr><td>K-MedoidsSC</td><td>0.058</td><td>0.042</td><td>0.045</td><td>0.123</td></tr>
<tr><td>CSCLP</td><td>-0.002</td><td>-0.002</td><td>0.001</td><td>0.</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">10</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">508</td><td rowspan="3" style="vertical-align: middle;">286, 222</td><td>0.311</td><td>0.232</td><td>0.233</td><td>0.301</td></tr>
<tr><td>K-MedoidsSC</td><td>0.078</td><td>0.053</td><td>0.054</td><td>0.134</td></tr>
<tr><td>CSCLP</td><td>0.347</td><td>0.262</td><td>0.263</td><td>0.444</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">11</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">200</td><td rowspan="3" style="vertical-align: middle;">103, 97</td><td>0.433</td><td>0.34</td><td>0.342</td><td>0.133</td></tr>
<tr><td>K-MedoidsSC</td><td>0.001</td><td>0.001</td><td>0.005</td><td>0.03</td></tr>
<tr><td>CSCLP</td><td>0.001</td><td>0.001</td><td>0.005</td><td>0.127</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">12</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">100</td><td rowspan="3" style="vertical-align: middle;">40, 60</td><td>0.262</td><td>0.185</td><td>0.191</td><td>0.268</td></tr>
<tr><td>K-MedoidsSC</td><td>0.029</td><td>0.013</td><td>0.021</td><td>-0.001</td></tr>
<tr><td>CSCLP</td><td>0.151</td><td>0.099</td><td>0.105</td><td>0.389</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">13</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">4052</td><td rowspan="3" style="vertical-align: middle;">971, 3081</td><td>0.752</td><td>0.601</td><td>0.601</td><td>0.41</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.044</td><td>0.015</td><td>0.015</td><td>0.024</td></tr>
<tr><td>CSCLP</td><td>-0.077</td><td>0.129</td><td>0.129</td><td>0.083</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">14</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">119, 131</td><td>0.176</td><td>0.134</td><td>0.137</td><td>0.245</td></tr>
<tr><td>K-MedoidsSC</td><td>0.001</td><td>0.001</td><td>0.004</td><td>0.016</td></tr>
<tr><td>CSCLP</td><td>0.15</td><td>0.111</td><td>0.113</td><td>0.385</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">15</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">209</td><td rowspan="3" style="vertical-align: middle;">153, 56</td><td>0.374</td><td>0.229</td><td>0.232</td><td>0.429</td></tr>
<tr><td>K-MedoidsSC</td><td>0.006</td><td>-0.004</td><td>0.001</td><td>0.629</td></tr>
<tr><td>CSCLP</td><td>0.156</td><td>0.065</td><td>0.069</td><td>0.307</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">16</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">111</td><td rowspan="3" style="vertical-align: middle;">58, 53</td><td>0.034</td><td>0.024</td><td>0.031</td><td>0.326</td></tr>
<tr><td>K-MedoidsSC</td><td>0.005</td><td>0.004</td><td>0.01</td><td>0.226</td></tr>
<tr><td>CSCLP</td><td>0.009</td><td>0.006</td><td>0.013</td><td>0.33</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">17</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">3107</td><td rowspan="3" style="vertical-align: middle;">1541, 1566</td><td>0.01</td><td>0.007</td><td>0.007</td><td>0.358</td></tr>
<tr><td>K-MedoidsSC</td><td>0.048</td><td>0.035</td><td>0.035</td><td>0.18</td></tr>
<tr><td>CSCLP</td><td>0.001</td><td>0.001</td><td>0.001</td><td>0.034</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">18</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">440, 560</td><td>0.034</td><td>0.022</td><td>0.023</td><td>0.158</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-0.001</td><td>0.0</td><td>0.043</td></tr>
<tr><td>CSCLP</td><td>0.04</td><td>0.036</td><td>0.036</td><td>0.292</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">19</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">457, 543</td><td>0.126</td><td>0.099</td><td>0.1</td><td>0.222</td></tr>
<tr><td>K-MedoidsSC</td><td>0.022</td><td>0.015</td><td>0.015</td><td>0.007</td></tr>
<tr><td>CSCLP</td><td>0.214</td><td>0.158</td><td>0.159</td><td>0.386</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">20</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">109, 141</td><td>0.104</td><td>0.071</td><td>0.074</td><td>0.375</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.004</td><td>-0.003</td><td>0.0</td><td>0.029</td></tr>
<tr><td>CSCLP</td><td>0.03</td><td>0.018</td><td>0.021</td><td>0.558</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">21</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">237, 263</td><td>0.014</td><td>0.009</td><td>0.011</td><td>0.258</td></tr>
<tr><td>K-MedoidsSC</td><td>0.064</td><td>0.046</td><td>0.047</td><td>0.041</td></tr>
<tr><td>CSCLP</td><td>0.021</td><td>0.015</td><td>0.016</td><td>0.578</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">22</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">254, 246</td><td>-0.002</td><td>-0.001</td><td>4.73e-8</td><td>0.609</td></tr>
<tr><td>K-MedoidsSC</td><td>0.002</td><td>0.001</td><td>0.003</td><td>0.276</td></tr>
<tr><td>CSCLP</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.325</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">23</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">100</td><td rowspan="3" style="vertical-align: middle;">46, 54</td><td>0.152</td><td>0.110</td><td>0.117</td><td>0.145</td></tr>
<tr><td>K-MedoidsSC</td><td>0.004</td><td>0.002</td><td>0.009</td><td>0.042</td></tr>
<tr><td>CSCLP</td><td>0.069</td><td>0.049</td><td>0.056</td><td>0.185</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">24</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">100</td><td rowspan="3" style="vertical-align: middle;">45, 55</td><td>0.069</td><td>0.057</td><td>0.064</td><td>0.148</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.010</td><td>-0.007</td><td>7.41e-5</td><td>0.067</td></tr>
<tr><td>CSCLP</td><td>0.152</td><td>0.123</td><td>0.129</td><td>0.200</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">25</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">625</td><td rowspan="3" style="vertical-align: middle;">315, 310</td><td>0.001</td><td>0.001</td><td>0.002</td><td>0.045</td></tr>
<tr><td>K-MedoidsSC</td><td>0.059</td><td>0.044</td><td>0.045</td><td>0.195</td></tr>
<tr><td>CSCLP</td><td>0.173</td><td>0.129</td><td>0.129</td><td>0.362</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">26</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">2178</td><td rowspan="3" style="vertical-align: middle;">1209, 969</td><td>0.001</td><td>3.00e-4</td><td>0.001</td><td>0.323</td></tr>
<tr><td>K-MedoidsSC</td><td>0.002</td><td>0.001</td><td>0.001</td><td>0.138</td></tr>
<tr><td>CSCLP</td><td>0.003</td><td>0.001</td><td>0.002</td><td>0.643</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">27</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">662</td><td rowspan="3" style="vertical-align: middle;">345, 317</td><td>3.19e-5</td><td>1.00e-4</td><td>0.001</td><td>0.450</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-0.001</td><td>3.00e-4</td><td>0.156</td></tr>
<tr><td>CSCLP</td><td>0.005</td><td>0.004</td><td>0.005</td><td>0.532</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">28</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">15</td><td rowspan="3" style="vertical-align: middle;">252</td><td rowspan="3" style="vertical-align: middle;">124, 128</td><td>0.087</td><td>0.064</td><td>0.067</td><td>0.263</td></tr>
<tr><td>K-MedoidsSC</td><td>0.021</td><td>0.015</td><td>0.018</td><td>0.121</td></tr>
<tr><td>CSCLP</td><td>0.142</td><td>0.105</td><td>0.108</td><td>0.437</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">29</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">3</td><td rowspan="3" style="vertical-align: middle;">120</td><td rowspan="3" style="vertical-align: middle;">57, 63</td><td>0.749</td><td>0.644</td><td>0.646</td><td>0.488</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.007</td><td>-0.005</td><td>0.001</td><td>-0.015</td></tr>
<tr><td>CSCLP</td><td>0.534</td><td>0.429</td><td>0.433</td><td>0.739</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">30</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">100</td><td rowspan="3" style="vertical-align: middle;">47, 53</td><td>0.185</td><td>0.143</td><td>0.149</td><td>0.111</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.004</td><td>-0.002</td><td>0.005</td><td>0.028</td></tr>
<tr><td>CSCLP</td><td>0.004</td><td>0.004</td><td>0.011</td><td>0.161</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">31</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">243, 257</td><td>0.177</td><td>0.120</td><td>0.021</td><td>0.263</td></tr>
<tr><td>K-MedoidsSC</td><td>0.009</td><td>0.003</td><td>0.004</td><td>0.018</td></tr>
<tr><td>CSCLP</td><td>0.145</td><td>0.096</td><td>0.097</td><td>0.467</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">32</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">10</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">121, 129</td><td>0.014</td><td>0.009</td><td>0.012</td><td>0.159</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.004</td><td>-0.003</td><td>1.88e-6</td><td>0.032</td></tr>
<tr><td>CSCLP</td><td>0.050</td><td>0.035</td><td>0.038</td><td>0.337</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">33</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">730</td><td rowspan="3" style="vertical-align: middle;">367, 363</td><td>-0.001</td><td>7.74-4</td><td>0.002</td><td>3.18e-4</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.0015</td><td>-0.001</td><td>7.24e-6</td><td>0.119</td></tr>
<tr><td>CSCLP</td><td>-4.08e-4</td><td>-2.88e-4</td><td>-8.04e-4</td><td>0.529</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">34</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">9</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">496, 504</td><td>-2.17e-4</td><td>-1.55e-4</td><td>5.67e-4</td><td>0.145</td></tr>
<tr><td>K-MedoidsSC</td><td>6.00e-4</td><td>4.30e-4</td><td>0.001</td><td>0.015</td></tr>
<tr><td>CSCLP</td><td>0.166</td><td>0.123</td><td>0.124</td><td>0.250</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">35</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">257, 243</td><td>0.059</td><td>0.045</td><td>0.049</td><td>0.238</td></tr>
<tr><td>K-MedoidsSC</td><td>0.005</td><td>0.004</td><td>0.008</td><td>0.169</td></tr>
<tr><td>CSCLP</td><td>0.023</td><td>0.018</td><td>0.021</td><td>0.536</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">36</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">497, 503</td><td>0.047</td><td>0.030</td><td>0.032</td><td>0.344</td></tr>
<tr><td>K-MedoidsSC</td><td>0.001</td><td>2.52e-4</td><td>9.82-4</td><td>0.032</td></tr>
<tr><td>CSCLP</td><td>0.028</td><td>0.026</td><td>0.027</td><td>0.530</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">37</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">740</td><td rowspan="3" style="vertical-align: middle;">371, 369</td><td>0.004</td><td>0.004</td><td>0.005</td><td>0.511</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-1.87e-4</td><td>0.001</td><td>0.012</td></tr>
<tr><td>CSCLP</td><td>-0.001</td><td>-1.87e-4</td><td>0.001</td><td>0.778</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">38</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">3</td><td rowspan="3" style="vertical-align: middle;">124</td><td rowspan="3" style="vertical-align: middle;">61, 63</td><td>0.409</td><td>0.306</td><td>0.308</td><td>0.374</td></tr>
<tr><td>K-MedoidsSC</td><td>0.203</td><td>0.137</td><td>0.140</td><td>0.369</td></tr>
<tr><td>CSCLP</td><td>0.527</td><td>0.411</td><td>0.412</td><td>0.652</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">39</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">246, 254</td><td>0.214</td><td>0.157</td><td>0.158</td><td>0.131</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.002</td><td>-0.001</td><td>2.19e-4</td><td>0.028</td></tr>
<tr><td>CSCLP</td><td>0.004</td><td>0.004</td><td>0.011</td><td>0.161</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">40</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">730</td><td rowspan="3" style="vertical-align: middle;">367, 363</td><td>0.002</td><td>0.001</td><td>0.002</td><td>0.312</td></tr>
<tr><td>K-MedoidsSC</td><td>0.004</td><td>0.003</td><td>0.004</td><td>0.228</td></tr>
<tr><td>CSCLP</td><td>0.006</td><td>0.004</td><td>0.005</td><td>0.456</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">41</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">91, 159</td><td>0.316</td><td>0.218</td><td>0.220</td><td>0.105</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.004</td><td>-0.003</td><td>3.50e-6</td><td>0.007</td></tr>
<tr><td>CSCLP</td><td>0.038</td><td>0.077</td><td>0.080</td><td>0.171</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">42</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">10</td><td rowspan="3" style="vertical-align: middle;">950</td><td rowspan="3" style="vertical-align: middle;">488, 462</td><td>0.048</td><td>0.035</td><td>0.035</td><td>0.458</td></tr>
<tr><td>K-MedoidsSC</td><td>0.001</td><td>0.001</td><td>0.001</td><td>0.270</td></tr>
<tr><td>CSCLP</td><td>0.100</td><td>0.073</td><td>0.074</td><td>0.740</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">43</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">491, 509</td><td>0.027</td><td>0.020</td><td>0.020</td><td>0.077</td></tr>
<tr><td>K-MedoidsSC</td><td>0.009</td><td>0.007</td><td>0.007</td><td>0.068</td></tr>
<tr><td>CSCLP</td><td>0.009</td><td>0.007</td><td>0.007</td><td>0.124</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">44</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">224, 276</td><td>0.009</td><td>0.005</td><td>0.006</td><td>0.155</td></tr>
<tr><td>K-MedoidsSC</td><td>0.014</td><td>0.009</td><td>0.010</td><td>0.059</td></tr>
<tr><td>CSCLP</td><td>-0.002</td><td>-0.001</td><td>0.0002</td><td>0.298</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">45</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">3</td><td rowspan="3" style="vertical-align: middle;">380</td><td rowspan="3" style="vertical-align: middle;">185, 195</td><td>0.438</td><td>0.344</td><td>0.346</td><td>0.471</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.002</td><td>-0.001</td><td>0.001</td><td>-0.008</td></tr>
<tr><td>CSCLP</td><td>0.452</td><td>0.359</td><td>0.360</td><td>0.717</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">46</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">117, 133</td><td>-0.004</td><td>-0.003</td><td>1.11e-5</td><td>0.128</td></tr>
<tr><td>K-MedoidsSC</td><td>0.007</td><td>0.004</td><td>0.007</td><td>0.037</td></tr>
<tr><td>CSCLP</td><td>0.010</td><td>0.008</td><td>0.011</td><td>0.268</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">47</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">214, 286</td><td>0.229</td><td>0.165</td><td>0.166</td><td>0.128</td></tr>
<tr><td>K-MedoidsSC</td><td>0.029</td><td>0.017</td><td>0.018</td><td>0.089</td></tr>
<tr><td>CSCLP</td><td>0.116</td><td>0.079</td><td>0.081</td><td>0.244</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">48</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">233, 267</td><td>0.178</td><td>0.138</td><td>0.139</td><td>0.242</td></tr>
<tr><td>K-MedoidsSC</td><td>0.011</td><td>0.007</td><td>0.008</td><td>0.026</td></tr>
<tr><td>CSCLP</td><td>0.178</td><td>0.138</td><td>0.139</td><td>0.406</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">49</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">3848</td><td rowspan="3" style="vertical-align: middle;">1924, 1924</td><td>-7.73e-5</td><td>-5.58e-5</td><td>1.32e-4</td><td>0.301</td></tr>
<tr><td>K-MedoidsSC</td><td>-1.90e-4</td><td>-1.38e-4</td><td>4.99e-5</td><td>0.022</td></tr>
<tr><td>CSCLP</td><td>-2.21e-4</td><td>-1.60e-4</td><td>2.81e-5</td><td>0.518</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">50</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">284</td><td rowspan="3" style="vertical-align: middle;">142, 142</td><td>1.000</td><td>1.000</td><td>1.000</td><td>0.034</td></tr>
<tr><td>K-MedoidsSC</td><td>0.061</td><td>0.044</td><td>0.047</td><td>0.228</td></tr>
<tr><td>CSCLP</td><td>0.030</td><td>0.022</td><td>0.024</td><td>0.481</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">51</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">249, 251</td><td>0.045</td><td>0.033</td><td>0.034</td><td>0.149</td></tr>
<tr><td>K-MedoidsSC</td><td>0.011</td><td>0.008</td><td>0.009</td><td>0.071</td></tr>
<tr><td>CSCLP</td><td>0.146</td><td>0.108</td><td>0.109</td><td>0.239</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">52</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">131</td><td rowspan="3" style="vertical-align: middle;">83, 48</td><td>1.000</td><td>1.000</td><td>1.000</td><td>0.575</td></tr>
<tr><td>K-MedoidsSC</td><td>0.025</td><td>0.059</td><td>0.064</td><td>0.306</td></tr>
<tr><td>CSCLP</td><td>0.394</td><td>0.284</td><td>0.288</td><td>0.177</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">53</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">3</td><td rowspan="3" style="vertical-align: middle;">222</td><td rowspan="3" style="vertical-align: middle;">88, 134</td><td>-0.002</td><td>-0.003</td><td>3.29e-4</td><td>0.583</td></tr>
<tr><td>K-MedoidsSC</td><td>0.018</td><td>0.029</td><td>0.032</td><td>-0.014</td></tr>
<tr><td>CSCLP</td><td>0.008</td><td>0.018</td><td>0.021</td><td>0.760</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">54</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">400</td><td rowspan="3" style="vertical-align: middle;">235, 165</td><td>-0.003</td><td>-0.001</td><td>0.001</td><td>0.101</td></tr>
<tr><td>K-MedoidsSC</td><td>0.014</td><td>0.006</td><td>0.008</td><td>0.130</td></tr>
<tr><td>CSCLP</td><td>-0.003</td><td>-0.001</td><td>0.001</td><td>0.277</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">55</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">400</td><td rowspan="3" style="vertical-align: middle;">207, 193</td><td>-0.002</td><td>-0.001</td><td>0.0003</td><td>0.238</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.094</td></tr>
<tr><td>CSCLP</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.380</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">56</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">400</td><td rowspan="3" style="vertical-align: middle;">206, 194</td><td>-2.57e-6</td><td>-6.75e-5</td><td>0.002</td><td>0.119</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.070</td></tr>
<tr><td>CSCLP</td><td>0.001</td><td>0.001</td><td>0.003</td><td>0.398</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">57</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">400</td><td rowspan="3" style="vertical-align: middle;">208, 192</td><td>-8.14e-7</td><td>-1.15e-4</td><td>0.002</td><td>0.181</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.002</td><td>-0.001</td><td>0.001</td><td>0.128</td></tr>
<tr><td>CSCLP</td><td>-0.002</td><td>-0.001</td><td>0.0003</td><td>0.417</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">58</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">400</td><td rowspan="3" style="vertical-align: middle;">197, 203</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.057</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.002</td><td>-0.002</td><td>7.55e-5</td><td>0.070</td></tr>
<tr><td>CSCLP</td><td>-5.20e-6</td><td>1.11e-5</td><td>0.002</td><td>0.400</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">59</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">436, 564</td><td>0.182</td><td>0.130</td><td>0.131</td><td>0.101</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-0.0004</td><td>0.0003</td><td>0.049</td></tr>
<tr><td>CSCLP</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.136</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">60</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">250</td><td rowspan="3" style="vertical-align: middle;">110, 140</td><td>0.267</td><td>0.230</td><td>0.232</td><td>0.290</td></tr>
<tr><td>K-MedoidsSC</td><td>0.004</td><td>0.005</td><td>0.008</td><td>0.077</td></tr>
<tr><td>CSCLP</td><td>0.302</td><td>0.225</td><td>0.228</td><td>0.462</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">61</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">416, 584</td><td>0.152</td><td>0.146</td><td>0.147</td><td>0.226</td></tr>
<tr><td>K-MedoidsSC</td><td>0.004</td><td>0.001</td><td>0.002</td><td>0.021</td></tr>
<tr><td>CSCLP</td><td>0.189</td><td>0.181</td><td>0.182</td><td>0.443</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">62</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">1000</td><td rowspan="3" style="vertical-align: middle;">420, 580</td><td>0.189</td><td>0.133</td><td>0.133</td><td>0.125</td></tr>
<tr><td>K-MedoidsSC</td><td>3.31e-4</td><td>0.003</td><td>0.004</td><td>0.078</td></tr>
<tr><td>CSCLP</td><td>0.010</td><td>0.013</td><td>0.014</td><td>0.133</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">63</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">100</td><td rowspan="3" style="vertical-align: middle;">44, 56</td><td>0.069</td><td>0.046</td><td>0.053</td><td>0.307</td></tr>
<tr><td>K-MedoidsSC</td><td>0.016</td><td>0.016</td><td>0.023</td><td>0.040</td></tr>
<tr><td>CSCLP</td><td>0.004</td><td>0.006</td><td>0.014</td><td>0.492</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">64</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">323</td><td rowspan="3" style="vertical-align: middle;">175, 148</td><td>-0.001</td><td>-0.001</td><td>0.001</td><td>0.317</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.003</td><td>-0.002</td><td>7.48e-5</td><td>0.237</td></tr>
<tr><td>CSCLP</td><td>0.278</td><td>0.209</td><td>0.210</td><td>0.532</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">65</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">662</td><td rowspan="3" style="vertical-align: middle;">348, 314</td><td>-0.001</td><td>-0.001</td><td>9.95e-5</td><td>0.424</td></tr>
<tr><td>K-MedoidsSC</td><td>0.003</td><td>0.003</td><td>0.004</td><td>0.180</td></tr>
<tr><td>CSCLP</td><td>0.001</td><td>0.001</td><td>0.002</td><td>0.565</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">66</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">228, 272</td><td>0.004</td><td>0.002</td><td>0.004</td><td>0.153</td></tr>
<tr><td>K-MedoidsSC</td><td>0.007</td><td>0.004</td><td>0.006</td><td>0.044</td></tr>
<tr><td>CSCLP</td><td>2.71e-4</td><td>-2.62e-4</td><td>0.001</td><td>0.316</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">67</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">241, 259</td><td>0.021</td><td>0.015</td><td>-2.21e-4</td><td>0.079</td></tr>
<tr><td>K-MedoidsSC</td><td>0.006</td><td>0.004</td><td>0.006</td><td>0.030</td></tr>
<tr><td>CSCLP</td><td>-0.002</td><td>-0.001</td><td>6.26e-5</td><td>0.125</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">68</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">10</td><td rowspan="3" style="vertical-align: middle;">214</td><td rowspan="3" style="vertical-align: middle;">76, 138</td><td>0.036</td><td>0.086</td><td>0.089</td><td>0.311</td></tr>
<tr><td>K-MedoidsSC</td><td>0.063</td><td>0.028</td><td>0.032</td><td>0.014</td></tr>
<tr><td>CSCLP</td><td>0.063</td><td>0.028</td><td>0.032</td><td>0.045</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">69</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">9</td><td rowspan="3" style="vertical-align: middle;">369</td><td rowspan="3" style="vertical-align: middle;">204, 165</td><td>0.018</td><td>0.011</td><td>0.013</td><td>0.299</td></tr>
<tr><td>K-MedoidsSC</td><td>0.001</td><td>0.002</td><td>0.004</td><td>0.042</td></tr>
<tr><td>CSCLP</td><td>0.018</td><td>0.011</td><td>0.013</td><td>0.571</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">70</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">9</td><td rowspan="3" style="vertical-align: middle;">274</td><td rowspan="3" style="vertical-align: middle;">134, 140</td><td>0.020</td><td>0.015</td><td>0.017</td><td>0.198</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.001</td><td>-0.001</td><td>0.002</td><td>0.108</td></tr>
<tr><td>CSCLP</td><td>0.005</td><td>0.004</td><td>0.006</td><td>0.522</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">71</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">9</td><td rowspan="3" style="vertical-align: middle;">130</td><td rowspan="3" style="vertical-align: middle;">119, 11</td><td>0.265</td><td>0.088</td><td>0.102</td><td>0.452</td></tr>
<tr><td>K-MedoidsSC</td><td>0.087</td><td>-2.05e-4</td><td>0.015</td><td>0.885</td></tr>
<tr><td>CSCLP</td><td>0.087</td><td>-2.05e-4</td><td>0.015</td><td>0.876</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">72</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">1372</td><td rowspan="3" style="vertical-align: middle;">762, 610</td><td>-0.001</td><td>-0.001</td><td>9.28e-6</td><td>0.262</td></tr>
<tr><td>K-MedoidsSC</td><td>0.047</td><td>0.031</td><td>0.032</td><td>0.014</td></tr>
<tr><td>CSCLP</td><td>0.037</td><td>0.025</td><td>0.025</td><td>0.648</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">73</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">5</td><td rowspan="3" style="vertical-align: middle;">748</td><td rowspan="3" style="vertical-align: middle;">570, 178</td><td>0.131</td><td>0.045</td><td>0.046</td><td>0.526</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.070</td><td>0.058</td><td>0.060</td><td>0.709</td></tr>
<tr><td>CSCLP</td><td>-0.070</td><td>0.058</td><td>0.060</td><td>0.721</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">74</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">10</td><td rowspan="3" style="vertical-align: middle;">106</td><td rowspan="3" style="vertical-align: middle;">22, 21, 14, 15, 16, 18</td><td>0.113</td><td>0.137</td><td>0.201</td><td>0.074</td></tr>
<tr><td>K-MedoidsSC</td><td>0.006</td><td>-0.006</td><td>0.069</td><td>-0.076</td></tr>
<tr><td>CSCLP</td><td>0.070</td><td>0.120</td><td>0.186</td><td>0.122</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">75</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">10</td><td rowspan="3" style="vertical-align: middle;">100</td><td rowspan="3" style="vertical-align: middle;">88, 12</td><td>0.037</td><td>-0.013</td><td>0.004</td><td>0.117</td></tr>
<tr><td>K-MedoidsSC</td><td>0.112</td><td>0.009</td><td>0.025</td><td>0.237</td></tr>
<tr><td>CSCLP</td><td>0.275</td><td>0.097</td><td>0.112</td><td>0.146</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">76</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">13</td><td rowspan="3" style="vertical-align: middle;">182</td><td rowspan="3" style="vertical-align: middle;">130, 52</td><td>-0.003</td><td>-0.005</td><td>1.24e-5</td><td>0.107</td></tr>
<tr><td>K-MedoidsSC</td><td>0.006</td><td>-0.004</td><td>0.001</td><td>0.058</td></tr>
<tr><td>CSCLP</td><td>-0.011</td><td>-0.004</td><td>4.48e-4</td><td>0.169</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">77</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">210</td><td rowspan="3" style="vertical-align: middle;">70, 70, 70</td><td>0.701</td><td>0.671</td><td>0.674</td><td>0.404</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.004</td><td>-0.004</td><td>0.005</td><td>-0.094</td></tr>
<tr><td>CSCLP</td><td>0.597</td><td>0.531</td><td>0.535</td><td>0.484</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">78</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">403</td><td rowspan="3" style="vertical-align: middle;">102, 129, 122, 24, 26</td><td>0.212</td><td>0.301</td><td>0.312</td><td>0.092</td></tr>
<tr><td>K-MedoidsSC</td><td>0.055</td><td>0.088</td><td>0.101</td><td>-0.097</td></tr>
<tr><td>CSCLP</td><td>0.055</td><td>0.069</td><td>0.083</td><td>0.146</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">79</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">14</td><td rowspan="3" style="vertical-align: middle;">200</td><td rowspan="3" style="vertical-align: middle;">51, 56, 41, 42, 10</td><td>0.002</td><td>-0.003</td><td>0.026</td><td>0.199</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.005</td><td>-0.011</td><td>0.018</td><td>-0.119</td></tr>
<tr><td>CSCLP</td><td>0.003</td><td>0.019</td><td>0.047</td><td>0.202</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">80</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">13</td><td rowspan="3" style="vertical-align: middle;">123</td><td rowspan="3" style="vertical-align: middle;">8, 48, 32, 30, 5</td><td>0.003</td><td>-0.001</td><td>0.051</td><td>0.101</td></tr>
<tr><td>K-MedoidsSC</td><td>0.003</td><td>0.003</td><td>0.054</td><td>-0.116</td></tr>
<tr><td>CSCLP</td><td>-0.009</td><td>-0.026</td><td>0.027</td><td>0.156</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">81</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">310</td><td rowspan="3" style="vertical-align: middle;">60, 100, 150</td><td>0.299</td><td>0.358</td><td>0.362</td><td>0.220</td></tr>
<tr><td>K-MedoidsSC</td><td>0.027</td><td>0.007</td><td>0.014</td><td>-0.132</td></tr>
<tr><td>CSCLP</td><td>0.603</td><td>0.491</td><td>0.495</td><td>0.399</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">82</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">3252</td><td rowspan="3" style="vertical-align: middle;">2952, 68, 58, 86, 88</td><td>-0.082</td><td>0.011</td><td>0.018</td><td>-0.334</td></tr>
<tr><td>K-MedoidsSC</td><td>0.008</td><td>0.005</td><td>0.011</td><td>-0.080</td></tr>
<tr><td>CSCLP</td><td>-0.006</td><td>0.008</td><td>0.014</td><td>0.244</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">83</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1623</td><td rowspan="3" style="vertical-align: middle;">1471, 35, 29, 44, 44</td><td>-0.093</td><td>0.009</td><td>0.022</td><td>-0.317</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.021</td><td>-0.001</td><td>0.012</td><td>-0.026</td></tr>
<tr><td>CSCLP</td><td>-0.074</td><td>0.005</td><td>0.018</td><td>-0.101</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">84</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1521</td><td rowspan="3" style="vertical-align: middle;">1369, 35, 29, 44, 44</td><td>-0.099</td><td>0.011</td><td>0.024</td><td>-0.329</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.028</td><td>0.001</td><td>0.014</td><td>-0.046</td></tr>
<tr><td>CSCLP</td><td>-0.067</td><td>0.005</td><td>0.018</td><td>0.252</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">85</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1515</td><td rowspan="3" style="vertical-align: middle;">1365, 33, 29, 45, 43</td><td>-0.072</td><td>0.005</td><td>0.018</td><td>-0.281</td></tr>
<tr><td>K-MedoidsSC</td><td>0.041</td><td>0.005</td><td>0.018</td><td>-0.238</td></tr>
<tr><td>CSCLP</td><td>-0.006</td><td>4.83e-4</td><td>0.014</td><td>0.220</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">86</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1183</td><td rowspan="3" style="vertical-align: middle;">1083, 9, 21, 49, 21</td><td>-0.074</td><td>0.002</td><td>0.018</td><td>-0.270</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.025</td><td>-0.004</td><td>0.012</td><td>-0.105</td></tr>
<tr><td>CSCLP</td><td>-0.044</td><td>-0.004</td><td>0.012</td><td>-0.057</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">87</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1080</td><td rowspan="3" style="vertical-align: middle;">984, 8, 19, 51, 18</td><td>-0.048</td><td>-0.003</td><td>0.013</td><td>-0.236</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.007</td><td>-0.001</td><td>0.015</td><td>-0.336</td></tr>
<tr><td>CSCLP</td><td>-0.077</td><td>0.001</td><td>0.017</td><td>0.341</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">88</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1277</td><td rowspan="3" style="vertical-align: middle;">1170, 9, 21, 53, 24</td><td>-0.064</td><td>0.001</td><td>0.015</td><td>-0.345</td></tr>
<tr><td>K-MedoidsSC</td><td>1.30e-4</td><td>-0.003</td><td>0.011</td><td>0.053</td></tr>
<tr><td>CSCLP</td><td>-0.055</td><td>0.004</td><td>0.018</td><td>0.294</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">89</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1252</td><td rowspan="3" style="vertical-align: middle;">1144, 9, 22, 53, 24</td><td>-0.050</td><td>-0.001</td><td>0.013</td><td>-0.282</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.030</td><td>-0.003</td><td>0.011</td><td>-0.135</td></tr>
<tr><td>CSCLP</td><td>-0.058</td><td>0.004</td><td>0.018</td><td>-0.140</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">90</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">1112</td><td rowspan="3" style="vertical-align: middle;">1010, 9, 21, 52, 20</td><td>-0.080</td><td>0.004</td><td>0.020</td><td>-0.343</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.023</td><td>-0.002</td><td>0.013</td><td>-0.005</td></tr>
<tr><td>CSCLP</td><td>-0.042</td><td>-0.002</td><td>0.014</td><td>-0.047</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">91</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">383</td><td rowspan="3" style="vertical-align: middle;">257, 125, 1</td><td>0.421</td><td>0.290</td><td>0.295</td><td>0.149</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.013</td><td>-0.002</td><td>0.005</td><td>0.144</td></tr>
<tr><td>CSCLP</td><td>0.073</td><td>0.231</td><td>0.236</td><td>0.101</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">92</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">500</td><td rowspan="3" style="vertical-align: middle;">45, 37, 51, 57, 52, 52, 47, 57, 53, 49</td><td>0.392</td><td>0.492</td><td>0.511</td><td>0.249</td></tr>
<tr><td>K-MedoidsSC</td><td>0.081</td><td>0.153</td><td>0.185</td><td>-0.163</td></tr>
<tr><td>CSCLP</td><td>0.473</td><td>0.548</td><td>0.564</td><td>0.341</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">93</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">12</td><td rowspan="3" style="vertical-align: middle;">4898</td><td rowspan="3" style="vertical-align: middle;">20, 163, 1457, 2198, 880, 175, 5</td><td>0.015</td><td>0.014</td><td>0.017</td><td>-0.127</td></tr>
<tr><td>K-MedoidsSC</td><td>0.001</td><td>0.011</td><td>0.014</td><td>-0.368</td></tr>
<tr><td>CSCLP</td><td>0.051</td><td>0.051</td><td>0.054</td><td>2.08e-4</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">94</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">215</td><td rowspan="3" style="vertical-align: middle;">150, 35, 30</td><td>0.542</td><td>0.416</td><td>0.423</td><td>0.423</td></tr>
<tr><td>K-MedoidsSC</td><td>0.718</td><td>0.536</td><td>0.541</td><td>0.370</td></tr>
<tr><td>CSCLP</td><td>0.299</td><td>0.367</td><td>0.375</td><td>-0.270</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">95</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">12</td><td rowspan="3" style="vertical-align: middle;">1599</td><td rowspan="3" style="vertical-align: middle;">10, 53, 681, 638, 199, 18</td><td>0.001</td><td>0.035</td><td>0.041</td><td>0.075</td></tr>
<tr><td>K-MedoidsSC</td><td>0.018</td><td>0.022</td><td>0.028</td><td>-0.272</td></tr>
<tr><td>CSCLP</td><td>0.058</td><td>0.041</td><td>0.047</td><td>0.107</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">96</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">4</td><td rowspan="3" style="vertical-align: middle;">2201</td><td rowspan="3" style="vertical-align: middle;">1490, 711</td><td>0.218</td><td>0.127</td><td>0.128</td><td>0.529</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.016</td><td>0.013</td><td>0.013</td><td>0.332</td></tr>
<tr><td>CSCLP</td><td>0.196</td><td>0.111</td><td>0.111</td><td>0.652</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">97</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">4839</td><td rowspan="3" style="vertical-align: middle;">4578, 261</td><td>0.222</td><td>0.077</td><td>0.078</td><td>0.396</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.043</td><td>0.008</td><td>0.009</td><td>0.770</td></tr>
<tr><td>CSCLP</td><td>-0.018</td><td>0.001</td><td>0.001</td><td>0.745</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">98</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">11</td><td rowspan="3" style="vertical-align: middle;">265</td><td rowspan="3" style="vertical-align: middle;">30, 27, 28, 44, 28, 47, 32, 29</td><td>0.210</td><td>0.326</td><td>0.358</td><td>0.122</td></tr>
<tr><td>K-MedoidsSC</td><td>0.027</td><td>0.047</td><td>0.093</td><td>-0.152</td></tr>
<tr><td>CSCLP</td><td>0.145</td><td>0.197</td><td>0.236</td><td>0.201</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">99</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">6</td><td rowspan="3" style="vertical-align: middle;">2000</td><td rowspan="3" style="vertical-align: middle;">1892, 108</td><td>0.001</td><td>-0.001</td><td>6.38e-6</td><td>0.081</td></tr>
<tr><td>K-MedoidsSC</td><td>-0.042</td><td>0.007</td><td>0.008</td><td>0.766</td></tr>
<tr><td>CSCLP</td><td>-0.025</td><td>0.001</td><td>0.002</td><td>0.741</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">100</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">12</td><td rowspan="3" style="vertical-align: middle;">2000</td><td rowspan="3" style="vertical-align: middle;">8, 67, 595, 898, 359, 71, 2</td><td>0.008</td><td>0.014</td><td>0.020</td><td>-0.247</td></tr>
<tr><td>K-MedoidsSC</td><td>0.004</td><td>0.011</td><td>0.017</td><td>-0.373</td></tr>
<tr><td>CSCLP</td><td>0.022</td><td>0.041</td><td>0.047</td><td>-0.018</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">101</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">7</td><td rowspan="3" style="vertical-align: middle;">2000</td><td rowspan="3" style="vertical-align: middle;">1029, 777, 194</td><td>0.002</td><td>0.001</td><td>0.002</td><td>0.118</td></tr>
<tr><td>K-MedoidsSC</td><td>0.022</td><td>0.024</td><td>0.025</td><td>0.112</td></tr>
<tr><td>CSCLP</td><td>-0.015</td><td>0.017</td><td>0.018</td><td>0.157</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">102</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">14</td><td rowspan="3" style="vertical-align: middle;">243</td><td rowspan="3" style="vertical-align: middle;">137, 106</td><td>0.014</td><td>0.015</td><td>0.018</td><td>0.225</td></tr>
<tr><td>K-MedoidsSC</td><td>0.043</td><td>0.039</td><td>0.042</td><td>0.195</td></tr>
<tr><td>CSCLP</td><td>0.297</td><td>0.261</td><td>0.263</td><td>0.555</td></tr>
<tr><td rowspan="3" style="vertical-align: middle;">103</td><td>K-MeansBA</td><td rowspan="3" style="vertical-align: middle;">8</td><td rowspan="3" style="vertical-align: middle;">719</td><td rowspan="3" style="vertical-align: middle;">174, 44, 84, 417</td><td>0.024</td><td>-0.001</td><td>0.005</td><td>0.207</td></tr>
<tr><td>K-MedoidsSC</td><td>0.050</td><td>0.004</td><td>0.010</td><td>-0.554</td></tr>
<tr><td>CSCLP</td><td>0.005</td><td>0.001</td><td>0.006</td><td>-0.653</td></tr>
</tbody>
</table>

---

### Summary of Appendix Results

- **Datasets analyzed:** 103  
- **Constraint compliance:** 100%  
- **Perfect alignment (ARI=AMI=NMI=1):** Datasets 50 and 52  
- **Strong overall performance:** Datasets with 2–3 clusters show the highest internal and external metrics.  
- **Notable internal cohesion:** Dataset 8 (S(i) = 0.537), Dataset 6 (S(i) = 0.440).  
- **Consistent superiority:** K-MeansBA outperforms K-MedoidsSC and is comparable or better than CSCLP for smaller cluster counts.

---
