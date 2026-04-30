# Character Tables of Symmetric Groups
## Contents
1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [Data](#data)
4. [Examples](#examples)
5. [Computations](#computations)
<!-- 5. Results
6. Conclusion -->

## Introduction
This is the repository for Lab of Geometry at Michigan (LoG(M): MATH 440) for the project Statistics of Character Tables of Symmetric Groups mentored by Prof. Sarah Peluse and Graduate Instructor Karthik Ganapathy. Our group members are Arnav Shah, Atharva Kulkarni, and Tony Zhang.

## Getting Started
Clone the repository on your laptop using the following command:
```
git clone https://github.com/TonyZhang2004/Character_Table_of_Symmetric_Groups.git
```
Change your working directory to the folder with our code as follows:
```
cd Character_Table_of_Symmetric_Groups
```
Install all pip requirements by running the following command:
```
pip install -r requirements.txt
```

Download the public data files from the Google Drive folder linked in the [Data](#data) section.

## Data
Full generated character-table data is too large to store directly in this GitHub repository. Public data files are available in this Google Drive folder:

[Public character-table data folder](https://drive.google.com/drive/folders/1A_DN63jodrUVwPwOMIMrhuG_VJAzheFu?usp=share_link)

For the full $S_{40}$ table, download the compressed CSV and place it in the repository root as `S40.csv.gz` or decompress it as `S40.csv`. The analysis notebook also checks `runs/S40/S40.csv.gz` and `runs/S40/S40.csv`.

To decompress the file while keeping the compressed copy:
```
gunzip -k S40.csv.gz
```

To verify the downloaded table without loading the full file into memory:
```
python3 verify_character_table.py --n 40 --csv-file S40.csv.gz
```

The verifier streams the CSV data, checks the first column against the hook-length formula, verifies the trivial and sign rows, checks dimension-square normalization, and tests selected row norms.

## Examples
The following examples demonstrate running some of our functions to copmute the character tables.
### Frobenius Formula ([frobenius.py](frobenius.py))
Import the function to compute the character table using the frobenius formula as follows:
```
import frobenius as fb
```
Note: The above command will only work if your file is in the same working directory as this repository. Please modify the above import statement according to your directory.

To compute the character table and print it, use the following code:
```
N = 6   # Size of Character Table
char_table = fb.get_character_table(N)
print(char_table)
```

### Murnghan-Nakayama Rule ([murnaghan_nakayama.py](murnaghan_nakayama.py))
Import the function to compute the character table using the frobenius formula as follows:
```
import murnaghan_nakayama as mn
```
Note: The above command will only work if your file is in the same working directory as this repository. Please modify the above import statement according to your directory.

To compute the character table and print it, you can use the following code:
```
N = 6   # Size of Character Table
char_table = mn.get_character_table(N)
print(char_table)
```

To compute the character values of, for instance, the column corresponding to the staircase partition of 6, you can use the following code:
```
N = 6   # Size of Character Table
col_bit_string = "010101"   # Abaci bit string representation of the column
char_table = mn.get_character_value_of_column(N, col_bit_string)
print(char_table)
```

For larger staircase columns, use the resumable streaming runner. The value of `N` must be triangular, i.e. `N = k(k + 1) / 2`.
```
python3 run_staircase_column.py --n 78 --output-dir runs/staircase
```

This writes a one-column CSV such as `runs/staircase/S78_staircase.csv` and logs progress without building the full table in memory.

For VM runs larger than 78, use:
```
bash run_staircase_on_vm.sh
```

This defaults to `K=13`, hence `N=91`. To run a larger triangular value, set `K`:
```
K=14 bash run_staircase_on_vm.sh  # N = 105
K=15 bash run_staircase_on_vm.sh  # N = 120
```

You can also save the character table in a file as follows:
```
N = 6   # Size of Character Table
character_table = mn.get_character_table(N, output_file_name="S6.csv")
```

Moreover, you can speed things up for calculating character values of large n if you have already calculated for smaller n by saving it in a memo. 
As an illustration, suppose you run the followinf code: 
```
N = 6   # Size of Character Table
character_table = mn.get_character_table(N, output_file_name="S6.csv", memo_file_name="memo.txt")
```
Now, if we want to calculate $S_{7}$, then we can speed it up by using the same memo as follows:
```
N = 7   # Size of Character Table
character_table = mn.get_character_table(N, output_file_name="S7.csv", memo_file_name="memo.txt")
```

### Heatmaps and Graphs ([char_table.py](char_table.py))
Please read the instructions given at the begining of each code block in our Jupyter notebook file named [char_table.py](char_table.py).

Download all required files that are mentioned right before the code block from the gdrive link provided [above](#getting-started).

You may also change some parameters as you require. The parameters which you may change are given at the begining of each code block. If you feel comfortable coding and want to modify our code to fit your requirements, you may change other parts of the code as well.

### Zero Density Analysis ([zero_density_analysis.ipynb](zero_density_analysis.ipynb))
The notebook `zero_density_analysis.ipynb` computes the density of zero entries in character tables of $S_n$. It computes small tables directly and streams the full $S_{40}$ CSV or compressed CSV file.

By default, the notebook looks for the full $S_{40}$ data at:
```
S40.csv.gz
S40.csv
runs/S40/S40.csv.gz
runs/S40/S40.csv
```

The output plot has $n$ on the x-axis and the density of zero entries on the y-axis.

## Computations
### [Frobenius Formula](https://en.wikipedia.org/wiki/Frobenius_formula)
Given an integer partition $\lambda = \lambda_1 + \lambda_2 + \cdots + \lambda_k$ of $n$, let $\chi^{\lambda}$ be the corresponding irreducible character of $S_n$, and let $\chi^{\lambda}_{\mu}$ be short for the value of $\chi^{\lambda}$ at any $g$ with cycle type $\mu$, denote $l_j  =\lambda_j + k - j$, and $i_j$ the number of times $j$ appears in $\mu$, so $\sum\limits_j i_jj = n$.

Then we have the following Frobenius's Formula: 
$\chi^{\lambda}_{\mu}$ =

$\text{coeff. of } x_{1}^{l_1} x_{2}^{l_2}\cdots x_{k}^{l_k}$ in $\Delta(x) P_{\mu}(x)$ 

, where $\Delta(x) = \prod\limits_{1 \leq i < j \leq k} (x_i-x_j)$ is the anti-symmetrizer and $P_\mu(x) = \prod\limits_j P_j(x_1,\cdots, x_k)^{i_j}$, where $P_j(x_1,\cdots, x_k) = x_1^j + \cdots + x_k^j$ is the $j$-th sum.

### [Murnaghan-Nakayama Rule](https://en.wikipedia.org/wiki/Murnaghan%E2%80%93Nakayama_rule)

Let $n$ and $t$ be positive integers, with $t \leq n$. Let $\sigma \in S_n$ be of the form $\sigma = \tau \cdot \rho$, where $\rho$ is a $t$-cycle, and $\tau$ is a permutation of $S_n$ with support disjoint from $\rho$. Let $\lambda$ be a partition of $n$.

Then
    $\chi^{\lambda}(\sigma) = \sum_{h \in \lambda, \, \ell(h) = t} (-1)^{\text{ht}(h)} \chi^{\lambda \backslash \text{bs}(h)}(\tau)$.


Above, $\chi_\lambda(\sigma)$ denotes the value of the character of the irreducible representation of $S_n$ corresponding to the partition $\lambda$, evaluated on the conjugacy class of $\sigma$, $\lambda \backslash \text{bs}(h)$ denotes the partition of $n - t$ obtained by removing the border strip $\text{bs}(h)$ from the Young diagram of $\lambda$, and $\chi_{\lambda \backslash \text{bs}(h)}(\tau)$ denotes the character value of the irreducible representation of $S_{n-t}$ corresponding to the partition $\lambda \backslash \text{bs}(h)$, evaluated on the conjugacy class of $\tau$.
<!-- 
## Results

## Conclusion -->
