# Character Tables of Symmetric Groups
## Contents
1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [Examples](#examples)
4. [Computations](#computations)
<!-- 5. Results
6. Conclusion -->

## Introduction
This is the repository for Lab of Geometry at Michigan (LoG(M): MATH 440) for the project Statistics of Character Tables of Symmetric Groups mentored by Prof. Sarah Peluse and Graduate Instructor Karthik Ganapathy. Our group members are Arnav Shah, Atharva Kulkarni, and Tony Zhang.

## Getting Started
Clone the repository on your laptop.

Gdrive link: [Click Here](https://drive.google.com/drive/folders/1J1zih494ypp2f18tC5mQWfqi8FodiMLS)  

Download all necessary files from the above link. Each function's description has the list of necessary files. Please read the [heatmaps and graphs](#heatmaps-and-graphs-char_tablepy) section for more information.

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

To compute the character table and print it, use the following code:
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

### Heatmaps and Graphs ([char_table.py](char_table.py))
Please read the instructions given at the begining of each code block in our Jupyter notebook file named [char_table.py](char_table.py).

Download all required files that are mentioned right before the code block from the gdrive link provided [above](#getting-started).

## Computations
### [Frobenius Formula](https://en.wikipedia.org/wiki/Frobenius_formula)
Given an integer partition $\lambda = \lambda_1 + \lambda_2 + \cdots + \lambda_k$ of $n$, let $\chi^{\lambda}$ be the corresponding irreducible character of $S_n$, and let $\chi^{\lambda}_{\mu}$ be short for the value of $\chi^{\lambda}$ at any $g$ with cycle type $\mu$, denote $l_j  =\lambda_j + k - j$, and $i_j$ the number of times $j$ appears in $\mu$, so $\sum\limits_j i_jj = n$, then we have the following Frobenius's Formula:

$\chi^{\lambda}_{\mu} = \text{coeff. of }  x_{1}^{l_1} x_{2}^{l_2}\cdots x_{k}^{l_k}$ in $\Delta(x) P_{\mu}(x)$ 

where $\Delta(x) = \prod\limits_{1 \leq i < j \leq k} (x_i-x_j)$ is the anti-symmetrizer and $P_\mu(x) = \prod\limits_j P_j(x_1,\cdots, x_k)^{i_j}$, where $P_j(x_1,\cdots, x_k) = x_1^j + \cdots + x_k^j$ is the $j$-th sum.

### [Murnaghan-Nakayama Rule](https://en.wikipedia.org/wiki/Murnaghan%E2%80%93Nakayama_rule)

Let $n$ and $t$ be positive integers, with $t \leq n$. Let $\sigma \in S_n$ be of the form $\sigma = \tau \cdot \rho$, where $\rho$ is a $t$-cycle, and $\tau$ is a permutation of $S_n$ with support disjoint from $\rho$. Let $\lambda$ be a partition of $n$.

Then
    $\chi^{\lambda}(\sigma) = \sum_{h \in \lambda, \, \ell(h) = t} (-1)^{\text{ht}(h)} \chi^{\lambda \backslash \text{bs}(h)}(\tau)$.


Above, $\chi_\lambda(\sigma)$ denotes the value of the character of the irreducible representation of $S_n$ corresponding to the partition $\lambda$, evaluated on the conjugacy class of $\sigma$, $\lambda \backslash \text{bs}(h)$ denotes the partition of $n - t$ obtained by removing the border strip $\text{bs}(h)$ from the Young diagram of $\lambda$, and $\chi_{\lambda \backslash \text{bs}(h)}(\tau)$ denotes the character value of the irreducible representation of $S_{n-t}$ corresponding to the partition $\lambda \backslash \text{bs}(h)$, evaluated on the conjugacy class of $\tau$.
<!-- 
## Results

## Conclusion -->
