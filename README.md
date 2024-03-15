# Character Tables of Symmetric Groups
## Contents
1. Introduction
2. Getting Started
3. Example
4. Computations
5. Results
6. Conclusion

## Introduction
This is the repository for Lab of Geometry at Michigan (MATH 440: LoG(M)) for the project Statistics of Character Tables of Symmetric Groups mentored by Prof. Sarah Peluse and Graduate Instructor Karthik Ganapathy.

## Getting Started
Gdrive link: [Click Here](https://drive.google.com/drive/folders/1J1zih494ypp2f18tC5mQWfqi8FodiMLS)  
Download necessary files from the above link. Each function's description has the list of necessary files.

## Example
### [Frobenius Formula](frobenius.py)
Change the variable `N` to change the size of the character table.

    N = 3  # Size of Character Table

Run:
    python frobenius.py

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

## Results

## Conclusion
