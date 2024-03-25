from typing import List
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.utilities.iterables import partitions

R, *x = ring("x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10", ZZ)


def convert_dict_partitions_to_list(n: int):
    """
    uses the builtin partition module and convert that to a python list
    """
    index = list(partitions(n))
    partitions_list = []
    for pair in index:
        tmp = []
        for key in pair:
            for _ in range(pair[key]):
                tmp.append(key)
        partitions_list.append(tmp)
    return partitions_list


def Delta(k: int):
    """
    The anti-symmetrizer in the Frobenius formula
    """
    return R.mul([x[i] - x[j] for j in range(2, k + 1) for i in range(1, j)])


def P(j: int, k: int):
    """
    The product of power sum in the Frobenius formula
    """
    # return the k-th power sum of the list of variables x
    return sum([x[i] ** j for i in range(1, k + 1)])


def frobenius(irrep: List[int], conj_class: List[int]):
    """
    calculate the product of the polynomial and get the character value using the monomial that corresponds to the Specht Module
    """
    kv = {}
    for num in conj_class:
        if num not in kv:
            kv[num] = conj_class.count(num)
    P_conj_class = R.mul([R.mul([P(i, len(irrep))
                         for _ in range(kv[i])]) for i in kv])
    P_conj_class = R.mul([P_conj_class, Delta(len(irrep))])
    monomial = R.mul(
        [x[i] ** (len(irrep) + irrep[i - 1] - i)
         for i in range(1, len(irrep) + 1)]
    )
    return P_conj_class.coeff(monomial)


def get_character_table(n: int):
    """
    get the final character table by computing each slot using the Frobenius formula.
    """
    partition_list = convert_dict_partitions_to_list(n)
    character_table = [
        [0 for _ in range(len(partition_list))] for _ in range(len(partition_list))
    ]

    for i, partition_i in enumerate(partition_list):
        for j, partition_j in enumerate(partition_list):
            character_table[i][len(partition_list) - j - 1] = frobenius(
                partition_i, partition_j
            )
    return character_table


if __name__ == "__main__":
    N = 6
    char_table = get_character_table(N)
    print(char_table)
