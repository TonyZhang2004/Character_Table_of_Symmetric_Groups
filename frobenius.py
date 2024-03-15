from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.utilities.iterables import partitions

# GLOBAL VARIABLES
N = 3  # Size of Character Table

R, *x = ring("x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10", ZZ)

def convert_dict_partitions_to_list(N):
    """
    TODO: What does this funciton do?
    """
    index = list(partitions(N))
    partitions_list = []
    for pair in index:
        tmp = []
        for key in pair:
            for _ in range(pair[key]):
                tmp.append(key)
        partitions_list.append(tmp)
    return partitions_list


def Delta(k):
    return R.mul([x[i] - x[j] for j in range(2, k + 1) for i in range(1, j)])


def P(j, k):
    # return the k-th power sum of the list of variables x
    return sum([x[i] ** j for i in range(1, k + 1)])


def character(irrep: list, conj_class: list):
    kv = {}
    for num in conj_class:
        if num not in kv:
            kv[num] = conj_class.count(num)
    P_conj_class = R.mul([R.mul([P(i, len(irrep)) for _ in range(kv[i])]) for i in kv])
    P_conj_class = R.mul([P_conj_class, Delta(len(irrep))])
    monomial = R.mul(
        [x[i] ** (len(irrep) + irrep[i - 1] - i) for i in range(1, len(irrep) + 1)]
    )
    return P_conj_class.coeff(monomial)


if __name__ == "__main__":
    partition_list = convert_dict_partitions_to_list()
    character_table = [
        [0 for _ in range(len(partition_list))] for _ in range(len(partition_list))
    ]

    for i in range(len(partition_list)):
        for j in range(len(partition_list)):
            character_table[i][len(partition_list) - j - 1] = character(
                partition_list[i], partition_list[j]
            )

    print(character_table)
