'''
Date: 2021-08-18 11:05:39
LastEditors: jiyuyang
LastEditTime: 2021-11-08 15:16:35
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import json
import re
import string
from os import PathLike
from typing import List, Sequence, Union

import numpy as np


def remove_empty(a: list) -> list:
    """Remove '' and [] in `a`"""
    while '' in a:
        a.remove('')
    while [] in a:
        a.remove([])


def skip_notes(line: str) -> str:
    """Delete comments lines with '#' or '//'

    :params line: line will be handled 
    """
    line = re.compile(r"#.*").sub("", line)
    line = re.compile(r"//.*").sub("", line)
    line = line.strip()
    return line


def list_elem2strip(a: List[str], ds=string.whitespace) -> List[str]:
    """Strip element of list with `str` type"""
    def list_strip(s):
        return s.strip(ds)
    return list(map(list_strip, a))


def search_sentence(file: PathLike, sentence: str) -> Union[bool, str]:
    """Search sentence in file

    :params file: file descriptor
    :params sentence: sentence will be searched
    """

    if isinstance(sentence, str):
        sentence = sentence.strip()
        for line in file:
            line = skip_notes(line).strip()
            if line == sentence:
                return line
    elif isinstance(sentence, list):
        sentence = list_elem2strip(sentence)
        for line in file:
            line = skip_notes(line).strip()
            if line in sentence:
                return line

    file.seek(0, 0)
    return False


def list_elem_2int(a: List[str]) -> List[int]:
    """Convert type of list element to int

    :params a: 1-D list
    """
    return list(map(int, a))


def list_elem_2float(a: List[str]) -> List[float]:
    """Convert type of list element to float

    :params a: 1-D list
    """
    return list(map(float, a))

# KPT


class Kpt:
    """K-points information"""

    def __init__(self, mode: str, numbers: list = [], special_k: list = [], offset: list = [0.0, 0.0, 0.0], klabel: list = []) -> None:
        """Initialize Kpt object

        :params mode: ‘Gamma’, ‘MP’ or 'Line'
        :params numbers: for ‘Gamma’ and ‘MP’ mode, list of three integers, for 'Line' mode, list of number of k points between two adjacent special k points.
        :params special_k: list of special k-points
        :params offset: offset of the k grid. Default: [0.0, 0.0, 0.0]
        """

        self.mode = mode
        self.numbers = numbers
        self.special_k = special_k
        self.offset = offset
        self.klabel = klabel

    def get_kpt(self) -> str:
        """Return the `KPT` file as a string"""

        line = []
        line.append("K_POINTS")
        if self.mode in ["Gamma", "MP"]:
            line.append("0")
            line.append(self.mode)
            line.append(" ".join(list_elem2str(self.numbers+self.offset)))
        elif self.mode == "Line":
            line.append(str(len(self.special_k)))
            line.append(self.mode)
            for i, k in enumerate(self.special_k):
                if self.klabel:
                    line.append(" ".join(list_elem2str(
                        k+[self.numbers[i]]))+'\t#'+self.klabel[i])
                else:
                    line.append(" ".join(list_elem2str(k+[self.numbers[i]])))

        return '\n'.join(line)

    def write_kpt(self, filename: str):
        """Write k-points file

        :params filename: absolute path of k-points file
        """

        with open(filename, 'w') as file:
            file.write(self.get_kpt())

    @property
    def full_kpath(self):
        """Uniform k-point path"""

        total_k = np.sum(self.numbers)
        spec_k_coor = np.array(self.special_k)
        interval = (np.roll(spec_k_coor, -1, axis=0) - spec_k_coor) / \
            np.reshape(self.numbers, (-1, 1))
        max_num = np.max(self.numbers)
        len_num = len(self.numbers)
        k_coor_span = np.zeros((len_num, max_num), dtype=np.float)
        X, Y, Z = np.split(spec_k_coor, 3, axis=1)
        i_X, i_Y, i_Z = np.split(interval, 3, axis=1)
        for i, j in enumerate(self.numbers):
            k_coor_span[i][:j] = np.arange(j)
        X = (i_X * k_coor_span + X.repeat(max_num, axis=1)).flatten()
        Y = (i_Y * k_coor_span + Y.repeat(max_num, axis=1)).flatten()
        Z = (i_Z * k_coor_span + Z.repeat(max_num, axis=1)).flatten()
        k_direct_coor = np.empty((3, total_k), dtype=np.float)
        k_direct_coor[0] = X[:total_k]
        k_direct_coor[1] = Y[:total_k]
        k_direct_coor[2] = Z[:total_k]

        return k_direct_coor.T

    @property
    def label_special_k(self):
        """Label special k-points based on `numbers` list"""

        index = np.cumsum(np.concatenate(([1], self.numbers), axis=0))[
            :len(self.special_k)]
        if self.klabel:
            return zip(self.klabel, index)
        else:
            return index


def read_kpt(kpt_file: PathLike) -> Kpt:
    """Read k-points file

    :params kpt_file: absolute path of k-points file
    """
    number = 0
    with open(kpt_file, "r") as file:
        if search_sentence(file, ["K_POINTS", "KPOINTS", "K"]):
            number = int(file.readline())
        mode = search_sentence(file, ["Gamma", "MP", "Line"])
        if mode in ["Gamma", "MP"]:
            line = skip_notes(file.readline()).split()
            numbers = list_elem_2int(line[:3])
            offset = list_elem_2float(line[3:])
            return Kpt(mode, numbers, special_k=[], offset=offset)
        elif mode == "Line":
            special_k = []
            numbers = []
            klabel = []
            for k in range(number):
                line = file.readline()
                if re.match("#", line):
                    continue
                else:
                    linesplit = line.split(maxsplit=4)
                special_k.append(list_elem_2float(linesplit[:3]))
                numbers.append(int(linesplit[3]))
                if len(linesplit) == 5:
                    klabel.append(linesplit[4].strip('#\n '))
            return Kpt(mode, numbers, special_k, offset=[], klabel=klabel)


def read_json(filename: PathLike) -> dict:
    """ Read json file and return dict

    :params filename: json file
    """
    with open(filename, 'r') as file:
        text = json.load(file)
    return text


def list_elem2str(a: Union[List[float], List[int]]) -> List[str]:
    """Convert type of list element to str

    :params a: 1-D list
    """
    return list(map(str, a))


def energy_minus_efermi(energy: Sequence, efermi: float) -> np.ndarray:
    """Return energy after subtracting the Fermi level

    :params efermi: Fermi level in unit eV
    """

    return np.array(energy)-efermi


angular_momentum_label = ['s', 'p', 'd', 'f', 'g']


def get_angular_momentum_label(l_index: int) -> str:
    """Atomic orbital angular momentum label from l_index

    :params l_index: 0 or 1 or 2 or 3 or 4
    """

    return angular_momentum_label[l_index]


angular_momentum_name = [
    ['$s$'],
    ['$p_x$', '$p_y$', '$p_z$'],
    ['$d_{3z^2-r^2}$', '$d_{xy}$', '$d_{xz}$', '$d_{x^2-y^2}$', '$d_{yz}$'],
    ['$f_{5z^2-3r^2}$', '$f_{5xz^2-xr^2}$', '$f_{5yz^2-yr^2}$',
        '$f_{zx^2-zy^2}$', '$f_{xyz}$', '$f_{x^3-3*xy^2}$', '$f_{3yx^2-y^3}$'],
    ['$g_1$', '$g_2$', '$g_3$', '$g_4$', '$g_5$',
        '$g_6$', '$g_7$', '$g_8$', '$g_9$']
]


def get_angular_momentum_name(l_index: int, m_index: int) -> str:
    """Atomic orbital angular momentum name from l_index and m_index

    :params l_index: 0 or 1 or 2 or 3 or 4
    :params m_index: 0 ... 2*l_index
    """

    return angular_momentum_name[l_index][m_index]
