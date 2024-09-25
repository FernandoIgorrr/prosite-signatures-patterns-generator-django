from enum import Enum
from typing import List, Optional

class Classification(Enum):
    HYDROFOBIC = 'hydrofobic'
    NONPOLAR_ALIPHATIC = 'Nonpolar, aliphatic'
    AROMATIC = 'aromatic'
    POLAR_UNCHARGED = 'polar, uncharged'
    POSITIVELY_CHARGED = 'positively charged'
    NEGATIVELY_CHARGED = 'negatively charged'
    CYSTEINES = 'cysteines'
    GLYCINES = 'glycines'
    PROLINES = 'prolines'
    NONE = 'none'


class L1(Enum):
    A = 'A'
    C = 'C'
    D = 'D'
    E = 'E'
    F = 'F'
    G = 'G'
    H = 'H'
    I = 'I'
    K = 'K'
    L = 'L'
    M = 'M'
    N = 'N'
    P = 'P'
    Q = 'Q'
    R = 'R'
    S = 'S'
    T = 'T'
    V = 'V'
    W = 'W'
    Y = 'Y'
    NONE = '0'
    GAP_DOT = '.'
    GAP_HIFEN = '-'


class L3(Enum):
    ALA = 'Ala'
    CYS = 'Cys'
    ASP = 'Asp'
    GLU = 'Glu'
    PHE = 'Phe'
    GLY = 'Gly'
    HIS = 'His'
    ILE = 'Ile'
    LYS = 'Lys'
    LEU = 'Leu'
    MET = 'Met'
    ASN = 'Asn'
    PRO = 'Pro'
    GLN = 'Gln'
    ARG = 'Arg'
    SER = 'Ser'
    THR = 'Thr'
    VAL = 'Val'
    TRP = 'Trp'
    TYR = 'Tyr'
    NONE = '0'
    GAP_DOT = '.'
    GAP_HIFEN = '-'


class Aminoacid:
    AMINOACIDS: List['Aminoacid'] = []

    def __init__(self, name: str, l1: L1, l3: L3, classification: Classification):
        self.name = name
        self.l1 = l1
        self.l3 = l3
        self.classification = classification

    @classmethod
    def get_aminoacid(cls, letter: Optional[str]) -> 'Aminoacid':
        if letter is None:
            return cls('Gap', L1.NONE, L3.NONE, Classification.NONE)

        aminoacid_found = next((a for a in cls.AMINOACIDS if a.l1 == letter), None)
        return aminoacid_found or cls('Gap', L1.GAP_DOT, L3.GAP_DOT, Classification.NONE)

    def equals(self, aminoacid: 'Aminoacid') -> bool:
        return (self.l1 == aminoacid.l1 and
                self.l3 == aminoacid.l3 and
                self.name == aminoacid.name and
                self.classification == aminoacid.classification)

    def equals_classification(self, aminoacid: 'Aminoacid') -> bool:
        return self.classification == aminoacid.classification

    def substitution_score(self, aminoacid: 'Aminoacid', matrix: dict) -> int:
        return matrix[self.l1.value][aminoacid.l1.value]

    def substitution_validation(self, aminoacid: 'Aminoacid', matrix: dict, threshold: int) -> bool:
        return self.substitution_score(aminoacid, matrix) >= threshold

    def show_aminoacid(self) -> str:
        return f'Nome: {self.name} | Símbolo: {self.l1.value} | Iniciais: {self.l3.value} | Classificação: {self.classification.value}'


# Inicialização da lista de aminoácidos
Aminoacid.AMINOACIDS = [
    Aminoacid('Alanine', L1.A, L3.ALA, Classification.HYDROFOBIC),
    Aminoacid('Isoleucine', L1.I, L3.ILE, Classification.HYDROFOBIC),
    Aminoacid('Leucine', L1.L, L3.LEU, Classification.HYDROFOBIC),
    Aminoacid('Methionine', L1.M, L3.MET, Classification.HYDROFOBIC),
    Aminoacid('Phenylalanine', L1.F, L3.PHE, Classification.HYDROFOBIC),
    Aminoacid('Tryptophan', L1.W, L3.TRP, Classification.HYDROFOBIC),
    Aminoacid('Valine', L1.V, L3.VAL, Classification.HYDROFOBIC),
    Aminoacid('Lysine', L1.K, L3.LYS, Classification.POSITIVELY_CHARGED),
    Aminoacid('Arginine', L1.R, L3.ARG, Classification.POSITIVELY_CHARGED),
    Aminoacid('Glutamate', L1.E, L3.GLU, Classification.NEGATIVELY_CHARGED),
    Aminoacid('Aspartate', L1.D, L3.ASP, Classification.NEGATIVELY_CHARGED),
    Aminoacid('Asparagine', L1.N, L3.ASN, Classification.POLAR_UNCHARGED),
    Aminoacid('Glutamine', L1.Q, L3.GLN, Classification.POLAR_UNCHARGED),
    Aminoacid('Serine', L1.S, L3.SER, Classification.POLAR_UNCHARGED),
    Aminoacid('Threonine', L1.T, L3.THR, Classification.POLAR_UNCHARGED),
    Aminoacid('Cysteine', L1.C, L3.CYS, Classification.CYSTEINES),
    Aminoacid('Glycine', L1.G, L3.GLY, Classification.GLYCINES),
    Aminoacid('Proline', L1.P, L3.PRO, Classification.PROLINES),
    Aminoacid('Histidine', L1.H, L3.HIS, Classification.AROMATIC),
    Aminoacid('Tyrosine', L1.Y, L3.TYR, Classification.AROMATIC),
    Aminoacid('None', L1.NONE, L3.NONE, Classification.NONE),
    Aminoacid('Gap', L1.GAP_DOT, L3.GAP_DOT, Classification.NONE),
    Aminoacid('Gap', L1.GAP_HIFEN, L3.GAP_HIFEN, Classification.NONE),
]