import numpy as np
import sys

def save_result():
    aaDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
    List = list(aaDict.keys())
    List.sort()
    return List
if __name__ == "__main__":
    args = sys.argv[1:]
    file = args[0]
    read_dictionary = np.load(file).item()
    List = save_result()
    for amino1 in List:
        for amino2 in List:
            print (sum(sum(read_dictionary[amino1][amino2])))

