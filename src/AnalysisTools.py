
from os import listdir
import pandas as pd
import numpy as np
import plotnine as pn

def plot_scan(data):
    """
    Assumes that you have cleaned and added aa class.
    """
    (
    pn.ggplot(data, pn.aes('position', 'I_sc', color='aa class'))
    + pn.geom_point()
    )

AA_CLASS = {
    'GLY':'nonpolar',
    'ALA':'nonpolar',
    'VAL':'nonpolar',
    'LEU':'nonpolar',
    'ILE':'nonpolar',
    'MET':'nonpolar',
    'PRO':'nonpolar',
    'PHE':'nonpolar',
    'TRP':'nonpolar',
    'SER':'polar',
    'THR':'polar',
    'TYR':'polar',
    'CYS':'polar',
    'ASN':'polar',
    'GLN':'polar',
    'ASP':'acidic',
    'GLU':'acidic',
    'LYS':'basic',
    'HIS':'basic',
    'ARG':'basic',
}

AA_CODE = {
    'GLY':'G',
    'ALA':'A',
    'VAL':'V',
    'LEU':'L',
    'ILE':'I',
    'MET':'M',
    'PRO':'P',
    'PHE':'F',
    'TRP':'W',
    'SER':'S',
    'THR':'T',
    'TYR':'Y',
    'CYS':'C',
    'ASN':'N',
    'GLN':'Q',
    'ASP':'D',
    'GLU':'E',
    'LYS':'K',
    'HIS':'H',
    'ARG':'R',
}

def assign_class(data):
    data['aa class'] = data['AA'].apply(lambda x: AA_CLASS[x])
    return data

LBIT_SEQUENCE = "VFTLEDFVGDWEQTAAYNLDQVLEQGGVSSLLQNLAVSVTPIQRIVRSGENALKIDIHVIIPYEGLSADQMAQIEEVFKVVYPVDDHHFKVILPYGTLVIDGVTPNMLNYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLITPDGSMLFRVTINS"

def res_LBit(n):
    return LBIT_SEQUENCE[n-1]

def add_labels(data):
    data['position'] = data['decoy'].apply(lambda x: int(x.split('_')[0][1:]))
    data['mutcode'] = data['decoy'].apply(lambda x: x.split('_')[0][1:]+x.split('_')[1][:3])
    data['AA'] = data['decoy'].apply(lambda x: x.split('_')[1][:3])
    assign_class(data)
    return data

def parseSmBits(date_time):
    """
    This probably does not work anymore.
    """
    kds = {'peptide86':0.7E-9,
           'peptide78':3.4E-9,
           'peptide79':8.5E-9,
           'peptide99':1.8E-7,
           'peptide128':2.8E-7,
           'native_test':0.9E-6,
           'peptide104':1.3E-6,
           'peptide101':2.5E-6,
           'peptide114':1.9E-4
           }

    l = listdir("./decoys/")
    fascs = []
    for f in l:
        if ".fasc" in f:
            if date_time in f:
                fascs.append(f)

    print(fascs)

    data = pd.read_json("./decoys/"+fascs[0], orient='records', lines=True)

    for f in fascs[1:]:
        d = pd.read_json("./decoys/"+f, orient='records', lines=True)
        data = data.append(d)

    print(data.shape)

    peptide_names = []
    for i, r in data.iterrows():
        s = r['filename'].split('-')
        s = s[0].split('/')
        name = s[2]
        peptide_names.append(name)

    data['peptide'] = peptide_names
    #print(data.head())

    peptides = set(data['peptide'])
    dfs = []

    for peptide in peptides:
        subdata = data[data['peptide']==peptide]
        subdata['kds'] = kds[peptide]

        stdev_total_score = np.std(subdata['total_score'])
        mean_total_score = np.mean(subdata['total_score'])
        print("{}: {} plus minus {}".format(peptide, mean_total_score, stdev_total_score))

        threshold_total_score = mean_total_score - 2 * stdev_total_score

        #dfs.append(subdata[subdata['total_score']<threshold_total_score])
        dfs.append(subdata.sort_values('total_score')[:10])

    significant_structures = pd.concat(dfs)
    print(significant_structures.shape)
    return significant_structures

def read_decoy_dir(directory, substr):
    jsons = listdir(directory)
    frames = []

    for f in jsons:
        if substr in f:
            df = pd.read_json(directory+f, orient='records', lines=True)
            df['filename'] = f
            frames.append(df)

    data = pd.concat(frames)
    
    return data

def get_top(df, cat_column, sort_column, n):
    cc = set(df[cat_column])

    topl = []

    for c in cc:
        sd = df[(df[cat_column]==c)].sort_values(sort_column)
        topl.append(sd[:n])

    top = pd.concat(topl)
    
    return top