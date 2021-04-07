import pandas as pd
import numpy as np
import os
periodic_table = {'H': 1.00797, 'HE': 4.0026, 'LI': 6.941, 'BE': 9.01218, 'B': 10.81, 'C': 12.011, 'N': 14.0067, 'O': 15.9994, 'F': 18.998403, 'NE': 20.179, 'NA': 22.98977, 'MG': 24.305, 'AL': 26.98154, 'SI': 28.0855, 'P': 30.97376, 'S': 32.06, 'CL': 35.453, 'AR': 39.948, 'K': 39.0983, 'CA': 40.08}
class process:
    def __init__(self, rmsd=False, energypoints=False):
        self.rmsd = rmsd
        self.energypoints = energypoints
        self.files = []
        
    def add(self, path):
        self.files.append(path)
        
    def extract_features(self):
        self.dataset = np.chararray((0,8 if self.rmsd and self.energypoints else 7 if self.rmsd or self.energypoints else 6), 32, True)
        for path in self.files: 
            if len(str(path)) > 5 and str(path)[-6:] == '.pdbqt': #Take just it .pdbqt
                pass
            else: #Take all .pdbqt into a folder
                for root, directories, files in os.walk(path, topdown=False):
                    if root == path: 
                        for name in files:
                            add = self.make_array(self.load(path, name), name)
                            self.dataset = np.append(self.dataset, add).reshape(self.dataset.shape[0] + add.shape[0], self.dataset.shape[1])
        if self.energypoints: self.dataset[:, 6] = self.dataset[:, 5].astype(np.float16) / self.dataset[:, 6].astype(np.float16)
        self.dataset[:, 2:] = self.dataset[:, 2:].astype(np.float16, copy=False)
        print(self.dataset[:2], self.dataset.shape)
            

    def load(self, path, name):
        pre_data = ''
        with open(os.path.join(path, name)) as fp:
            pre_data = fp.read().split('MODEL')
            del pre_data[0]
            fp.close()
        return pre_data


    def make_array(self, data, name):
        pos_data = np.chararray((len(data), self.dataset.shape[1]), 32, True)
        pos_data[:,0] = np.array([name.split('.pdbqt')[0]]*len(data), dtype=str)
        pos_data[:,1] = np.array([x for x in range(1, len(data)+1)], dtype=str)
        for i, ival in enumerate(data):
            result = self.take_position(ival)
            if len(result[1]) != 0:
                pos_data[i, 5] = result[0]
                centermass = self.cm(result[1])
                if self.energypoints: 
                    pos_data[i, 2:5] = centermass[:-1]
                    pos_data[i, 6] = centermass[-1]
                else: pos_data[i, 2:5] = centermass
        return pos_data
    
    
    def take_position(self, every_result):
        array = []
        energy = every_result.split('\n')[1].split()[3]
        for i in every_result.split('\n'):
            if 'HETATM' in i or 'ATOM' in i:
                inside = i.split()
                inside_numeric = []
                for x in inside:
                    try: inside_numeric.append(float(x))
                    except: pass
                array.append([self.find(''.join([x for x in inside[2] if not x.isdigit()]).upper())] + inside_numeric[2:5])
        return energy, np.array(array)
    
    def find(self, v):
        if v in periodic_table.keys(): return v
        else: return self.find(v[:-1])
        
        
    def cm(self, array):
        total_mass = np.sum([periodic_table[x] for x in array[:,0]])
        x = np.array([periodic_table[x[0]]*float(x[1]) for x in array[:,:2]], dtype=np.float16)
        y = np.array([periodic_table[x[0]]*float(x[1]) for x in array[:,[0,2]]], dtype=np.float16)
        z = np.array([periodic_table[x[0]]*float(x[1]) for x in array[:,[0,3]]], dtype=np.float16)
        if self.energypoints: 
            num_non_h = np.count_nonzero(array[:,0] != 'H')
            return np.array([sum(x)/total_mass, sum(y)/total_mass, sum(z)/total_mass, num_non_h], dtype=str)
        return np.array([sum(x)/total_mass, sum(y)/total_mass, sum(z)/total_mass], dtype=str)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    