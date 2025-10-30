#!/usr/bin/env python
import os
import numpy as np
import codecs
from pymatgen.core.periodic_table import Element

class electride:
    """
    Parse electride from a directory containing ELFCAR, PARCHG and vasprun.xml
    """

    def __init__(self, path = './',
        vaspxml = 'vasprun.xml',
        cmd = 'bader ',
        ELF_min = 0.60):

        self.label = ['e0025', 'e05', 'e10', 'band0', 'band1']
        self.volume = [0]*len(self.label)
        self.error = None
        self.ELF_maxima = []

        # Save current directory and change to path for all file operations
        original_dir = os.getcwd()
        path = os.path.abspath(path)
        
        if not os.path.exists(os.path.join(path, 'ELFCAR')):
            self.error = 'ELFCAR not exist'
            return
        
        # Change to structure directory for all file operations
        os.chdir(path)
        
        try:
            ELFCAR_path = 'ELFCAR'
            PARCHG_path = 'PARCHG'
            clean_cmd1 = 'rm -f ELFCAR-m bader_log ACF.dat BCF.dat AVF.dat'
            clean_cmd2 = 'rm -f PARCHG-m bader_log ACF.dat BCF.dat AVF.dat'

            cell, coor, radii, grid = self.Read_ELFCAR(ELFCAR_path, ELF_min)
            os.system(cmd + 'ELFCAR-m > bader_log')
            if not os.path.exists('BCF.dat'):
                self.error = 'bader error in parsing ELFCAR'
            else:
                self.ELF_maxima, self.fac1 = self.Read_BCF('BCF.dat', radii, cell)
            os.system(clean_cmd1)

            if len(self.ELF_maxima) > 0:
                self.cell = cell
                for i, label in enumerate(self.label):
                    fpath = PARCHG_path + '-' + label
                    if os.path.exists(fpath):
                        self.ModifyPARCHG(self.ELF_maxima, fpath)
                        os.system(cmd + 'PARCHG-m -vac off> bader_log')

                        if os.path.exists('ACF.dat'):
                            charge = self.Read_ACF('ACF.dat')
                            os.system(clean_cmd2)
                            self.volume[i] = self.parse_charge(charge)
                        else:
                            self.error = f'Bader analysis failed for {fpath}'
        finally:
            # Always return to original directory
            os.chdir(original_dir)

    def parse_charge(self, charge):
        totalv = 0.0
        for i in range(len(self.ELF_maxima)):
            if charge[i][0] > 0.01:
                totalv += charge[i][2]
        volume = totalv/np.linalg.det(self.cell)*100
        return volume

    @staticmethod
    def Read_ACF(filename):
        f = open(filename, 'rb')
        input_content = f.readlines()
        f.close()
        chg = []
        for i in range(len(input_content)):
            s = input_content[i].split()
            if s[0].isdigit():
               a=[float(f) for f in s]
               chg.append(a[4:])
        return np.array(chg)

    @staticmethod
    def Read_BCF(filename, radii, cell):
        f = open(filename, 'rb')
        input_content = f.readlines()
        f.close()
        pos = []
        count = 0
        cell = np.linalg.inv(cell)
        for i in range(len(input_content)):
            s = input_content[i].split()
            if s[0].isdigit():
               a=[float(f) for f in s]
               if 1.2*radii[int(a[-2])-1] < a[-1]:
                  count += 1
                  if len(pos) < 4:
                     pos.append(np.dot(np.array(a[1:4]), cell))
        if len(pos) > 0:
           fac = count/len(pos)
        else:
           fac = 1.0
        return np.array(pos), fac

    @staticmethod
    def ModifyPARCHG(pos, filename):
        """
        This module modifies PARCHG for bader analysis
        """
        f1  = codecs.open(filename, 'rb', encoding='utf-8')
        input_content = f1.readlines()
        f1.close()
        f2 = open('PARCHG-m', 'w')
        f2.writelines(input_content[:5])
        f2.write(' H '+input_content[5])
        f2.write(str(len(pos)) + ' ' + input_content[6])
        f2.write(input_content[7])
        for coor in pos:
            f2.write('%10.6f %9.6f %9.6f\n' % (coor[0], coor[1], coor[2]))
        f2.writelines(input_content[8:])
        f2.close()

    @staticmethod
    def Read_ELFCAR(filename, ELF_min):
        """
        This module reads ELFCAR
        """
        f = open(filename, 'rb')
        f1 = open('ELFCAR-m', 'w')
        input_content = f.readlines()
        f.close()
        count = 0
        cell = []
        coor = []
        ELF_raw = []
        N_atoms = 0
        grid = []
        for line in input_content:
            line=str(line,'utf-8')
            count = count + 1
            if count < 3:
               f1.write(line)
            elif (count>=3) and (count<=5):
               cell.append([float(f) for f in line.split()])
               f1.write(line)
            elif count==6:
               f1.write(line)
               symbol = line.split()
            elif count==7:
               f1.write(line)
               numIons = [int(f) for f in line.split()]
               N_atoms = sum(numIons)
               f1.write('Direct\n')
            elif (count>=9) and (count<9+N_atoms):
               f1.write(line)
               coor.append([float(f) for f in line.split()])
            elif count == 10+N_atoms:
               f1.write('\n')
               f1.write(line)
               grid = [int(f) for f in line.split()]
            elif count > 10+N_atoms:
               ELF_raw = line.split()
               for i, f in enumerate(ELF_raw):
                   if float(f)<ELF_min:
                      f = '0.00000'
                   if i==0:
                      f1.write('%8s' % (f))
                   else:
                      f1.write('%12s' % (f))
               f1.write('\n')
        f1.close()

        radii = np.array([])
        for ele in range(len(symbol)):
            # Use atomic_radius from pymatgen (in Angstroms)
            radius = Element(symbol[ele]).atomic_radius
            if radius is None:
                radius = 1.0  # Default fallback for elements without atomic_radius
            else:
                radius = float(radius)  # Convert FloatWithUnit to plain float
            tmp = radius * np.ones(numIons[ele])
            radii = np.append(radii, tmp)
        cell = np.array(cell)
        return cell, coor, radii, grid

def get_subdir(a_dir):
   return sorted([name for name in os.listdir(a_dir)
           if os.path.isdir(os.path.join(a_dir, name))
           and name not in ['mp_mattersim_cache', '.git', '__pycache__']])

if __name__ == "__main__":
    from optparse import OptionParser
    import pandas as pd
    from tabulate import tabulate

    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="dir", default='./',
                      help="directory containing PARCHG and ELFCAR", metavar="dir")
    parser.add_option("-p",  "--pdir", dest='pdir',
                  help="by parent directory")
    parser.add_option("-b", "--bader-exe", dest='bader_exe', default='bader',
                  help="path to bader executable (default: bader)")
    parser.add_option("-t", "--threshold", dest='threshold', type='float', default=0.6,
                  help="ELF threshold value (default: 0.6)")

    (options, args) = parser.parse_args()
    total_dir = []
    if options.pdir is None:
       total_dir.append(options.dir)
    else:
       total_dir =  get_subdir(options.pdir)
       os.chdir(options.pdir)

    col_name = {'formula':[]}
    for label in ['e0025','e05', 'e10', 'band0', 'band1']:
        col_name[label] = []

    for subdir in total_dir:
        print(subdir)
        test = electride(subdir, cmd=options.bader_exe + ' ', ELF_min=options.threshold)
        
        # Extract structure name from path (e.g., "Li10B1N4_s001" from full path)
        # Handles both "path/to/Li10B1N4_s001/ELF" and "Li10B1N4_s001"
        import re
        struct_name = subdir
        match = re.search(r'([A-Z][a-z]?\d+(?:[A-Z][a-z]?\d+)*_s\d+)', subdir)
        if match:
            struct_name = match.group(1)
        else:
            # Fallback: just use basename
            struct_name = os.path.basename(subdir.rstrip('/'))
        
        col_name['formula'].append(struct_name)
        
        # Print error if any
        if test.error:
            print(f"  Warning: {test.error}")

        for i, label in enumerate(['e0025','e05', 'e10', 'band0', 'band1']):
            col_name[label].append(test.volume[i])

    df = pd.DataFrame(col_name)
    print(tabulate(df, headers='keys', tablefmt='psql'))
    df.to_csv('electride_analysis.csv')
