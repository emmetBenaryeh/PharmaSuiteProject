import sys
from rdkit import Chem
from rdkit.Chem import MolFromMolFile
from rdkit.Chem import MolFromMolBlock
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
#from rdkit.Chem.AllChem import CalcCrippenDescriptors
from rdkit.Chem.Draw import IPythonConsole

class SdfReader:
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname,'r',encoding='utf-8')
    
    def readSDF(self):
        """ 
        Return a list representing a molecule 
        """
        molList = []
        
        with self.doOpen() as fileH:
            name = ""
            molBlock = ""
            currentMol = None
            attribute = []

            for line in fileH:
                if line.startswith("F"):
                    if currentMol:  # Append previous molecule if exists
                        currentMol.append(attribute)
                        molList.append(currentMol)

                    name = line[1:].rstrip()
                    molBlock = ""
                    molBlock += line
                    while not line.startswith("M"):
                        line = next(fileH)
                        molBlock += line
                    currentMol = [name, molBlock]

                elif line.startswith(">"):
                    print("attribute"+line)
                    attribute.append(line.rstrip())
                    line = next(fileH)
                    print(line)
                    attribute.append(line.rstrip())
                    currentMol.append(attribute)
                    print("appended")
                    attribute = []
                elif line.startswith("$"):
                    if currentMol:  # Append previous molecule if exists
                        molList.append(currentMol)
                    currentMol = None
                    attribute = []

        return molList
    def writeSdf(self,input, outfile):
        """
        Writes an sdf file from a list of molecule data 
        """
        listIn = list(input)
        with open(outfile, 'w') as fileH:
            for item in listIn:
                fileH.write(item[1])
                fileH.write('\n')

                for attribute in item[2:]:
                    fileH.write(attribute[0])
                    fileH.write('\n')
                    fileH.write(attribute[1])
                    fileH.write('\n\n')
                fileH.write('$$$$\n')
    def returnVals(self,list):
        clogP = None
        hAcceptors = None
        hDonors = None
        molMass = None
        numAtoms = None
        rotBonds = None
        tpsa = None
        logP = None
        molRefract = None
        dockScore = None 
        molecule = MolFromMolBlock(list[1])
        for attribute in list[2:]:
            if "clogP" in attribute[0]:
                clogP = attribute[1]
            if "cceptor" in attribute[0]:
                hAcceptors = attribute[1]
            if "onor" in attribute[0]:
                hDonors = attribute[1]
            if "MW" in attribute[0]:
                molMass = attribute[1]
            if "TPSA" in attribute[0]:
                tpsa = attribute[1] 
            if "RotBonds" in attribute[0]:
                rotBonds = attribute[1]
            if "Human Fructose Aldolase Docking" in attribute[0]:
                dockScore = attribute[1]
        numAtoms = molecule.GetNumAtoms()
        logP,molRefract = AllChem.CalcCrippenDescriptors(molecule)
        return [clogP, hAcceptors, hDonors, molMass, numAtoms, rotBonds, tpsa, logP, molRefract, dockScore]