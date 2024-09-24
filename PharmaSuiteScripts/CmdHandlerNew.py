import subprocess
import shutil
from SDFREADERNEW import SdfReader
from openbabel import pybel
from openbabel import openbabel as OB
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromMolBlock
from rdkit.Chem import MolFromMolFile
from rdkit.Chem import AddHs
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit.Chem import rdDetermineBonds
from dataImager import sortAndPlot
import os
import re
import sys

class Cmdhandler():
    def __init__(self,library,enzyme,name):
        print("init used")
        self.library = library #sdf file with ligands
        self.enzyme = enzyme #pdbqt of trimmed enzyme
        self.name = name #options? idk what tf this will be 
    def getMols(self):
        print("getMols used")
        Sdf = SdfReader(self.library)
        for i,value in enumerate(Sdf.readSDF()):
            filename = os.path.join("C:/adt/Lib/site-packages/AutoDockTools/DockingAnalysis/LigandsRaw", f"LigandEmmet_{i}.sdf")
            mol = Chem.MolFromMolBlock(value[1])
            template = Chem.MolToSmiles(mol,kekuleSmiles=True)
            template2 = Chem.MolFromSmiles(template)
            mol = AllChem.AssignBondOrdersFromTemplate(template2,mol)
            mol = Chem.rdmolops.AddHs(mol,addCoords=True)
            AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
            with open(filename,'w') as filename:
                filename.write(Chem.MolToMolBlock(mol))
    def convert_sdf_to_pdb(self,sdf_file, pdb_file):
        print("convert Sdf to pdb used")
        try:
            # Load molecules from SDF file
            mol = next(pybel.readfile("sdf", sdf_file))

            # Write molecule to PDBQT file
            mol.write("pdb", pdb_file, overwrite=True)

            print("Conversion successful sdf-pdb")
        except Exception as e:
            print("Error:", e)
    def getBindingAffinity(self,file):
        with open(file,'r') as infile:
            line = infile.readline()
            line = infile.readline()
        match = re.search(r"\s+([-]?\d+\.?\d*)", line)
        if match:
            first_number = float(match.group())
            return first_number
        else:
            return "error"
    def run_vina_for_ligands(self,working_directory, vina_executable,):
        # Change to the working directory
        os.chdir(working_directory)

        # Loop through all .pdbqt files in the working directory
        for ligand_file in os.listdir():
            if ligand_file.endswith(".pdbqt"):
                if ligand_file.startswith("LigandEmmet_"):
                    print(f"Running vina.exe for {ligand_file}")
                    # Construct the full path to the ligand file
                    ligand_path = os.path.join(working_directory, ligand_file)
                    # Define the output file path
                    output_file = os.path.splitext(ligand_path)[0] + "Out.pdbqt"
                    # Run vina.exe with subprocess
                    subprocess.run([vina_executable, "--ligand", ligand_path, "--out", output_file, "--config", self.enzyme])
    def CleanUpFile(self,directory):
        #Move Files to Clean directories
        for file in os.listdir(directory):
            if file.startswith("LigandEmmet"):
                if file.endswith("Out.pdbqt"):
                    shutil.move(os.path.join("C:/adt/Lib/site-packages/AutoDockTools", file),"C:/Program Files/PharmaSuiteScripts/Output")
                if file.endswith("M.pdb"):
                    shutil.move(os.path.join("C:/adt/Lib/site-packages/AutoDockTools", file),"C:/adt/Lib/site-packages/AutoDockTools/DockingAnalysis/LigandsPDB")
                if file.endswith("M.pdbqt"):
                    shutil.move(os.path.join("C:/adt/Lib/site-packages/AutoDockTools", file),"C:/adt/Lib/site-packages/AutoDockTools/DockingAnalysis/LigandsPDBQT")
def main():
    """
    if len(sys.argv) != 4:
        print("Usage: python main.py infile enzymefile name")
        sys.exit(1)
    infile = sys.argv[1]
    enzymefile = sys.argv[2]
    name = sys.argv[3]
    """

    infile = 'FructoseShortTest.sdf'
    enzymefile = "C:/adt/Lib/site-packages/AutoDockTools/configNew.txt"
    name = "Fructose Test"
    
    data = Cmdhandler(infile,enzymefile,name)
    reader = SdfReader(infile)
    CList = reader.readSDF()
    data.getMols()
    
    print("start convert to pdb")
    for i, filename in enumerate(os.listdir("C:/adt/Lib/site-packages/AutoDockTools/DockingAnalysis/LigandsRaw")):
        print("working on" + filename)
        # Construct the input and output file paths
        input_file = os.path.join("C:/adt/Lib/site-packages/AutoDockTools/DockingAnalysis/LigandsRaw", filename)
        output_file = os.path.join("C:/adt/Lib/site-packages/AutoDockTools", f'{os.path.splitext(filename)[0]}M.pdb')            
        data.convert_sdf_to_pdb(input_file, output_file)
        print(input_file +" converted to " + output_file)
    print("sdf converted to pdb /n start convert to pdbqt")
    
    batchPrepLigand = "C:/adt/Lib/site-packages/AutoDockTools/PrepareLigand2.bat"
    os.chdir("C:/adt/Lib/site-packages/AutoDockTools/")
    subprocess.run([batchPrepLigand], shell=True)

    data.run_vina_for_ligands("C:/adt/Lib/site-packages/AutoDockTools","C:/adt/Lib/site-packages/AutoDockTools/vina.exe")
    data.CleanUpFile("C:/adt/Lib/site-packages/AutoDockTools")
    #move output data to same file as script
    for index, item in enumerate(CList):
        for file in os.listdir("C:/Program Files/PharmaSuiteScripts/Output"):
            filePath = os.path.join("C:/Program Files/PharmaSuiteScripts/Output",file)
            name2 = "LigandEmmet_"+str(index)+"MOut.pdbqt"
            if str(name2) in file:
                score = str(data.getBindingAffinity(filePath))
                item.append(["> <Human Fructose Aldolase Docking>",score])
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    tempFile = data.name+'_wDock.sdf'
    reader.writeSdf(CList,tempFile)
    plotter = sortAndPlot(tempFile,data.name)
    plotter.generateDockingScorePlot
    plotter.generateRO5Histogram()
    plotter.generateRO5Histogram()
if __name__ == "__main__" :
    main()