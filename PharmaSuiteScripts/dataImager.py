from SDFREADERNEW import SdfReader
from sortData import sortData


import numpy as np
import matplotlib.pyplot as plt 
class sortAndPlot:
    def __init__(self,infile,name):
        #initialize variables, multiple zones for dock score bc histogram
        self.runSortData = sortData(infile)
        self.validListofList = self.runSortData.checkvalidity()
        
        #number of drugs
        self.numOfCandidates = len(self.validListofList)
        self.name = name
        self.clogTot = 0
        self.hAcceptorsTot = 0
        self.hDonorsTot = 0
        self.molMassTot = 0
        self.numAtomsTot = 0
        self.rotBondsTot = 0
        self.tpsaTot = 0
        self.logPTot = 0
        self.molRefractTot = 0
        self.dockScoreZ1Tot = 0
        self.dockScoreZ2Tot = 0
        self.dockScoreZ3Tot = 0
        self.ghosePassTot = 0
        self.rule3Tot = 0
        
        self.get_the_tots()
        
    def get_the_tots(self):
        #checks each field in every daughter list to make an approriate sum for each class variable
        for i in self.validListofList:
            if i[0] == 1:
                self.clogTot += 1
            if i[1] == 1:
                self.hAcceptorsTot += 1
            if i[2] == 1:
                self.hDonorsTot += 1
            if i[3] == 1:
                self.molMassTot += 1
            if i[4] == 1:
                self.numAtomsTot += 1
            if i[5] == 1:
                self.rotBondsTot += 1
            if i[6] == 1:
                self.tpsaTot += 1
            if i[7] == 1:
                self.logPTot += 1
            if i[8] == 1:
                self.molRefractTot += 1
                
            #specifies which bar in dock score bar graph should be incremented
            if (i[9][1] is not None):
                if (i[9][0] == 1) or (i[9][0] == 0):
                    #change to all possible docking scores since we want then to be in the plot still 
                    if (float(i[9][1]) > -3.6) and (float(i[9][1]) <= 0):
                        self.dockScoreZ1Tot += 1
                    if (float(i[9][1]) > -7.2) and (float(i[9][1]) <= -3.6):
                        self.dockScoreZ2Tot += 1
                    if (float(i[9][1]) >= -11) and (float(i[9][1]) <= -7.2):
                        self.dockScoreZ3Tot += 1
                        
            if i[10] == 1:
                self.ghosePassTot += 1
            if i[11] == 1:
                self.rule3Tot += 1
    
    def generateRO5Histogram(self): 
        """ Plot the number of candidates that are RO5 compliant based on the categories"""
        # create the dataset
        data = {'Partition Coefficient': self.clogTot,
               'H bond Acceptors': self.hAcceptorsTot,
               'H bond Donors': self.hDonorsTot,
               'Molecular Weight': self.molMassTot,
               'Number of Atoms': self.numAtomsTot,
               'Number of Rotatable Bonds': self.rotBondsTot,
               'Total Polar Surface Area': self.tpsaTot,
               'Lipophilicity': self.logPTot,
               'Molar Refractivity': self.molRefractTot}
        category = list(data.keys())
        values = list(data.values())
        
        fig = plt.figure(figsize = (10, 5))
        font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}
        
        #generate the histogram
        plt.bar(category, values, color ='blue', width = 0.7)
        
        plt.xlabel("Catergory of RO5 Compliance")
        plt.ylabel("Number of Drug Candidates")
        plt.title("Overall Assessment of Screening Library Based on RO5 Compliance")
        plt.savefig(self.name+"RO5Histogram.png")
        plt.show()

        
    def generateVariantHistogram(self): 
        """Plot the number of candidates that are RO3 Compliant and Ghose Filter Compliant"""
        # create the dataset
        data = {'Rule of 3 Compliant': self.rule3Tot,
               'Ghose Filter Compliant': self.ghosePassTot}
        category = list(data.keys())
        values = list(data.values())
        
        fig = plt.figure(figsize = (10, 5))
        
        #generate the histogram
        plt.bar(category, values, color ='pink', width = 0.7)
        
        plt.xlabel("Catergory of Variant Tests")
        plt.ylabel("Number of Drug Candidates")
        plt.title("Overall Assessment of Screening Library Based on Variant Compliance")
        plt.savefig(self.name+"VariantHistogram.png")
        plt.show()

    def generateDockingScorePlot(self): 
        """ Plot the range of the docking scores"""
        #create the dataset
        data = {'0 to -3.6': self.dockScoreZ1Tot,
               '-3.6 to -7.2': self.dockScoreZ2Tot,
               '-7.2 to -11': self.dockScoreZ3Tot}
        #generate the histogram
        category = list(data.keys())
        values = list(data.values())
        
        fig = plt.figure(figsize = (10, 5))
        plt.bar(category, values, color ='green', width = 0.7)
        font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}
        
        plt.xlabel("Range of Docking Scores")
        plt.ylabel("Number of Drug Candidates")
        plt.title("Overall Assesment of Docking Scores in Library")
        plt.savefig(self.name+'DockingScorePlot.png')
        plt.show()