from SDFREADERNEW import SdfReader

class sortData: 
    """ Sort through the input SDF file,
    output drug candidates based on their R05 variant compliance"""
    """ Take returned data from SDF reader class, 
    parse and sort data according to RO5 variant categories"""
    def __init__(self,infile): 
        """Initialize categories of RO5 variant compliance"""
        self.reader=SdfReader(infile) # I don't know if this is the appropriate way to call the file into the sdf reader class
        self.masterList=self.reader.readSDF()
        self.sortedinput=[]
        
        #not sure where to run the function so thought init would be good
        #self.checkvalidity()
    
    def checkvalidity(self):
        self.sortedinput=[]
        for i in self.masterList:
            temp = list(self.reader.returnVals(i)) #going through readSDF list of lists
            temp2 = list(self.reader.returnVals(i)) #for ghose filter and rule of 3 compliance
            print("list")
            print(i[2:])
            print("temp")
            print(temp)
            
            #we might need to change the undefined values for the attributes to None, otherwise it will mess up
            #checking clogP validity because the valid range includes 0. Also the heat map for dock score would
            #be skewed by any drug without a dock score as they would be set to 0 instead of invalid
            
            if temp[0] is not None:
                if (float(temp[0]) >= -0.4) and (float(temp[0]) <= 5.6): #check if clogP is valid
                    temp[0] = 1 #valid
                else:
                    temp[0] = 0
            else:
                temp[0] = 0 #invalid
                
            if temp[1] is not None:
                if (float(temp[1]) <= 3): #check for hAcceptor validity
                    temp[1] = 1
                else:
                    temp[1] = 0
            else:
                temp[1] = 0
            
            if temp[2] is not None:
                if (float(temp[2]) <= 3):
                    temp[2] = 1
                else:
                    temp[2] = 0
            else:
                temp[2] = 0
                
            if temp[3] is not None:
                if (float(temp[3]) >= 180) and (float(temp[3]) <= 480): #check for molecular weight validity
                    temp[3] = 1
                else:
                    temp[3] = 0
            else:
                temp[3] = 0
                
            if temp[4] is not None:
                if (float(temp[4]) >= 20) and (float(temp[4]) <= 70): #check for number of atoms validity 
                    temp[4] = 1
                else:
                    temp[4] = 0
            else:
                temp[4] = 0
                
            if temp[5] is not None:
                if (float(temp[5]) <= 3): #check rotatable bond value validity
                    temp[5] = 1
                else:
                    temp[5] = 0
            else:
                temp[5] = 0
            
            if temp[6] is not None:
                if (float(temp[6]) <= 140): #check tpsa validity
                    temp[6] = 1
                else:
                    temp[6] = 0
            else:
                temp[6] = 0
                
            if temp[7] is not None:#the comparison values need to be changed because the ghose filter is less specific than this filter
                if (float(temp[7]) >= -1) and (float(temp[7]) <= 6): #check lipo validity
                    temp[7] = 1
                else:
                    temp[7] = 0
            else:
                temp[7] = 0
                
            if temp[8] is not None:
                if (float(temp[8]) >= 40) and (float(temp[8]) <= 130): #check molecular refractivity validity
                    temp[8] = 1
                else:
                    temp[8] = 0
            else:
                temp[8] = 0
            
            if temp[9] is not None:
                if (float(temp[9]) <= -5.7): #check docking score validity and storing value for heat map
                    temp[9] = [1, temp[9]]
                else:
                    temp[9] = [0, temp[9]]
            else:
                temp[9] = [0, temp[9]]
            
            #for ghose and rule of 3 make sure to check if None first
            #check to see if it passses ghose filter
            if (temp2[3] is not None) and (temp2[7] is not None) and (temp2[4] is not None) and (temp2[8] is not None):
                if (float(temp2[3]) >= 160) and (float(temp2[3]) <= 480) and (float(temp2[7]) >= -0.4) and (float(temp2[7]) <= 5.6) and (float(temp2[4]) >= 20) and (float(temp2[4]) <=70) and (float(temp2[8]) >= 40) and (float(temp2[8]) <= 130):
                    temp.append(1)
                else:
                    temp.append(0)
            else:
                temp.append(0)
                
            #check for rule of 3    
            if (temp2[7] is not None) and (temp2[2] is not None) and (temp2[1] is not None) and (temp2[3] is not None) and (temp2[5] is not None):
                if (float(temp2[7]) <= 3) and (float(temp2[2]) <= 3) and (float(temp2[1]) <= 3) and (float(temp2[3]) < 300) and (float(temp2[5]) <= 3):
                    temp.append(1)
                else:
                    temp.append(0)
            else:
                temp.append(0)
            
            
            
            
            self.sortedinput.append(temp) #the valid/invalid list added to list of list with same index as list from SDFREADER
        return(self.sortedinput)