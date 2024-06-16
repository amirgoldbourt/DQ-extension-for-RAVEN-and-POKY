import tkinter
from tkinter import filedialog
import urllib.request
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin


#### Object Oriented Programming - defined class of amino acid (intended to allow the user add new features more conveniently) ####
class Amino_acid:
    def __init__(self, amino_acid_and_location):
        self.amino_acid_and_location=amino_acid_and_location
        amino_acid=""
        i=0
        while amino_acid_and_location[i].isalpha():
            amino_acid+=amino_acid_and_location[i]
            i+=1
        self.AA_identity = amino_acid
        self.AA_Location = str(amino_acid_and_location[i:])

    def AA_ID(self):
        if self.AA_identity == "ALA" or  self.AA_identity=="Ala" or self.AA_identity=="A":
            return "A"
        if self.AA_identity == "CYS" or  self.AA_identity=="Cys" or self.AA_identity=="C":
            return "C"
        if self.AA_identity == "ASP" or  self.AA_identity=="Asp" or self.AA_identity=="D":
            return "D"
        if self.AA_identity == "GLU" or  self.AA_identity=="Glu" or self.AA_identity=="E":
            return "E"
        if self.AA_identity == "PHE" or  self.AA_identity=="Phe" or self.AA_identity=="F":
            return "F"
        if self.AA_identity == "GLY" or  self.AA_identity=="Gly" or self.AA_identity=="G":
            return "G"
        if self.AA_identity == "HIS" or  self.AA_identity=="His" or self.AA_identity=="H":
            return "H"
        if self.AA_identity == "ILE" or  self.AA_identity=="Ile" or self.AA_identity=="I":
            return "I"
        if self.AA_identity == "LYS" or  self.AA_identity=="Lys" or self.AA_identity=="K":
            return "K"
        if self.AA_identity == "LEU" or  self.AA_identity=="Leu" or self.AA_identity=="L":
            return "L"
        if self.AA_identity == "MET" or  self.AA_identity=="Met" or self.AA_identity=="M":
            return "M"
        if self.AA_identity == "ASN" or  self.AA_identity=="Asn" or self.AA_identity=="N":
            return "N"
        if self.AA_identity == "PRO" or  self.AA_identity=="Pro" or self.AA_identity=="P":
            return "P"
        if self.AA_identity == "GLN" or  self.AA_identity=="Gln" or self.AA_identity=="Q":
            return "Q"
        if self.AA_identity == "ARG" or  self.AA_identity=="Arg" or self.AA_identity=="R":
            return "R"
        if self.AA_identity == "SER" or  self.AA_identity=="Ser" or self.AA_identity=="S":
            return "S"
        if self.AA_identity == "THR" or  self.AA_identity=="Thr" or self.AA_identity=="T":
            return "T"
        if self.AA_identity == "VAL" or  self.AA_identity=="Val" or self.AA_identity=="V":
            return "V"
        if self.AA_identity == "TRP" or  self.AA_identity=="Trp" or self.AA_identity=="W":
            return "W"
        if self.AA_identity == "TYR" or  self.AA_identity=="Tyr" or self.AA_identity=="Y":
            return "Y"
        if self.AA_Location == None or self.AA_Location == "-" or self.AA_Location == " ":
            return "-"

    def __repr__(self):
        return f"{self.AA_ID()}{self.AA_Location}"

    def chemical_property(self):
        if self.AA_ID() in ["G","A","P","V","L","I","M"]:
            return "Nonpolar aliphatic"
        if self.AA_ID() in ["F","Y","W"]:
            return "Aromatic"
        if self.AA_ID() in ["S","T","C","N","Q"]:
            return "Polar Uncharged"
        if self.AA_ID() in ["K","H","R"]:
            return "Positively Charged"
        if self.AA_ID() in ["D","E"]:
            return "Negatively Charged"

    def __eq__(self, other):
        if isinstance(other, Amino_acid):
            if f"{self.AA_identity}{self.AA_Location}" == f"{other.AA_identity}{other.AA_Location}":
                return True
        return False


# STAGE 1- get an assignment file in RAVEN format and convert it into an assignment table suitable for INADEQUATE experiment processing:
# step 1- convert the assignment file into a matrix
# step 2- create a dictionary of the possible atoms (C,Ca,Cb,....) and their associated index
# step 3- functions that create new line according to INADEQUATE saving data technique (need the number of amino acid in the protein)
# step 4- function that returns new matrix after adjustments according to the chosen function in step 4
# step 5- function that merge matrices form step 1 (the original) and from step 4 (after adjustment) and write it into a new file (destination file is required)
# step 6- function that generalize all the above functions into one. also adjust step 1 in case there is no assignment to an amino acid.


##########################################################################################################################################
##########################################################################################################################################

# function that receives an assignment file and converts it into a matrix
# delete empty line
def get_Assignment_str_matrix(file_name):
    f=open(file_name)
    list=f.readlines()
    f.close()
    row_number=len(list)
    mat=[]
    for i in range (row_number):
        mat.append([None])

    for j in range (row_number):
        list[j]=list[j].replace("\t"," ")
        list[j]=list[j].rstrip(" \n")
        while "  " in list[j]:
            list[j] = list[j].replace("  ", " ")
        list[j]=list[j].title()
        if j>=1:
            list[j]=list[j].strip()
        mat[j]=list[j].split(" ")

    for k in range (1,row_number):
        mat[k][0]=Amino_acid(mat[k][0])

    for m in range (len(mat[0])):
        mat[0][m]=mat[0][m].title()
    return mat


# function that creates a dictionary whose key is the name of the atom and its value is the corresponding index (for the first line)
def get_Dictionary(mat):
    D={}
    for i in range(len(mat[0])):
        D[mat[0][i]]=mat[0].index(mat[0][i])
    return D


# list of functions that get an assignment line and return a new line that saves data according to possible J-based INADEQUATE criteria
def A_adjustment (o_l,D,num_of_AA):  #o_l=original_list
    list=o_l[:]
    list[D.get("Cb")]="-"
    if o_l[D.get("Ca")] not in ["*","-"] and o_l[D.get("Cb")] not in ["*","-"] :
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")])+float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")] not in ["*","-"]  and o_l[D.get("Ca")] not in ["*","-"] :
        list[D.get("C")]=str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def C_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cb")] = "-"
    if o_l[D.get("Cb")] not in ["*","-"]  and o_l[D.get("Cb")] not in ["*","-"] :
        list[D.get("Cb")] = str(round(float(o_l[D.get("Ca")]) + float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")] = "-"
    if o_l[D.get("C")] not in ["*","-"]  and o_l[D.get("Ca")] not in ["*","-"] :
        list[D.get("C")] = str(round(float(o_l[D.get("C")]) + float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")] = "-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def D_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cg")]="-"
    if o_l[D.get("Cb")] not in ["*","-"]  and o_l[D.get("Cg")] not in ["*","-"] :
        list[D.get("Cb")]=str(round(float(o_l[D.get("Cb")])+float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")] not in ["*","-"]  and o_l[D.get("Cb")] not in ["*","-"] :
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")])+float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")] not in ["*","-"]  and o_l[D.get("Ca")] not in ["*","-"] :
        list[D.get("C")] =str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def E_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cd")]="-"
    if o_l[D.get("Cg")] not in ["*","-"] and o_l[D.get("Cd")]not in ["*","-"] :
        list[D.get("Cg")]=str(round(float(o_l[D.get("Cg")])+float(o_l[D.get("Cd")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")] not in ["*","-"] and o_l[D.get("Cg")] not in ["*","-"] :
        list[D.get("Cb")]=str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"] :
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")]) + float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"] :
        list[D.get("C")] =str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def F_adjustment (o_l,D,num_of_AA): #In rings it starts with the carbon with the "highest" letter, and going clockwise until completing the ring
    list = o_l[:]
    if o_l[D.get("Cz")]not in ["*","-"] and o_l[D.get("Ce2")]not in ["*","-"] :
        list[D.get("Cz")]=str(round(float(o_l[D.get("Cz")])+float(o_l[D.get("Ce2")]),2))
    else:
        list[D.get("Cz")]="-"
    if o_l[D.get("Ce2")]not in ["*","-"] and o_l[D.get("Cd2")]not in ["*","-"] :
        list[D.get("Ce2")] = str(round(float(o_l[D.get("Ce2")]) + float(o_l[D.get("Cd2")]),2))
    else:
        list[D.get("Ce2")]="-"
    if o_l[D.get("Cd2")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"] :
        list[D.get("Cd2")] = str(round(float(o_l[D.get("Cd2")]) + float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cd2")]="-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd1")]not in ["*","-"] :
        list[D.get("Cg")] = str(round(float(o_l[D.get("Cg")])+float(o_l[D.get("Cd1")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cd1")]not in ["*","-"] and o_l[D.get("Ce1")]not in ["*","-"] :
        list[D.get("Cd1")] = str(round(float(o_l[D.get("Cd1")])+ float(o_l[D.get("Ce1")]),2))
    else:
        list[D.get("Cd1")]="-"
    if o_l[D.get("Ce1")]not in ["*","-"] and o_l[D.get("Cz")]not in ["*","-"] :
        list[D.get("Ce1")]=str(round(float(o_l[D.get("Ce1")])+ float(o_l[D.get("Cz")]),2))
    else:
        list[D.get("Ce1")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"] :
        list[D.get("Cb")]=str(round(float(o_l[D.get("Cb")])+float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"] :
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")])+float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"] :
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def G_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"] :
        list[D.get("C")]=str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def H_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Ce1")] = "-"
    list[D.get("Cd2")] = "-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd2")]not in ["*","-"]:
        list[D.get("Cg")] = str(round(float(o_l[D.get("Cg")])+float(o_l[D.get("Cd2")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")]=str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def I_adjustment (o_l,D,num_of_AA): 
    list = o_l[:]
    if o_l[D.get("Cg2")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Cg2")] = str(round(float(o_l[D.get("Cg2")]) + float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Cg2")]="-"
    list[D.get("Cd1")] = "-"
    if o_l[D.get("Cg1")] not in ["*","-"] and o_l[D.get("Cd1")] not in ["*","-"]:
        list[D.get("Cg1")] = str(round(float(o_l[D.get("Cg1")])+float(o_l[D.get("Cd1")]),2))
    else:
        list[D.get("Cg1")]="-"
    if o_l[D.get("Cb")] not in ["*","-"] and o_l[D.get("Cg1")] not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg1")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")] not in ["*","-"] and o_l[D.get("Cb")] not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")] not in ["*","-"] and o_l[D.get("Ca")] not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def K_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Ce")]="-"
    if o_l[D.get("Cd")] not in ["*","-"] and o_l[D.get("Ce")]not in ["*","-"]:
        list[D.get("Cd")]=str(round(float(o_l[D.get("Cd")])+float(o_l[D.get("Ce")]),2))
    else:
        list[D.get("Cd")]="-"
    if o_l[D.get("Cg")]not in ["*","-"]  and o_l[D.get("Cd")]not in ["*","-"]:
        list[D.get("Cg")]=str(round(float(o_l[D.get("Cg")])+float(o_l[D.get("Cd")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"] :
        list[D.get("Cb")]=str(round(float(o_l[D.get("Cb")])+float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"]  and o_l[D.get("Cb")]not in ["*","-"] :
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"] :
        list[D.get("C")] =str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def L_adjustment (o_l,D,num_of_AA): 
    list = o_l[:]
    if o_l[D.get("Cd2")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"] :
        list[D.get("Cd2")] = str(round(float(o_l[D.get("Cd2")]) + float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cd2")]="-"
    list[D.get("Cd1")] = "-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd1")]not in ["*","-"]:
        list[D.get("Cg")] = str(round(float(o_l[D.get("Cg")])+ float(o_l[D.get("Cd1")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def M_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Ce")]="-"
    list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")]=str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")]=str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] =str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def N_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cg")] = "-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="- "
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def P_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cd")] = "-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd")]not in ["*","-"]:
        list[D.get("Cg")] = str(round(float(o_l[D.get("Cg")]) + float(o_l[D.get("Cd")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] =str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] =str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def Q_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cd")] = "-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd")]not in ["*","-"]:
        list[D.get("Cg")] =str(round(float(o_l[D.get("Cg")])+ float(o_l[D.get("Cd")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] =str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] =str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def R_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cz")]= "-"
    list[D.get("Cd")]= "-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd")]not in ["*","-"]:
        list[D.get("Cg")]= str(round(float(o_l[D.get("Cg")])+ float(o_l[D.get("Cd")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")]= str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")]= str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"] :
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def S_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cb")]= "-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")]= str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def T_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    list[D.get("Cg2")]= "-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg2")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg2")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def V_adjustment (o_l,D,num_of_AA): 
    list = o_l[:]
    if o_l[D.get("Cg2")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Cg2")] = str(round(float(o_l[D.get("Cg2")]) + float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Cg2")]="-"
    list[D.get("Cg1")] = "-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg1")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")]) +float(o_l[D.get("Cg1")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def W_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    if o_l[D.get("Ch2")]not in ["*","-"] and o_l[D.get("Cz2")]not in ["*","-"]:
        list[D.get("Ch2")] = str(round(float(o_l[D.get("Ch2")]) + float(o_l[D.get("Cz2")]),2))
    else:
        list[D.get("Ch2")]="-"
    if o_l[D.get("Cz2")]not in ["*","-"] and o_l[D.get("Ce2")]not in ["*","-"]:
        list[D.get("Cz2")] = str(round(float(o_l[D.get("Cz2")])+ float(o_l[D.get("Ce2")]),2))
    else:
        list[D.get("Cz2")]="-"
    if o_l[D.get("Ce2")]not in ["*","-"] and o_l[D.get("Cd2")]not in ["*","-"]:
        list[D.get("Ce2")] = str(round(float(o_l[D.get("Ce2")]) + float(o_l[D.get("Cd2")]),2))
    else:
        list[D.get("Ce2")]="-"
    if o_l[D.get("Cd2")]not in ["*","-"] and o_l[D.get("Ce3")]not in ["*","-"]:
        list[D.get("Cd2")] = str(round(float(o_l[D.get("Cd2")])+ float(o_l[D.get("Ce3")]),2))
    else:
        list[D.get("Cd2")]="-"
    if o_l[D.get("Ce3")]not in ["*","-"] and o_l[D.get("Cz3")]not in ["*","-"]:
        list[D.get("Ce3")] = str(round(float(o_l[D.get("Ce3")])+ float(o_l[D.get("Cz3")]),2))
    else:
        list[D.get("Ce3")]="-"
    if o_l[D.get("Cz3")]not in ["*","-"] and o_l[D.get("Ch2")]not in ["*","-"]:
        list[D.get("Cz3")] = str(round(float(o_l[D.get("Cz3")])+ float(o_l[D.get("Ch2")]),2))
    else:
        list[D.get("Cz3")]="-"
    if o_l[D.get("Cd1")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cd1")] = str(round(float(o_l[D.get("Cd1")]) + float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cd1")]="-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd2")]not in ["*","-"]:
        list[D.get("Cg")] = str(round(float(o_l[D.get("Cg")])+ float(o_l[D.get("Cd2")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list

def Y_adjustment (o_l,D,num_of_AA):
    list = o_l[:]
    if o_l[D.get("Cz")]not in ["*","-"] and o_l[D.get("Ce2")]not in ["*","-"]:
        list[D.get("Cz")] = str(round(float(o_l[D.get("Cz")]) + float(o_l[D.get("Ce2")]),2))
    else:
        list[D.get("Cz")]="-"
    if o_l[D.get("Ce2")]not in ["*","-"] and o_l[D.get("Cd2")]not in ["*","-"]:
        list[D.get("Ce2")] = str(round(float(o_l[D.get("Ce2")]) + float(o_l[D.get("Cd2")]),2))
    else:
        list[D.get("Ce2")]="-"
    if o_l[D.get("Cd2")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cd2")] = str(round(float(o_l[D.get("Cd2")]) + float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cd2")]="-"
    if o_l[D.get("Cg")]not in ["*","-"] and o_l[D.get("Cd1")]not in ["*","-"]:
        list[D.get("Cg")] = str(round(float(o_l[D.get("Cg")])+ float(o_l[D.get("Cd1")]),2))
    else:
        list[D.get("Cg")]="-"
    if o_l[D.get("Cd1")]not in ["*","-"] and o_l[D.get("Ce1")]not in ["*","-"]:
        list[D.get("Cd1")] = str(round(float(o_l[D.get("Cd1")]) + float(o_l[D.get("Ce1")]),2))
    else:
        list[D.get("Cd1")]="-"
    if o_l[D.get("Ce1")]not in ["*","-"] and o_l[D.get("Cz")]not in ["*","-"]:
        list[D.get("Ce1")] = str(round(float(o_l[D.get("Ce1")])+float(o_l[D.get("Cz")]),2))
    else:
        list[D.get("Ce1")]="-"
    if o_l[D.get("Cb")]not in ["*","-"] and o_l[D.get("Cg")]not in ["*","-"]:
        list[D.get("Cb")] = str(round(float(o_l[D.get("Cb")])+ float(o_l[D.get("Cg")]),2))
    else:
        list[D.get("Cb")]="-"
    if o_l[D.get("Ca")]not in ["*","-"] and o_l[D.get("Cb")]not in ["*","-"]:
        list[D.get("Ca")] = str(round(float(o_l[D.get("Ca")])+ float(o_l[D.get("Cb")]),2))
    else:
        list[D.get("Ca")]="-"
    if o_l[D.get("C")]not in ["*","-"] and o_l[D.get("Ca")]not in ["*","-"]:
        list[D.get("C")] = str(round(float(o_l[D.get("C")])+ float(o_l[D.get("Ca")]),2))
    else:
        list[D.get("C")]="-"
    list[0]=f"&{o_l[0].AA_ID()}{int(o_l[0].AA_Location)+num_of_AA}"
    return list


# function that gets an assignment matrix and returns a new matrix after the data adjustments
def make_adjustment(mat, D, num_of_AA):
    new_mat=mat[:][:]
    for i in range(1,len(mat)):
        if mat[i][0].AA_ID() == "A":
            new_mat[i] = A_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "C":
            new_mat[i] = C_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "D":
            new_mat[i] = D_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "E":
            new_mat[i] = E_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "F":
            new_mat[i] = F_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "G":
            new_mat[i] = G_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "H":
            new_mat[i] = H_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "I":
            new_mat[i] = I_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "K":
            new_mat[i] = K_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "L":
            new_mat[i] = L_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "M":
            new_mat[i] = M_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "N":
            new_mat[i] = N_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "P":
            new_mat[i] = P_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "Q":
            new_mat[i] = Q_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "R":
            new_mat[i] = R_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "S":
            new_mat[i] = S_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "T":
            new_mat[i] = T_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "V":
            new_mat[i] = V_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "W":
            new_mat[i] = W_adjustment(mat[i], D,num_of_AA)
        elif mat[i][0].AA_ID() == "Y":
            new_mat[i] = Y_adjustment(mat[i], D,num_of_AA)
    return new_mat


# function that creates a new assignment file that contains both original and INADEQUATE data.
def write_merge_file (flie_path,original_mat,adjust_mat):
    combine_mat=original_mat+adjust_mat[1:]
    text=""
    for i in range(len(combine_mat)):
        for j in range (len(combine_mat[0])-1):
            text+=f"{combine_mat[i][j]}\t"
            if j==0:
                text+="\t"
        text+=combine_mat[i][-1]+"\n"
    f=open(flie_path,"w")
    f.write(text)
    return


# function that gets a sequence and an assignment file and creates two new files suitable for RAVEN.
def adjustment_assignment_table_and_seq (shift_assignment_fill,seq_fil,path):
    # write a new sequence file with double the amino acids
    f = open(seq_fil, "r")
    txt = f.read()
    txt = txt.strip("\n")
    new_str=""
    for i in txt:
        if i!=" ":
            new_str+=i
    step0 = len(new_str)
    new_str+=new_str+"\n"
    seq_name=""
    i=-5
    while seq_fil[i]!="/":
        seq_name+=seq_fil[i]
        i-=1
    f_new = open(f"{path}/{seq_name[::-1]} duplicated.seq", "w")
    f_new.write(new_str)
    # get the assignment matrix
    step0_1=get_Assignment_str_matrix(shift_assignment_fill)
    # adjust the matrix: if there is no amino acid assignment then change the "-" sign to the supposed amino acid
    step1=[step0_1[0]]
    i=1
    j=0
    while i<len(step0_1):

        a = step0_1[i][0].AA_ID()
        b = new_str[j]

        if step0_1[i][0].AA_ID() =="-":
            step1 += [step0_1[i]]
            step1[-1][0]=Amino_acid(new_str[j]+str(j+1))
            j += 1
            i += 1
        elif step0_1[i][0].AA_ID()!=new_str[j]:
            add_line=["-"]*len(step1[0])
            add_line[0]=new_str[j]+str(j+1)
            add_line[0]=Amino_acid(add_line[0])
            step1+=[add_line]
            j+=1
        elif step0_1[i][0].AA_ID()==new_str[j]:
            step1+=[step0_1[i]]
            i+=1
            j+=1
    # in case that there is an unassigned amino acid at the end of the sequence and there is no "-" sign in the assignment table, add a new line to the matrix with the corresponding amino acid
    while j<step0:
        add_line = ["-"] * len(step1[0])
        add_line[0] = new_str[j] + str(j + 1)
        add_line[0] = Amino_acid(add_line[0])
        step1 += [add_line]
        j += 1
    # continue with the process and create a new assignment file according to J-based INADEQUATE criteria
    step2=get_Dictionary(step1)
    step3=make_adjustment(step1,step2,step0)
    table_name = ""
    i = -8
    while shift_assignment_fill[i] != "/":
        table_name += shift_assignment_fill[i]
        i -= 1
    step4 = write_merge_file(f"{path}/{table_name[::-1]} DQ adjastment.shifts", step1,step3)
    return (f"{path}/{table_name[::-1]} DQ adjastment.shifts",f"{path}/{seq_name[::-1]} duplicated.seq")


######################################################################################################################################################
######################################################################################################################################################

# STAGE 2- getting a full match file and reduce it to the possible correlations that can appear in a J-based INADEQUATE experiment. change the name of the stored data to the correct one.
# step 1- convert the match file into a matrix
# step 2- delete different amino acid correlation
# step 3- functions that change the name of carbon to the two carbons correlation it represents (for each amino acid)
# step 4- function that gets a matrix and changes the name of the represented INADEQUATE experiment data carbon into the correct  carbons-couple correlation
# step 5- functions that modify the match data according to isotopic labeling scheme
# step 6- function that deletes peaks that are unlikely to appear (with lower ILP than criterion)
# step 7- function that changes the overall value according to the new mat
# step 8- functions that sort the match file according to the amino acid sequence or according to the value in the peak list
# step 9- function that combines all the above function into one.
# get match file, sequence file and isotopic labeling scheme and create a new match file that contains only valid J-based INADEQUATE options after fitting to the isotopic labeling scheme


# convert match file into a matrix
def matches_as_matrix(file_location):
    f=open(file_location)
    lists=f.readlines()
    f.close()
    mat=[]
    for i in range(len(lists)):
        mat.append([None])
    for j in range (len(lists)):
        mat[j]=lists[j].split("\t")
    return mat


# functions that get the number of digits in the length of the matched carbon location
def get_first_aa_len (lst):
    firs_aa_len = 0
    j = 1
    while lst[1][j].isdigit():
        firs_aa_len += 1
        j += 1
    return firs_aa_len

def get_second_aa_len (lst):
    separate_index = lst[1].index("-")
    k = separate_index + 2
    second_aa_len = 0
    while lst[1][k].isdigit():
        second_aa_len += 1
        k += 1
    return  second_aa_len


# functions that get a line in a match matrix and return the name of the associated carbons
def get_first_carbon_ID(lst):
    firs_aa_len = get_first_aa_len(lst)
    carbon1 = ""
    n = 0
    while lst[1][1 + firs_aa_len + n] != "(":
        carbon1 += lst[1][1 + firs_aa_len + n]
        n += 1
    return carbon1

def get_sec_carbon_ID(lst):
    separate_index = lst[1].index("-")
    second_aa_len = get_second_aa_len(lst)
    n = 0
    carbon2 = ""
    while lst[1][separate_index + 2 + second_aa_len + n] != "(":
        carbon2 += lst[1][separate_index + 2 + second_aa_len + n]
        n += 1
    return carbon2


# delete different amino acid correlations
def same_aa_correlation_mat (mat,num_of_AA):
    new_mat=[]
    new_mat.append(mat[0])

    for i in range(1,len(mat)):
        if mat[i]!=['\n']:
            # finding the number of the index that is needed for the location of the first amino acid in the match
            firs_aa_len=get_first_aa_len(mat[i])
            separd_index=mat[i][1].index("-")
            secend_aa_len=get_second_aa_len(mat[i])
                # add to the new matrix only if the correlation is between two carbons in the same amino acid
            if mat[i][1][0]==mat[i][1][separd_index+1]:
                # add to the new matrix only if the correlation is between two carbons in the same amino acid
                    if int(mat[i][1][1:firs_aa_len+1])==int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len])+num_of_AA :
                        new_mat.append(mat[i])
    return new_mat



# list of functions that change the name of a carbon to the two carbons correlation it represents (for each amino acid)
def A_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    return str

def C_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    return str

def D_rename(str):
    if str == "C":
        str = "CCa"
    elif str == "Ca":
        str = "CaCb"
    elif str=="Cb":
        str="CbCg"
    return str

def E_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd"
    return str

def F_rename(str):
    if str== "C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd1"
    elif str=="Cd1":
        str="Cd1de1"
    elif str=="Ce1":
        str="Ce1Cz"
    elif str=="Cz":
        str="CzCe2"
    elif str=="Ce2":
        str="Ce2Cd2"
    elif str=="Cd2":
        str="Cd2Cg"
    return str

def G_rename (str):
    if str=="C":
        str="CCa"
    return str

def H_rename(str):
    if str== "C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd2"
    return str

def I_rename(str):
    if str== "C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg1"
    elif str=="Cg1":
        str="Cg1Cd1"
    elif str=="Cg2":
        str="CbCg2"
    return str

def K_rename (str):
    if str== "C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd"
    elif str=="Cd":
        str="CdCe"
    return str

def L_rename (str):
    if str== "C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd1"
    elif str=="Cd2":
        str="CgCd2"
    return str

def M_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    return str

def N_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    return str

def P_rename (str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd"
    return str

def Q_rename (str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd"
    return str

def R_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd"
    return str

def S_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    return str

def T_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg2"
    return str

def V_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg1"
    elif str=="Cg2":
        str="CbCg2"
    return str

def W_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd2"
    elif str== "Cd1":
        str="CgCd1"
    elif str=="Cd2":
        str="Cd2Ce3"
    elif str=="Ce3":
        str="Ce3Cz3"
    elif str=="Cz3":
        str="Cz3Ch2"
    elif str== "Ch2":
        str="Cz2Ch2"
    elif str=="Cz2":
        str="Ce2Cz2"
    elif str=="Ce2":
        str="Cd2Ce2"
    return  str

def Y_rename(str):
    if str=="C":
        str="CCa"
    elif str=="Ca":
        str="CaCb"
    elif str=="Cb":
        str="CbCg"
    elif str=="Cg":
        str="CgCd1"
    elif str=="Cd1":
        str="Cd1Ce1"
    elif str=="Ce1":
        str="Ce1Cz"
    elif str=="Cz":
        str="CzCe2"
    elif str=="Ce2":
        str="Ce2Cd2"
    elif str=="Cd2":
        str="CgCd2"
    return str


# function that gets a matrix with data represented by the INADEQUATE data storage technique and transforms it
def rename_mtches_and_short(mat):
    new_mat=[]
    new_mat.append(mat[0])
    for i in range(1,len(mat)):
        separd_index = mat[i][1].index("-")
        firs_aa_len = get_first_aa_len(mat[i])
        secend_aa_len = get_second_aa_len(mat[i])
        carbon1= get_first_carbon_ID(mat[i])
        carbon2= get_sec_carbon_ID(mat[i])
        # delete any interaction that does not exist between two adjacent carbons
        adjacent_carbon_options=[["a","b"],["b","g"],["b","g1"],["b","g2"],["g1","d1"], ["g","d"],["g","d1"] , ["g","d2"], ["d","e"] , ["d1","e1"], ["d2","e2"], ["d2","e3"], ["e","z"] , ["e1","z"], ["e2","z"], ["e2","z2"], ["e3","z3"], ["e2","z"], ["z","h"], ["z2","h2"], ["z3","h2"]]
        if (carbon1.title() in ["C","Ca"] and carbon2.title() in ["C","Ca"] ) or (carbon1.upper()!="C" and carbon2!="C" and  any(all(pair in option for pair in [carbon1.title()[1:],carbon2.title()[1:]])for option in adjacent_carbon_options)):
        # making the adjustments for the supposed carbon between the two
            if mat[i][1][0]=="A":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+A_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+A_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0]=="C":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+C_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+C_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "D":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+D_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+D_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "E":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+E_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+E_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "F":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+F_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+F_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "G":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+G_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+G_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "H":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+H_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+H_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "I":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+I_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+I_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "K":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+K_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+K_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "L":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+L_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+L_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "M":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+M_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+M_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "N":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+N_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+N_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "P":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+P_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+P_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "Q":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                  str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+Q_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                  #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+Q_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "R":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+R_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+R_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "S":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+S_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+S_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "T":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+T_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+T_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "V":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+V_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+V_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0] == "W":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+W_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+W_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            elif mat[i][1][0]=="Y":
                if int(mat[i][1][1:1+firs_aa_len])>int(mat[i][1][separd_index+2:separd_index+2+secend_aa_len]):
                    str=mat[i][1][0]+mat[i][1][separd_index+2:separd_index+2+secend_aa_len]+Y_rename(carbon1.title())+mat[i][1][1+firs_aa_len+len(carbon1):]
                #else:
                    #str=mat[i][1][0:separd_index+2]+mat[i][1][1:1+firs_aa_len]+Y_rename(carbon2.title())+mat[i][1][separd_index+2+secend_aa_len+len(carbon2):]

            new_mat.append(mat[i])
            new_mat[-1][1]=str
    return  new_mat





# Functions that delete information that does not fit the isotopic labeling scheme.
# The isotopic labeling probability according to:
# V.A. Higman, J. Flinders, M. Hiller, S. Jehle, S. Markovic, S. Fiedler, B.J. van Rossum, H. Oschkinat-
# Assigning large proteins in the solid state: a MAS NMR resonance assignment strategy using selectively and extensively 13C-labelled proteins
# https://www.sciencedirect.com/science/article/pii/S0079656518300116#b0620
def uniform_adj(mat):
    uniform_labeling_mat = [["aa", "C", "Ca", "Cb", "Cg", "Cg1", "Cg2", "Cd", "Cd1", "Cd2", "Ce", "Ce1", "Ce2", "Ce3", "Cz", "Cz2", "Cz3","Ch2", "N", "Nd1", "Nd2", "Ne", "Ne1", "Ne2", "Nz", "Nh1", "Nh2"],
                            ["A", 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["R", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,0.0, 0.0, 0.0, 1.0, 1.0],
                            ["D", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["N", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["C", 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["E", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["Q", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 1.0, 0.0, 0.0, 0.0],
                            ["G", 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["H", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,0.0, 1.0, 0.0, 0.0, 0.0],
                            ["I", 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["L", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["K", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 1.0, 0.0, 0.0],
                            ["M", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["F", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["P", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["S", 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["T", 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["W", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,1.0, 0.0, 0.0, 0.0, 0.0],
                            ["Y", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0],
                            ["V", 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0]]
    aa_option = ["A", "R", "D", "N", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    new_mat = []
    new_mat.append(mat[0])
    new_mat[0][-1] = new_mat[0][-1][:-1]
    new_mat[0].append("ILP\n")

    for i in range(1,len(mat)):
        new_mat.append(mat[i])
        new_mat[i][-1]=new_mat[i][-1][:-1]
        aa_in_mat = aa_option.index(new_mat[i][1][0])+1
        carbon_1_in_mat ="C"
        n = new_mat[i][1].find("C")+1
        while new_mat[i][1][n] != "C":
            carbon_1_in_mat += new_mat[i][1][n]
            n += 1
        carbon_2_in_mat =""
        while new_mat[i][1][n] != "(":
            carbon_2_in_mat+=new_mat[i][1][n]
            n+=1
        carbon_1_id=uniform_labeling_mat[0].index(carbon_1_in_mat)
        carbon_2_id=uniform_labeling_mat[0].index(carbon_2_in_mat)
        probability=round(float(uniform_labeling_mat[aa_in_mat][carbon_1_id])*float(uniform_labeling_mat[aa_in_mat][carbon_2_id]),2)
        new_mat[i].append(str(probability))
        new_mat[i][-1]=new_mat[i][-1]+"\n"
    new_mat.append(["\n"])
    return new_mat



def isotopic_1_3_Glycerol_adj (mat):
    ###################################################
    ###################################################
    ###################################################
    ###################################################




    glycerol_1_3_labeling_mat = [["aa", "C", "Ca", "Cb", "Cg","Cg1","Cg2", "Cd","Cd1","Cd2", "Ce","Ce1","Ce2","Ce3", "Cz","Cz2","Cz3","Ch2", "N", "Nd1", "Nd2", "Ne", "Ne1", "Ne2", "Nz", "Nh1", "Nh2"],
                                 ["A", 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["R", 0.40, 0.80, 0.60, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0],
                                 ["D", 0.50, 0.67, 0.83, 0.33, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["N", 0.50, 0.67, 0.83, 0.33, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["C", 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["E", 0.40, 0.80, 0.60, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["Q", 0.40, 0.80, 0.60, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                                 ["G", 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["H", 1.00, 0.00, 1.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                                 ["I", 0.50, 0.67, 0.00, 0.00, 0.83, 1.00, 0.00, 0.33, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["L", 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["K", 0.75, 0.42, 0.92, 0.33, 0.00, 0.00, 0.83, 0.00, 0.00, 0.42, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                                 ["M", 0.50, 0.67, 0.83, 0.33, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["F", 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.50, 0.50, 0.00, 1.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["P", 0.40, 0.80, 0.60, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["S", 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["T", 0.50, 0.67, 0.83, 0.00, 0.00, 0.33, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["W", 1.00, 0.00, 1.00, 0.50, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 1.00, 0.00, 1.00, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                                 ["Y", 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.50, 0.50, 0.00, 1.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                 ["V", 1.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

    aa_option=["A","R","D","N","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    new_mat=[]
    new_mat.append(mat[0])
    new_mat[0][-1]=new_mat[0][-1][:-1]
    new_mat[0].append("ILP\n")
    for i in range (1,len(mat)):
        aa_type=mat[i][1][0]
        correlation_id = ""
        n=mat[i][1].find("C")
        while mat[i][1][n]!="(":
            correlation_id+=mat[i][1][n]
            n+=1
        #continue only if the correlation fits the chosen labeling scheme
        if aa_type in ["F","Y"]:
            if correlation_id in ["Cd1Ce1","Ce2Cd2","Ce1Cz","CzCe2"]:
                new_mat.append(mat[i])
                new_mat[-1][-1]=new_mat[-1][-1][:-1]
                new_mat[-1].append("1.00\n")
        elif aa_type== "W":
            if correlation_id in ["CbCg","CgCd1"]:
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.50\n")
            elif correlation_id in ["Ce2Cz2","Cz2Ch2"]:
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("1.00\n")
        elif aa_type=="H":
            if correlation_id == "CbCg":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.50\n")
            elif correlation_id== "CgCd2":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.25\n")
        elif aa_type in ["R","Q","E","P"]:
            if correlation_id == "CCa":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.20\n")
            elif correlation_id== "CaCb":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.40\n")
            elif correlation_id=="CbCg":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.60\n")
        elif aa_type in ["N","D","M","T"]:
            if correlation_id == "CCa":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.17\n")
            elif correlation_id == "CaCb":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.50\n")
            elif correlation_id in ["CbCg","CbCg2"] :
                    new_mat.append(mat[i])
                    new_mat[-1][-1] = new_mat[-1][-1][:-1]
                    new_mat[-1].append("0.17\n")
        elif aa_type == "K":
            if correlation_id in ["CCa","CgCd"]:
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.17\n")
            elif correlation_id == "CaCb":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.33\n")
            elif correlation_id in ["CbCg","CdCe"]:
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.25\n")
        elif aa_type == "I":
            if correlation_id in ["CCa","Cg1Cd1"]:
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.17\n")
    new_mat.append(["\n"])
    return new_mat


def isotopic_2_Glycerol_adj (mat):
    glycerol_2_labeling_mat = [["aa", "C", "Ca", "Cb", "Cg","Cg1","Cg2", "Cd","Cd1","Cd2", "Ce","Ce1","Ce2","Ce3", "Cz","Cz2","Cz3","Ch2", "N", "Nd1", "Nd2", "Ne", "Ne1", "Ne2", "Nz", "Nh1", "Nh2"],
                               ["A", 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["R", 0.60, 0.20, 0.40, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0],
                               ["D", 0.50, 0.33, 0.17, 0.67, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["N", 0.50, 0.33, 0.17, 0.67, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["C", 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["E", 0.60, 0.20, 0.40, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["Q", 0.60, 0.20, 0.40, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                               ["G", 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["H", 0.00, 1.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                               ["I", 0.50, 0.33, 1.00, 0.00, 0.17, 0.00, 0.00, 0.67, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["L", 1.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["K", 0.25, 0.58, 0.08, 0.67, 0.00, 0.00, 0.17, 0.00, 0.00, 0.58, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                               ["M", 0.50, 0.33, 0.17, 0.67, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["F", 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["P", 0.60, 0.20, 0.40, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["S", 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["T", 0.50, 0.33, 0.17, 0.00, 0.00, 0.67, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["W", 0.00, 1.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 1.0,0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                               ["Y", 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               ["V", 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

    aa_option=["A","R","D","N","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    new_mat = []
    new_mat.append(mat[0])
    new_mat[0][-1] = new_mat[0][-1][:-1]
    new_mat[0].append("ILP\n")
    for i in range(1,len(mat)):
        aa_type = mat[i][1][0]
        correlation_id = ""
        n = mat[i][1].find("C")
        while mat[i][1][n] != "(":
            correlation_id += mat[i][1][n]
            n += 1
        # continue only if the correlation fits the chosen labeling scheme
        if aa_type =="W":
            if correlation_id =="CgCd2":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.25\n")
        elif aa_type== "H":
            if correlation_id == "CgCd2":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.25\n")
        elif aa_type=="V":
            if correlation_id == "CaCb":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("1.00\n")
        elif aa_type == "L" :
            if correlation_id == "CbCg":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("1.00\n")
        elif aa_type == "I":
            if correlation_id == "CaCb":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.33\n")
            elif correlation_id== "CbCg1":
                new_mat.append(mat[i])
                new_mat[-1][-1] = new_mat[-1][-1][:-1]
                new_mat[-1].append("0.17\n")
    new_mat.append(["\n"])
    return new_mat

def custom_isotopic_adj (mat,custom_location_file):
    aa_option = ["A", "R", "D", "N", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    f = open(custom_location_file)
    list = f.readlines()
    user_mat = [["aa", "C", "Ca", "Cb", "Cg","Cg1","Cg2", "Cd","Cd1","Cd2", "Ce","Ce1","Ce2","Ce3", "Cz","Cz2","Cz3","Ch2", "N", "Nd1", "Nd2", "Ne", "Ne1", "Ne2", "Nz", "Nh1", "Nh2"]]
    for i in range(len(list)):
        list[i] = list[i].rstrip(" \n")
        list[i] = list[i].split(" ")
        list[i].insert(0,aa_option[i])
        user_mat.append(list[i])

    for i in range (1,len(user_mat)):
        for j in range(1, len(user_mat[0])):
          if user_mat[i][j]=="-1.0":
              user_mat[i][j]="0"

    new_mat = []
    new_mat.append(mat[0])
    new_mat[0][-1] = new_mat[0][-1][:-1]
    new_mat[0].append("ILP\n")

    for i in range(1, len(mat)):
        new_mat.append(mat[i])
        new_mat[i][-1] = new_mat[i][-1][:-1]
        aa_in_mat = aa_option.index(new_mat[i][1][0]) + 1
        carbon_1_in_mat = "C"
        n = new_mat[i][1].find("C") + 1
        while new_mat[i][1][n] != "C":
            carbon_1_in_mat += new_mat[i][1][n]
            n += 1
        carbon_2_in_mat = ""
        while new_mat[i][1][n] != "(":
            carbon_2_in_mat += new_mat[i][1][n]
            n += 1
        carbon_1_id = user_mat[0].index(carbon_1_in_mat)
        carbon_2_id = user_mat[0].index(carbon_2_in_mat)
        probability = round(
            float(user_mat[aa_in_mat][carbon_1_id]) * float(user_mat[aa_in_mat][carbon_2_id]),2)
        new_mat[i].append(str(probability))
        new_mat[i][-1] = new_mat[i][-1] + "\n"
    new_mat.append(["\n"])
    return new_mat

# function that deletes improbable (with ILP lower than the criterion) assignments
def criterion_for_ILP (mat,min_Probability):
    new_mat=[]
    new_mat.append(mat[0])
    for i in range(1,len(mat)-1):
        if float(mat[i][-1]) >= min_Probability:
            new_mat.append(mat[i])
    new_mat.append(["\n"])
    return new_mat


# function that changes the overall value according to the new ambiguity
def overall_ambiguity (mat):
    lst_of_peak=[]
    for i in range(1, len(mat)-1):
        lst_of_peak.append(mat[i][0])
    for j in range(1, len(mat)-1):
        mat[j][4]=f"{lst_of_peak.count(mat[j][0])}I0S0M0L0X"
        mat[j][3]="I"
    return mat



# function that arranges the matches according to criteria: protein sequence or peak value

def arrange_by_seq (mat):
    new_mat=[]
    new_mat.append(mat[0])
    d={}
###create a dictionary
    for i in range(1,len(mat)-1):
        location=""
        j=1
        while mat[i][1][j]!="C":
            location+=mat[i][1][j]
            j+=1
        relative_carbons= get_first_carbon_ID(mat[i]) + get_sec_carbon_ID(mat[i])
        location_and_carbon= (int(location),relative_carbons)
        d[location_and_carbon]=mat[i]

    lst_of_key=sorted(d.keys())
    for i in range(len(lst_of_key)):
        new_mat.append(d.get(lst_of_key[i]))
    new_mat.append(["\n"])
    return new_mat


def arrange_by_peak (mat):
    new_mat=[]
    new_mat.append(mat[0])
    d={}
###create a dictionary
    for i in range(1,len(mat)-1):
        separate_index = mat[i][0].index("-")
        aa_location=""
        n = 1
        while mat[i][1][n] != "C":
            aa_location += mat[i][1][n]
            n += 1
        carbons_corelation= get_first_carbon_ID(mat[i])+get_sec_carbon_ID(mat[i])
        key=(float(mat[i][0][:separate_index]),float(mat[i][0][separate_index+1:]),int(aa_location),carbons_corelation)
        d[key]=mat[i]

    lst_of_key=sorted(d.keys())
    for i in range(len(lst_of_key)):
        new_mat.append(d.get(lst_of_key[i]))
    new_mat.append(["\n"])
    # cluster the ambiguity assignment together
    for k in range (1,len(new_mat)-1):
        if new_mat[k][0]==new_mat[k+1][0]:
            new_mat[k][-1]=new_mat[k][-1][:-1]
    return new_mat


# function that generalizes all the above functions into one and creates an INADEQUATE match file.
# sequence file is necessary because the function needs the length of the protein in order to convert the stored data into the real representation of it.
def write_adjust_matches_matrix(matche_file,seq_file,isotopic_tagging,custom_labeling_file=None):
    f=open(seq_file,"r")
    txt=f.read()
    step0 = 0
    for k in txt:
        if k.isupper():
            step0+=1
    step1=matches_as_matrix(matche_file)
    step2=same_aa_correlation_mat(step1,step0)
    step3 = rename_mtches_and_short(step2)
    if isotopic_tagging=="uniform":
        step4=uniform_adj(step3)
    elif isotopic_tagging=="1,3-Glycerol":
        step4= isotopic_1_3_Glycerol_adj(step3)
    elif isotopic_tagging=="2-Glycerol":
        step4= isotopic_2_Glycerol_adj(step3)
    elif isotopic_tagging== "Custom":
        step4= custom_isotopic_adj(step3,custom_labeling_file)
    step5= criterion_for_ILP(step4,min_Probability=0.20)
    step6= overall_ambiguity(step5)
    step7=arrange_by_peak(step6)
    str=""
    for i in range(len(step7)-1):
        for j in range(len(step7[0])-1):
            str+=step6[i][j]+"\t"
        str+=step6[i][-1]
        str+="\n"
    f=open(f"{matche_file[:-4]}_correlation_naming.list","w")
    f.write(str)
    return f"{matche_file[:-4]}_correlation_naming.list"


######################################################################################################################################################
######################################################################################################################################################


# STAGE 3- convert the match file into the cross peak analysis format

def match_list_into_peak_list (match_file):
    # convert the match file into a matrix
    f=open(match_file)
    lists_of_txt=f.readlines()
    f.close()
    mat=[]
    for i in range (len(lists_of_txt)):
        if lists_of_txt[i]!="\n":
            mat.append(lists_of_txt[i].strip("\n").split("\t"))

    str = "CC\nAssignment\tw1\tw2\n\n"

    for i in range(1,len (mat)):
        firs_aa_len = get_first_aa_len(mat[i])
        carbon1 = get_first_carbon_ID(mat[i])
        carbon2 = get_sec_carbon_ID(mat[i])
        # finding the peak value
        separd_index_on_real_peak = mat[i][0].index("-")
        carbon1_shift=mat[i][0][:separd_index_on_real_peak]
        carbon2_shift=mat[i][0][separd_index_on_real_peak+1:]

        str+=f"{mat[i][1][:1+firs_aa_len]}:{carbon1}-{carbon2}\t{carbon1_shift}\t{carbon2_shift}\n"

    new_file=open(match_file[:-5]+" for cross peak.list","w")
    new_file.write(str)
    return match_file[:-5]+" for cross peak.list"


######################################################################################################################################################
######################################################################################################################################################

# STAGE 4- getting the assignment file from BMRB website or convert an existing one from BMRB format to RAVEN format

# function that gets BMRB URL and writes the assignment file according to it
def get_BMRB_DATA (entry_number,save_in_location):
    url = f"https://bmrb.io/data_library/summary/index.php?bmrbId={entry_number}"

    # Send a GET request to the URL
    response = requests.get(url)

    # Parse the HTML content
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all the links on the page
    links = soup.find_all('a', href=True)

    # Loop through the links to find the one containing 'assigned_chemical_shifts'
    for link in links:
        if 'assigned_chemical_shifts' in link['href']:
            # Get the URL of the link
            url = urljoin(url, link['href'])

    website = urllib.request.urlopen(url)
    data=website.read().decode('utf-8')
    data=data.splitlines()
    for l in range(len(data)):
        data[l]=data[l].strip()

    desired_prefix = "_Assigned_chem_shift_list.Entry_ID"
    index = None

    for i, my_str in enumerate(data):
        if my_str.startswith(desired_prefix):
            index = i
            break
    name_in_line=data[index]
    while "  " in name_in_line:
        name_in_line = name_in_line.replace("  ", " ")
    name_in_line = name_in_line.split(" ")
    protein_name = name_in_line[1]

    for l in range(len(data)):
        data[l]=data[l].strip(" ")
    i=data.index("_Atom_chem_shift.Assigned_chem_shift_list_ID")+2
    new_data=""
    while data[i]!="stop_":
        new_data+= data[i]+"\n"
        i+=1

    j=data.index("_Atom_chem_shift.ID")
    headline=[]
    while data[j] != '':
        headline.append(data[j])
        j+=1
    reduce_headline=""
    for n in range (len (headline)):
        reduce_headline+=headline[n][headline[n].index(".")+1:] +" "


    f=open(save_in_location+f"/{protein_name} BMRB.shifts","w")
    f.write(reduce_headline + "\n"+ new_data)
    return save_in_location+f"/{protein_name} BMRB.shifts"



# function that converts form BMRB format to RAVEN format

def BMRM_to_assignment (BMRB_file):
    f=open(BMRB_file,"r")
    list_of_txt=f.readlines()
    f.close()

    #convert the match file into a matrix
    for i in range(len(list_of_txt)):
        list_of_txt[i]=list_of_txt[i].replace("\t"," ")
        list_of_txt[i]=list_of_txt[i].strip(" ")
        while "  " in list_of_txt[i]:
            list_of_txt[i]=list_of_txt[i].replace("  "," ")
        list_of_txt[i]=list_of_txt[i].split(" ")

    headline = list_of_txt[0]
    list_of_txt = list_of_txt[1:]
    # find the number of amino acids in the protein (that is assigned)
    num_of_aa=0
    for i in range(len(list_of_txt)):
        if int(list_of_txt[i][headline.index("Auth_seq_ID")])>num_of_aa:
            num_of_aa=int(list_of_txt[i][headline.index("Auth_seq_ID")])

    # create a matrix skeleton in Raven format and fill it with the BMRB data


    new_list = [["\t","N","C","Ca","Cb","Cg","Cg1","Cg2","Cd","Cd1","Cd2","Ce","Ce1","Ce2","Ce3","Cz","Cz2","Cz3","Ch2","Nd1","Nd2","Ne","Ne1","Ne2","Nz","Nh1","Nh2"]]
    for i in range (num_of_aa):
        new_list.append(["-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-"])
    for i in range (len(list_of_txt)):
        new_list[int(list_of_txt[i][headline.index("Auth_seq_ID")])][0]=list_of_txt[i][headline.index("Auth_comp_ID")]+list_of_txt[i][headline.index("Auth_seq_ID")]
        if list_of_txt[i][headline.index("Atom_type")] in ["C","c","N"]:
            new_list[int(list_of_txt[i][headline.index("Auth_seq_ID")])][new_list[0].index(list_of_txt[i][headline.index("Auth_atom_ID")].title())]=list_of_txt[i][headline.index("Val")]

    # create a new assignmt file in Raven format
    str=""
    for i in range(len(new_list)):
        for j in range (len(new_list[0])):
            if i>=1 and j==0:
                str+=new_list[i][j]+"\t\t"
            else:
                str+=new_list[i][j]+"\t"
        str+="\n"

    f=open(BMRB_file[:BMRB_file.index(".")]+" RAVEN format"+BMRB_file[BMRB_file.index("."):],"w")
    f.write(str)
    return BMRB_file[:BMRB_file.index(".")]+" RAVEN format"+BMRB_file[BMRB_file.index("."):]


######################################################################################################################################################
######################################################################################################################################################


# STAGE 5- functions that help the user select files in a convenient way

def print_name (str):
    new_str=""
    i=-1
    while str[i]!= "/":
        new_str+=str[i]
        i-=1
    return  f"\033[94m{new_str[::-1]}\033[0m"


def choose_file():
    root = tkinter.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    return  filedialog.askopenfilename()

def choose_location():
    root = tkinter.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    return  filedialog.askdirectory(title="Select a Directory")



# STAGE 6- Compare between two crosspeak files: original peak list and assigned one. Return a cross peak file with unassigned peaks:
# step 1- function that converts the cross peak file into a matrix (without headline)
# step 2- function that gets two matrices (original and assigned) and return a matrix with peak that appears only in the original
# step 3- function that deals with rounding issue with 5 as the last digit
# step 4- function that gets a cross peak matrix and a path to save and write cross peak file in that location
# step 5- function that generalizes all the above functions into one.

# step 1- function that converts the cross peak file into a matrix (without headline)
def cross_peak_mat (cross_peak_file):
    f= open (cross_peak_file)
    lists = f.readlines()
    mat=[]
    for i in range (3,len(lists)):
        lists[i]=lists[i].replace("\t", " ")
        lists[i] = lists[i].rstrip(" \n")
        while "  " in lists[i]:
            lists[i] = lists[i].replace("  ", " ")
        lists[i]=lists[i].lstrip(" ")
        lists[i]=lists[i].split(" ")
        mat.append(lists[i])
    return mat

# step 2- function that gets two matrices (original and assigned) and returns a matrix with peak that appears only in the original
def mat_of_unassigned_peak (original_mat,assigned_mat):
    mat=[]
    d_original_round={}
    d_assigned = {}
    for i in range(len(original_mat)):
        d_original_round[("{:.2f}".format(round(float(original_mat[i][1]),2)),"{:.2f}".format(round(float(original_mat[i][2]),2)))]=original_mat[i]
    for j in range(len(assigned_mat)):
        d_assigned[("{:.2f}".format(round(float(assigned_mat[j][1]),2)),"{:.2f}".format(round(float(assigned_mat[j][2]),2)))]=assigned_mat[j]

    for k in d_original_round.keys():
        if k not in d_assigned.keys():
            mat.append(d_original_round[k])

    plus=[]
    minus=[]
    plus_minos=[]
    minus_plus=[]
    first_plus=[]
    first_minus=[]
    sec_plus=[]
    sec_minus=[]


    for i in range(len(assigned_mat)):
        plus.append(("{:.3f}".format(float(assigned_mat[i][1])+0.005),"{:.3f}".format(float(assigned_mat[i][2])+0.005)))
        minus.append (("{:.3f}".format(float(assigned_mat[i][1])-0.005),"{:.3f}".format(float(assigned_mat[i][2])-0.005)))
        plus_minos.append(("{:.3f}".format(float(assigned_mat[i][1])+0.005),"{:.3f}".format(float(assigned_mat[i][2])-0.005)))
        minus_plus.append(("{:.3f}".format(float(assigned_mat[i][1])-0.005),"{:.3f}".format(float(assigned_mat[i][2])+0.005)))
        first_plus.append(("{:.3f}".format(float(assigned_mat[i][1])+0.005),"{:.3f}".format(float(assigned_mat[i][2]))))
        first_minus.append(("{:.3f}".format(float(assigned_mat[i][1])-0.005),"{:.3f}".format(float(assigned_mat[i][2]))))
        sec_plus.append(("{:.3f}".format(float(assigned_mat[i][1])),"{:.3f}".format(float(assigned_mat[i][2])+0.005)))
        sec_minus.append(("{:.3f}".format(float(assigned_mat[i][1])),"{:.3f}".format(float(assigned_mat[i][2])-0.005)))


    mat_2=[]
    for k in range(len(mat)):
        tup=(mat[k][1],mat[k][2])
        if tup in plus:
            mat_2.append(mat[k])
        elif tup in minus:
            mat_2.append(mat[k])
        elif tup in plus_minos:
            mat_2.append(mat[k])
        elif tup in minus_plus:
            mat_2.append(mat[k])
        elif tup in first_plus:
            mat_2.append(mat[k])
        elif tup in first_minus:
            mat_2.append(mat[k])
        elif tup in sec_plus:
            mat_2.append(mat[k])
        elif tup in sec_minus:
            mat_2.append(mat[k])

    new_mat=[]
    for i in mat:
        if i not in mat_2:
            new_mat.append(i)

    return new_mat

# step 3- function that deals with rounding issue with 5 as the last digit
def round_iso (o_mat, a_mat):
    matches={}
    round_1_plus_2 = {}
    round_1_minus_2 = {}
    plus_1_round_2 = {}
    minus_1_round_2 = {}

    for j in range(len(a_mat)):
        matches[(a_mat[j][1],a_mat[j][2])]= a_mat[j][1:]

    for i in range(len(o_mat)):
        w1= o_mat[i][1]
        w2= o_mat[i][2]
        if w1[-1]=="5" or w2[-1]=="5":
            round_1_plus_2 [("{:.2f}".format(round(float(o_mat[i][1]),2)), "{:.2f}".format(float(o_mat[i][2]) + 0.005))]=o_mat[i]
            round_1_minus_2[("{:.2f}".format(round(float(o_mat[i][1]),2)), "{:.2f}".format(float(o_mat[i][2]) - 0.005))]=o_mat[i]
            plus_1_round_2 [("{:.2f}".format(float(o_mat[i][1])+0.005), "{:.2f}".format(round(float(o_mat[i][2]),2)))]=o_mat[i]
            minus_1_round_2[("{:.2f}".format(float(o_mat[i][1])-0.005), "{:.2f}".format(round(float(o_mat[i][2]),2)))]=o_mat[i]

    mat=[]
    for k in round_1_plus_2.keys():
        if k in matches.keys():
            mat.append(round_1_plus_2[k])
    for k in round_1_minus_2.keys():
        if k  in matches.keys():
            mat.append(round_1_minus_2[k])
    for k in plus_1_round_2.keys():
        if k  in matches.keys():
            mat.append(plus_1_round_2[k])
    for k in minus_1_round_2.keys():
        if k  in matches.keys():
            mat.append(minus_1_round_2[k])

    new_mat = []
    for i in o_mat:
        if i not in mat:
            new_mat.append(i)

    return  new_mat


# step 4- function that gets a cross peak matrix and a path to save and write a cross peak file in that location
def write_unassigned_file (mat,path,name):
    if mat==[]:
        print (f"The two lists are identical and contain the same information.")
        return None
    text_of_file="CC\nAssignment\tw1\tw2\n\n"
    for i in range(len(mat)):
        for j in range (len(mat[0])):
            text_of_file += mat[i][j]+"\t"
        text_of_file += "\n"
    f=open(f"{path}/{name}.list","w")
    f.write(text_of_file)
    return f"{name}.list"


# step 5- function that generalizes all the above function into one.
def compare_original_to_assigned (original,assigned,place_to_save):
    original_mat=cross_peak_mat(original)
    assigned_mat=cross_peak_mat(assigned)
    step_2= mat_of_unassigned_peak(original_mat,assigned_mat)
    step_3= round_iso(step_2, assigned_mat)

    original_file_name=""
    i = original.find(".") - 1
    while original[i] != "/":
        original_file_name += original[i]
        i -= 1
    original_file_name=original_file_name[::-1]

    step_4=write_unassigned_file(step_3,place_to_save,f"unassigned peak in {original_file_name}")
    if step_4!= None:
        print(f"I have created a peak list file with all the unassigned peaks named: \033[94m{step_4}\033[0m\n")
    return




# STAGE 7- function that gets a BMRB protein entry number and creates an INADEQUATE assignment file suitable for POKY.
# step 1- function that gets the BMRB protein entry number and returns a text of the assignment part and the saving name for the file
# step 2- function that converts the assignment text into an assignment table
# step 3- function that converts the full assignment table into a reduced table that POKY can read
# step 4- function that creates the INADEQUATE assignment table
# step 5- function that generalizes all the above functions


def assignment_text_from_BMRB (entry_number):
    url = f"https://bmrb.io/data_library/summary/index.php?bmrbId={entry_number}"

    # Send a GET request to the URL
    response = requests.get(url)

    # Parse the HTML content
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all the links on the page
    links = soup.find_all('a', href=True)

    # Loop through the links to find the one containing 'assigned_chemical_shifts'
    for link in links:
        if 'assigned_chemical_shifts' in link['href']:
            # Get the URL of the link
            url = urljoin(url, link['href'])

    website = urllib.request.urlopen(url)
    data=website.read().decode('utf-8')
    data=data.splitlines()

    for l in range(len(data)):
        data[l] = data[l].strip()
    desired_prefix = "_Assigned_chem_shift_list.Entry_ID"
    index = None
    for i, my_str in enumerate(data):
        if my_str.startswith(desired_prefix):
            index = i
            break
    name_in_line = data[index]
    while "  " in name_in_line:
        name_in_line = name_in_line.replace("  ", " ")
    name_in_line = name_in_line.split(" ")

    j = data.index("_Atom_chem_shift.ID")
    headline = []
    while data[j] != '':
        headline.append(data[j])
        j += 1
    for n in range(len(headline)):
        headline[n] = headline[n][headline[n].index(".") + 1:]

    j+=1
    new_data=""
    while data[j]!="stop_":
        new_data+= data[j]+"\n"
        j+=1

    return (new_data,headline)

def bmrb_assignment_as_mat(str,headline):
    str=str.split("\n")
    mat = [headline]
    for i in range(len(str)-1):
        while "  " in str[i]:
            str[i]=str[i].replace("  "," ")
        mat.append(str[i].split(" "))
    return mat

def relevant_part_in_BMRB_mat (mat):
    new_mat=[]
    for i in range(1,len(mat)):
        if mat[i]!= [""] and mat[i][mat[0].index("Atom_isotope_number")] == "13" :
            new_mat.append([mat[i][mat[0].index("Comp_ID")]+mat[i][mat[0].index("Seq_ID")],mat[i][mat[0].index("Atom_ID")].title(),mat[i][mat[0].index("Atom_isotope_number")]+mat[i][mat[0].index("Atom_type")],mat[i][mat[0].index("Val")]])
            new_mat[-1][0]=Amino_acid(new_mat[-1][0])
    name=mat[2][mat[0].index("Entry_ID")]
    return (new_mat,name)

def INADEQUATE_adj_for_POKY_file(mat):
    new_mat=mat[:]
    for first_c in range(len(mat)-1):
        aa_type=mat[first_c][0].AA_ID()
        carbon1=mat[first_c][1].title()
        for sec_c in range(first_c+1, len(mat)):
            carbon2=mat[sec_c][1].title()

            x=mat[first_c][0]
            z=mat[sec_c][0]
            if mat[first_c][0]==mat[sec_c][0]:

                if aa_type == "G":
                    adjacent_carbon_options = [["C","Ca"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type in ["A", "C", "S"]:
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type in ["F","Y"]:
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"],["Cg","Cd1"],["Cg","Cd2"],["Cd1","Ce1"],["Cd2","Ce2"],["Ce1","Cz"],["Ce2","Cz"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "W":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"],["Cg","Cd1"],["Cg","Cd2"],["Cd2","Ce2"],["Cd2","Ce3"],["Ce2","Cz2"],["Ce3","Cz3"],["Cz2","Ch2"],["Cz3","Ch2"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "H":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"],["Cg","Cd2"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "V":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg1"],["Cb","Cg2"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "L":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"],["Cg","Cd1"],["Cg","Cd2"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "I":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg1"],["Cb","Cg2"],["Cg1","Cd1"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "K":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"],["Cg","Cd"],["Cd","Ce"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])


                if aa_type in ["R","E","Q","P"]:
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"],["Cg","Cd"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type in ["D","N","M"]:
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])

                if aa_type == "T":
                    adjacent_carbon_options = [["C","Ca"],["Ca","Cb"],["Cb","Cg2"]]
                    for option in adjacent_carbon_options:
                        if carbon1 in option and carbon2 in option:
                            new_mat.append([mat[first_c][0],carbon1+"&"+carbon2,mat[first_c][2],str(round(float(mat[first_c][3])+float(mat[sec_c][3]),2))])


    return new_mat

def write_BMRB_assignment_for_POKY (mat,name,location):
    text="Group\tAtom\tNuc\tShift\n"
    for i in range(len(mat)):
        for j in range(len(mat[0])-1):
            text+=f"{mat[i][j]}\t"
        text+= mat[i][-1]+"\n"
    path= f"{location}/POKY {name} DQ resonances.list"
    f=open(path,"w")
    f.write(text)
    return  f"\033[94mPOKY {name} DQ resonances.list\033[0m"

def entry_number_to_POKY_assignment_table (entry_number, folder):
    step1= assignment_text_from_BMRB (entry_number)
    step2= bmrb_assignment_as_mat (step1[0],step1[1])
    step3= relevant_part_in_BMRB_mat (step2)
    step4= INADEQUATE_adj_for_POKY_file (step3[0])
    step5= write_BMRB_assignment_for_POKY(step4, f"{step3[1]} protein", folder)
    return step5



def SQ_to_DQ (SQ_file,DQ_saving,saving_name):
    f = open(SQ_file)
    lists = f.readlines()
    f.close()

    mat=[]
    for i in range(len(lists)):
        while "  " in lists[i]:
            lists[i] = lists[i].replace("  ", " ")

        line= lists[i].strip(" \n\t").split(" ")
        if len(line) == 1:
            line=line[0].split("\t")
        if line != [""] :
            if i != 0:
                line[0] = Amino_acid(line[0])
                mat.append(line)

    DQ_res=INADEQUATE_adj_for_POKY_file (mat)

    text="Group\tAtom\tNuc\tShift\tSDev\t#\n"
    for n in range(len(DQ_res)):
        for m in range(len(DQ_res[n])-1):
            text += f"{DQ_res[n][m]}\t"
        text += f"{DQ_res[n][-1]}\n"

    path = f"{DQ_saving}/{saving_name}.list"
    f = open(path, "w")
    f.write(text)

    return






print("\033[1mDouble Quantum Assignment Table Creator\033[0m")
print("This Python script's goal is to make the necessary adjustments required for the use of RAVEN or POKY to assign proteins in J-based double quantum NMR experiment.")
first_choice=input("\nPlease choose the letter that corresponds to the option you desire:\na. Create an assignment table and a sequence file for the use of RAVEN Cross peaks assignment tool suitable for J-based double quantum experiment\nb. Process the received matches from Raven Cross peaks assignment tool\nc. Convert the matches format to the Cross peaks analysis format\nd. Find unassigned peaks after the assignment process\ne. Create a J-based double quantum assignment table suitable for POKY\n\n")

if first_choice=="a":
    print("You have selected option a - create assignment and sequence files suitable for RAVEN in J-based double quantum NMR experiment.\nPlease select the \033[91msequence\033[0m file of the studied protein:")
    seq_location = choose_file()
    print(f"\nyou have selected: {print_name(seq_location)}")

    source_for_assignment=input("\nPlease choose the \033[91mdatasource\033[0m for the assignment table:\na. Assignment table according to the BMRB entry number of the protein\nb. Assignment table file in BMRB format\nc. Assignment table file in RAVEN format\n\n")

    if source_for_assignment=="a":
        print("Please select the \033[91mpath\033[0m to the location that you would like to save the assignment file in\n")
        assignment_path =choose_location()
        print(f"You have selected: \033[94m{assignment_path}\033[0m\n")
        entry_num=input("Please enter the BMRB \033[91mEntry Number\033[0m of the studied protein\n\n")
        BMRB_name=get_BMRB_DATA(entry_num, assignment_path)
        assignment_table=BMRM_to_assignment(BMRB_name)
        two_file=adjustment_assignment_table_and_seq(assignment_table, seq_location,assignment_path)
        print (f"I have created 4 files for your use:\n(1) Duplicated sequence of amino acids named: {print_name(two_file[1])}\n(2) Assigned chemical shifts file according to the BMRB format named: {print_name(BMRB_name)}\n(3) Assigned chemical shifts file in the RAVEN format named: {print_name(assignment_table)}\n(4) Assigned chemical shifts with the adjustments for J-based double quantum experiment named: {print_name(two_file[0])}\n\nIn order to use the data correctly please load files (1) and (4) and the adequate peaklist into RAVEN.\nPick only the 'all correlations' option in the Output Configuration & Run.\n\n")


    if source_for_assignment=="b":
        print("Please choose the \033[91mBMRB\033[0m format shifts file\n")
        BMRB_file=choose_file()
        assignment_table=BMRM_to_assignment(BMRB_file)
        i = -1
        location=BMRB_file[:]
        while location[i] != "/":
            location = location[:i]
        two_file=adjustment_assignment_table_and_seq(assignment_table, seq_location,location)
        print(f"I have created 3 files for your use:\n(1) Duplicated sequence of amino acids named: {print_name(two_file[1])}\n(2) Assigned chemical shifts file in RAVEN format named: {print_name(assignment_table)}\n(3) Assigned shifts with the adjustments for J-based double quantum experiment named: {print_name(two_file[0])}\n\nIn order to use the data correctly please load files (1) and (3) and the adequate peaklist into RAVEN.\nPick only the 'all correlations' option in the Output Configuration & Run.\n\n")


    if source_for_assignment=="c":
        print("Please choose the \033[91massignment\033[0m shifts file\n")
        assignment_table=choose_file()
        i = -1
        location=assignment_table[:]
        while location[i] != "/":
            location = location[:i]
        two_file=adjustment_assignment_table_and_seq(assignment_table, seq_location,location)
        print(f"I have created 2 files for your use:\n(1) Duplicated sequence of amino acids named: {print_name(two_file[1])}\n(2) Assigned shifts with the adjustments for J-based double quantum experiment named: {print_name(two_file[0])}\n\nIn order to use the data correctly please load the created files and the adequate peaklist into RAVEN.\nPick only the 'all correlations' option in the Output Configuration & Run.\n\n")


    continue_yes_no=input("Would you like to continue and make the adjustments for the received matches file? (enter yes or no)\n")

    if continue_yes_no=="yes":
        print("Please select the \033[91mmatches\033[0m file:\n")
        match_file=choose_file()
        print(f"You have selected: {print_name(match_file)}\n")
        isotopic_tagging = input("Please select the \033[91mIsotope Labeling Scheme\033[0m that fits the sample:\na. Uniform\nb. 1,3-Glycerol\nc. 2-Glycerol\nd. Custom\n\n")
        if isotopic_tagging == "a":
            isotopic_tagging = "uniform"
        elif isotopic_tagging == "b":
            isotopic_tagging = "1,3-Glycerol"
        elif isotopic_tagging == "c":
            isotopic_tagging = "2-Glycerol"
        elif isotopic_tagging == "d":
            isotopic_tagging = "Custom"
        print(f"You have selected: \033[94m{isotopic_tagging}\033[0m Labeling\n")
        custom_isotopic_file=None
        if isotopic_tagging == "Custom":
            print("Please select the \033[91mCustom Isotope Labeling\033[0m file")
            custom_isotopic_file = choose_file()
            print(f"You have selected: {print_name(custom_isotopic_file)}\n\n")
        adjusted_match=write_adjust_matches_matrix(match_file,seq_location,isotopic_tagging,custom_isotopic_file)
        print(f"I created a new matches file containing only the possible 13C-13C correlations according to the isotope labeling probability and with the correct representation: {print_name(adjusted_match)}")
        continue_yes_no = input("Would you like to continue and make the adjustments for the cross peaks analysis tool? (enter yes or no)\n")
        if continue_yes_no=="yes":
            croos_peak_file=match_list_into_peak_list(adjusted_match)
            print(f"\nI created a new file in the cross peaks analysis format named: {print_name(croos_peak_file)}\ngood luck")
        else:
            print ("ok\ngood luck")
    else:
        print("ok\ngood luck")


if first_choice=="b":
    print("You have selected option b - Convert a matches file into the corresponding J-based double quantum corelations it represent.\nPlease select the \033[91msequence\033[0m file of the studied protein:\n")
    seq_location = choose_file()
    print(f"\nYou have selected: {print_name(seq_location)}\n\nPlease select the desired \033[91mmatches\033[0m file:")
    match_file=choose_file()
    print(f"You have selected: {print_name(match_file)}\n")
    isotopic_tagging = input("Please select the \033[91mIsotope Labeling Scheme\033[0m that fits the sample:\na. Uniform\nb. 1,3-Glycerol\nc. 2-Glycerol\nd. Custom\n\n")
    if isotopic_tagging=="a":
        isotopic_tagging= "uniform"
    elif isotopic_tagging=="b":
        isotopic_tagging = "1,3-Glycerol"
    elif isotopic_tagging == "c":
        isotopic_tagging = "2-Glycerol"
    elif isotopic_tagging =="d":
        isotopic_tagging = "Custom"
    print (f"You have selected: \033[94m{isotopic_tagging}\033[0m Labeling\n")
    custom_isotopic_file = None
    if isotopic_tagging == "Custom":
        print("Please select the \033[91mCustom Isotope Labeling\033[0m file")
        custom_isotopic_file=choose_file()
        print(f"You have selected: {print_name(custom_isotopic_file)}\n\n")
    adjusted_match = write_adjust_matches_matrix(match_file,seq_location,isotopic_tagging,custom_isotopic_file)
    print(f"I created a new matches file containing only the possible 13C-13C correlations according to the isotope labeling probability and with the correct representation: {print_name(adjusted_match)}")
    continue_yes_no = input("Would you like to continue and make the adjustment for the cross peaks analysis tool? (enter yes or no)\n")
    if continue_yes_no == "yes":
        croos_peak_file=match_list_into_peak_list(adjusted_match)
        print(f"I created a new file in the cross peaks analysis format named: {print_name(croos_peak_file)}\n\ngood luck")
    else:
        print("ok\ngood luck")



if first_choice=="c":
    print("You have selected option c - having a matches file and convert it into the Cross peaks analysis format.\nPlease select the \033[91mmatches\033[0m file:\n\n")
    adjusted_match=choose_file()
    print(f"\nYou have selected: {print_name(adjusted_match)}\n\n")
    croos_peak_file=match_list_into_peak_list(adjusted_match)
    print(f"I created a new file in the cross peaks analysis format named: {print_name(croos_peak_file)}\n\ngood luck")


if first_choice=="d":
    print("You have selected option d - find unassigned peaks after the assignment process\nThis option compares the original peak list and the one received after the assignment process, both in analysis format, and creates a new file which contains only the unassigned peaks\n")
    print ("Please select the \033[91moriginal\033[0m peak list:\n")
    original_peak_list=choose_file()
    print(f"You have selected: {print_name(original_peak_list)}\n")
    print ("Please select the \033[91massigned \033[0m peak list:\n")
    associated_peak_list=choose_file()
    print(f"You have selected: {print_name(associated_peak_list)}\n")
    print("Please select the \033[91mfolder\033[0m you would like to save the unassigned peaks in:\n")
    comparison_location=choose_location()
    print(f"You have selected: \033[94m{comparison_location}\033[0m\n")
    compare_original_to_assigned(original_peak_list, associated_peak_list, comparison_location)
    print("good luck")




if first_choice=="e":
    print(f"You have selected option e - create a J-based double quantum assignment table that can be uploaded into POKY.\nPlease choose your reference single quantum \033[91mresonance data\033[0m:")
    POKY_choice = input("\na. Resonances from BMRB database\nb. Existing resonances file\n\n")
    if POKY_choice == "a":
        print ("Please choose the \033[91mfolder\033[0m you would like to save the file in:\n")
        folder = choose_location()
        print(f"You have selected: \033[94m{folder}\033[0m\n")
        entry_number = input(f"Please type the BMRB \033[91mEntry Number\033[0m of the studied protein\n")
        POKY_assignment_table = entry_number_to_POKY_assignment_table(entry_number, folder)
        print( f"\nI created a double quantum resonances file that can be uploaded into POKY named:\n{POKY_assignment_table}\ngood luck")

    if POKY_choice == "b":
        print("Please choose the \033[91mexisting resonances\033[0m file:\n")
        SQ_resonance=choose_file()
        print(f"You have selected: {print_name(SQ_resonance)}\n")

        name = ""
        i = -1
        while SQ_resonance[i] != "/":
            name += SQ_resonance[i]
            i -= 1
        name = "double quantum resonance " + name[::-1]
        DQ_saving = SQ_resonance[:i]

        SQ_to_DQ(SQ_resonance,DQ_saving,name)
        print( f"I created a double quantum resonances file that can be uploaded into POKY named:\n\033[94m {name}\033[0m\ngood luck")






