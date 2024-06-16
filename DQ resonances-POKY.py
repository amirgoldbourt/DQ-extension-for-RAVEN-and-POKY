import tkinter
from tkinter import filedialog
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import urllib.request



def choose_saving_location():
    # This is a function that open a dialog box that allows the user to select the saving location
    root = tkinter.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    return  filedialog.askdirectory(title="Select a Directory")

def choose_file():
    root = tkinter.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    return  filedialog.askopenfilename()

def get_data_sheet_URL (entry_number):
    # This is a function that get from the user the BMRB entry number of the study protein and return its datasheet URL
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
            link_url = urljoin(url, link['href'])
            return link_url

    # If no entry was found
    print("No entry number was found.")
    return None



def assignment_text_from_BMRB (URL):
    # This is a function that get the datasheet URL and return a string of the shifts and a list of the headline

    # get list of lines
    website = urllib.request.urlopen(URL)
    data=website.read().decode('utf-8')
    data=data.splitlines()

    # "clean" the lines
    for l in range(len(data)):
        data[l] = data[l].strip(" ")

    # list of headlines
    j = data.index("_Atom_chem_shift.ID")
    headline = []
    while data[j] != '':
        headline.append(data[j])
        j += 1
    for n in range(len(headline)):
        headline[n] = headline[n][headline[n].index(".") + 1:]

    # str of shifts data
    j+=1
    shifts= ""
    while data[j]!="stop_":
        shifts+= data[j] + "\n"
        j+=1

    return (shifts, headline)



def bmrb_assignment_as_mat(shifts,headline):
    # This is a function that get the string of the shifts and the list of headline and combin them into matrix
    shifts=shifts.split("\n")
    mat = [headline]
    for i in range(len(shifts) - 1):
        while "  " in shifts[i]:
            shifts[i]=shifts[i].replace("  ", " ")
        mat.append(shifts[i].split(" "))
    while mat[-1] == [""]:
        mat=mat[:-1]
    return mat



def relevant_part_in_BMRB_mat (mat):
    # This is a function that reduce the original shits matrix and leave only the relevant parts for POKY resonances list

    # Relevant parts:
    new_mat=[]
    for i in range(1,len(mat)):
        if mat[i][mat[0].index("Atom_isotope_number")] == "13" :
            new_mat.append([mat[i][mat[0].index("Comp_ID")]+mat[i][mat[0].index("Seq_ID")],mat[i][mat[0].index("Atom_ID")].title(),mat[i][mat[0].index("Atom_isotope_number")]+mat[i][mat[0].index("Atom_type")],mat[i][mat[0].index("Val")]])

    # Uniform marking of amino acids by a single letter
    amino_acid_marks=[["A", "Ala"], ["C", "Cys"], ["D", "Asp"], ["E", "Glu"], ["F", "Phe"], ["G", "Gly"], ["H", "His"], ["I", "Ile"], ["K", "Lys"], ["L", "Leu"], ["M", "Met"], ["N", "Asn"], ["P", "Pro"], ["Q", "Gln"], ["R", "Arg"], ["S", "Ser"], ["T", "Thr"], ["V", "Val"], ["W", "Trp"], ["Y", "Tyr"]]
    for j in range (len(new_mat)):
        aa=""
        k=0
        while new_mat[j][0][k].isalpha():
            aa+=new_mat[j][0][k]
            k+=1
        location=new_mat[j][0][k:]
        for option in amino_acid_marks:
            if aa.title() in option:
                new_mat[j][0]=option[0]+location

    return new_mat


def INADEQUATE_adj_for_POKY(mat):
    # This function add to the single quantum resonances the Double Quantum ones
    new_mat=mat[:]
    for first_c in range(len(mat)-1):
        carbon1=mat[first_c][1].title()
        for sec_c in range(first_c+1, len(mat)):
            carbon2=mat[sec_c][1].title()

            # Consider only carbon that belong to the same amino acid
            if mat[first_c][0]==mat[sec_c][0]:
                aa_type = mat[first_c][0][0]

                # Consider only carbon that bound together - all the possible options
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
    # This function write the resonances list (including the Double Quantum) into a list file saving in the chosen location
    text="Group\tAtom\tNuc\tShift\n"
    for i in range(len(mat)):
        for j in range(len(mat[0])-1):
            text+=f"{mat[i][j]}\t"
        text+= mat[i][-1]+"\n"
    path= f"{location}/{name} DQ resonances.list"
    f=open(path,"w")
    f.write(text)
    return  f"\033[94m{name} DQ resonances.list\033[0m"



def entry_number_to_POKY_assignment_table (entry_number, folder):
    # This function use all of the above functions and with the input of BMRB entry number and the saving location create a new Double Quantum resonances file for POKY
    step1=get_data_sheet_URL (entry_number)
    if step1 != None:
        step2= assignment_text_from_BMRB (step1)
        step3=bmrb_assignment_as_mat(step2[0],step2[1])
        step4=relevant_part_in_BMRB_mat(step3)
        step5=INADEQUATE_adj_for_POKY(step4)
        step6= write_BMRB_assignment_for_POKY(step5, entry_number, folder)
    return step6




def SQ_to_DQ (SQ_file,DQ_saving,saving_name):
    f = open(SQ_file)
    lists = f.readlines()
    f.close()

    mat = []
    for i in range(len(lists)):
        while "  " in lists[i]:
            lists[i] = lists[i].replace("  ", " ")

        line = lists[i].strip(" \n\t").split(" ")
        if len(line) == 1:
            line=line[0].split("\t")
        if line != [""]:
            if i != 0:
                mat.append(line)

    DQ_res = INADEQUATE_adj_for_POKY(mat)

    text = "Group\tAtom\tNuc\tShift\tSDev\t#\n"
    for n in range(len(DQ_res)):
        for m in range(len(DQ_res[n]) - 1):
            text += f"{DQ_res[n][m]}\t"
        text += f"{DQ_res[n][-1]}\n"

    path = f"{DQ_saving}/{saving_name}"
    f = open(path, "w")
    f.write(text)

    return

print("\033[1mPOKY Double Quantum Resonances Converter\033[0m\n")
print("This Python script will create a double quantum assignment table base on single quantum resonances data either from BMRB database or from an existing resonance file in POKY format")
print("The double quantum assignment table could then be uploaded into POKY which allows a convenient analysis of double quantum experiment like INADEQUATE\n")

print("Please choose your single quantum resonances data source:")
data_source=input("a. Resonances from BMRB database\nb. Resonances from a single quantum experiment in POKY format\n")

if data_source=="a":
    print(f"Please choose the folder you would like to save the file in:")
    folder = choose_saving_location()
    print(f"You have selected: \033[94m{folder}\033[0m\n")
    entry_number = input(f"Please type the \033[91e BMRB entry number\033[0m of the studied protein\n")
    name=entry_number_to_POKY_assignment_table(entry_number, folder)
    print (f"Double quantum resonances file named {name} has been created")

if data_source =="b":
    print(f"Please choose the single quantum resonances file:")
    SQ_file=choose_file()
    print(f"You have selected: \033[94m{SQ_file}\033[0m\n")
    name=""
    i=-1
    while SQ_file[i] != "/":
        name+=SQ_file[i]
        i-=1
    name= "Double quantum resonance " + name[::-1]
    DQ_saving=SQ_file[:i]
    SQ_to_DQ(SQ_file, DQ_saving, name)


