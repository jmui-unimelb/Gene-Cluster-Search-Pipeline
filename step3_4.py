# ********************************************************
# *   ----------------------------------------------     *
# * Bruna Moreira da Silva                               *
# * bmoreiradasi@student.unimelb.edu.au                  *
# * Last modification :: 24/06/2020                      *
# *   ----------------------------------------------     *
# ********************************************************

## Please read carefully:

## This script uses the previously multiple outputs (.csv) of Script 2, one for each of the enzymes in the pathway.
## These all have the same format – three columns: 'Nucleotide Accession', 'Protein', 'Organism'.
## This script will combine the multiple input files into one file (.csv) by removing duplicates for 'nucleotide accession' column,
## and removing the column 'Protein'.
## This script will check whether each nucleotide accession no. was in each input file or not, and will create a new column for each inputfile (same name).
## The output will be named as: 'Total_acc_output.csv'

## This script was developed with Python version 3.8.2
## Please check the necessary Libraries before running this script. 
## To run this script, open the terminal ('cmd' in Windows), go to the directory where the script and the input files are located, 
## and then type: python3 step3_4.py inputfilename1.csv inputfilename2.csv ...
## Where: 'inputfilename.csv' is the name of the csv file (output file from script 2), so please consider the correct name to replace this example.

# --------------------------------------------------------------------------------------------------------------------------------

##Libraries to import:
import os
import sys
import pandas as pd


#Create the main function
def main():

    # Count and print the number of input files passed in cmd
    print("Number of files passed as input: ", (len(sys.argv) - 1))
    
    files = []
    fnames = {}
    
    for filename in sys.argv[1:]:

        inputfile = pd.read_csv(filename, header=0, usecols=[0, 2], index_col = None)
        ifile = pd.read_csv(filename, header=0, usecols=[0], index_col = None)
        col = ifile["Nucleotide Accession"]
        print("The total number of Accession number rows in the file: ", filename, "is: ", len(inputfile.index))
        files.append(inputfile)
        fnames.update({filename:col})

    fnames_frame = pd.DataFrame.from_dict(fnames, orient = 'columns')

    #Concatenate all the input files and remove the duplicates for 'Nucleotide Accession no.'
    total = pd.concat(files, axis = 0, ignore_index = True)
    total.drop_duplicates(subset = "Nucleotide Accession", keep = 'first', inplace = True)

    check = []


    #Check whether the 'nucleotide accession' is in which inputfile and insert a column with 1 (yes) or 0 (no):
    for cols in fnames_frame.columns:

        check = total['Nucleotide Accession'].isin(fnames_frame[cols]).astype(int)
        total[cols] = check


    #Sum the number of 1(yes) that each nucleotide accession has, row by row, and create a new column for it
    column_list = list(total)
    column_list.remove('Nucleotide Accession')
    column_list.remove('Organism')
    total["Sum"] = total[column_list].sum(axis=1)
    
    #Save the output file and print the number of rows within it.
    print("The total number of Accession number rows in the result output file is: ", len(total.index))   
    total.to_csv("Total_acc_output.csv", index=False)

   
if __name__ == '__main__':
    main()
   




