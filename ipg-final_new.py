# ********************************************************
# *   ----------------------------------------------     *
# * Bruna Moreira da Silva                               *
# * bmoreiradasi@student.unimelb.edu.au                  *
# * Last modification :: 10/06/2020                      *
# *   ----------------------------------------------     *
# ********************************************************

## Please read carefully:
## This script uses the E-Utilities APIs in order to download a specific Identical Protein Groups from protein NCBI databse.
## Info: "For any series of more than 100 requests, do this at weekends or outside USA peak times."

## This script uses the previously output file (.csv) as input and only consider the data from the column "Subject_Num".
## This script was developed with Python version 3.8.2
## Please check the necessary Libraries before running this script. 
## To run this script, open the terminal ('cmd' in Windows), open the directory where the script is located and also the 
## input file, and then type: python3 ipg-final_new.py inputfilename.csv
## Where: 'inputfilename.csv' is the name of the csv input file, so please consider the correct name to replace this example.
## Remember that the input files should be the output from script number 1, or following the same template header.
## This script will return a new csv file from Identical Protein Group Database, with 3 columns: 'Nucleotide Accession', 'Protein'
## and 'Organism'. The new file name will be as: 'inputfilename_cds-acc_output.csv'



# --------------------------------------------------------------------------------------------------------------------------------

##Libraries to import:
import os
import sys
import pandas as pd
import requests



#Inform email, Api_Key value for NCBI
email = "email_here@example"
akey = "XXX"
dbase = 'protein'

#Create the main function
def main(protein_csv):

    #Read the input file, consider only the second column named as "Subject_Num", and convert into a list separated by comma:
    inputfile = pd.read_csv(protein_csv, header=0, usecols=[1])
    pserie = inputfile['Subject_Num']
    plist = pserie.to_list()
    id_list = ",".join(plist)
    print("The total number of proteins in the input file is: ", len(inputfile.index))
    print("Please wait, this should take a while...Consider running this script at the recommended times from NCBI's rule.")


    #Assemble the Efetch URL using HTTP requests
    retrival = "ipg"    #The retrival mode of protein databse in this case is the identical protein group format
    mode = "text"
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    info = {'db':dbase, 'id':id_list, 'rettype':retrival, 'retmode':mode}
    r = requests.post(url=base, data=info)
    r_read = str(r.text)
    
    #Save the downloaded file from NCBI IPG Databse originally:
    with open(fname[:-4] + "_ipgdata.txt", "w") as file:
        file.write(r_read)
        file.close()
   
    #Open the original .txt downloaded and convert to csv table with desired columns:
    handler = open(fname[:-4] + "_ipgdata.txt", "r")  
    output = pd.read_csv(handler, header=0, sep='\t', engine='python', usecols=[2, 6, 8])
    print("The following result lines have no 'Nucleotide Accession' information: \n", output[output.isnull().any(axis=1)])
    output2 = output.dropna(subset = ['Nucleotide Accession'])
    print("The total number of Accession numbers in the output file is: ", len(output2.index))          

    #Save the output file, using the same input file name, adding "_cds-acc_output.csv", in the same directory:        
    output2.to_csv(fname[:-4]+"_cds-acc_output.csv", index=False)
    

   
if __name__ == '__main__':
    protein_list = sys.argv[1]
    fname = os.path.basename(protein_list)
    main(protein_list)
   




