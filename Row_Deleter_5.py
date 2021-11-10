# ********************************************************
# *   ----------------------------------------------     *
# * Janice Mui                                        *
# * jmui@student.unimelb.edu.au                  *
# * Last modification :: 11/09/2020                      *
# *   ----------------------------------------------     *
# ********************************************************

#this script takes the output of the script step3_4.py or NZ_Duplicate_Remover.py(or any other .csv with the same format) as input
#this script will delete any row in the file that comes from the organism Agrobacterium tumefaciens (ie. any line with the value in column Organism = Agrobacterium tumefaciens)
#this is needed even though A. tumefaciens sequences were excluded from the original BLASTp search, ipg-final_new.py can still retrieve genome sequences from A. tumefaciens 



import os
import sys
import pandas as pd

def main(list_csv):
    input_list = pd.read_csv(list_csv, header=0)
    print("The total number of rows in the file is: ", len(input_list))

    df = input_list.loc[input_list['Organism'].str.contains('Agrobacterium tumefaciens')].index
    output = input_list.drop(df)

    print("The total number of rows in the output file is: ", len(output))

    output.to_csv('No_Atumefaciens_output.csv', index=False)
   
if __name__ == '__main__':
    input_list = sys.argv[1]
    fname = os.path.basename(input_list)
    main(input_list)
