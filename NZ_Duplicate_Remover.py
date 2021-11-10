# ********************************************************
# *   ----------------------------------------------     *
# * Janice Mui                                        *
# * jmui@student.unimelb.edu.au                  *
# * Last modification :: 06/10/2020                      *
# *   ----------------------------------------------     *
# ********************************************************

#this script takes the output of the script step3_4.py or Row_Deleter_5.py(or any other .csv with the same format)
#some of the sequences with the prefix NZ_ are redundant - NZ_XXX.1 and XXX.1 are the same genome sequence
#this script will remove the redundant sequences with the prefix NZ (NZ_XXX.1 is removed, XXX.1 is retained)

import os
import sys
import pandas as pd

def main(list_csv):

    input_list = pd.read_csv(list_csv, header=0)
    print("The total number of rows in the file is: ", len(input_list))

    df1 = input_list.loc[input_list['Nucleotide Accession'].str.contains('NZ_')]
    df2 = pd.concat([input_list, df1]).drop_duplicates(keep=False)

    col1 = df1['Nucleotide Accession'].tolist()
    col2 = df2['Nucleotide Accession'].tolist()

    for x in col1:
        if any(substring in x for substring in col2) is False:
            col2.append(x)

    output = input_list[input_list['Nucleotide Accession'].isin(col2)]

    print("The total number of rows in the output file is: ", len(output))

    output.to_csv('No_Redundant_NZ.csv', index=False)
   
if __name__ == '__main__':
    input_list = sys.argv[1]
    fname = os.path.basename(input_list)
    main(input_list)
