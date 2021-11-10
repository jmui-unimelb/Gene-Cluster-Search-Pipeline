# ********************************************************
# *   ----------------------------------------------     *
# * Bruna Moreira                                        *
# * bmoreiradasi@student.unimelb.edu.au                  *
# * Last modification :: 06/10/2020                      *
# *   ----------------------------------------------     *
# ********************************************************

#To run this script:
#(1) Open the terminal window and go to the directory where the script and the input(.csv) file is located. 
#(1.1) Now type: python3, followed by the script name, the input file name (.csv). 
#(1.2) An example to illustrate: python3 evalue_script.py hittable.csv 0.0002
#(5) An output file (.csv) will be generated, at the same directory, with the filtered results (values <= informed threshold)
#(6) It will be also informed at the terminal screen the Total amount of filtered results.

# --------------------------------------------------------------------------------------------------------------------------------

#This part of the script import the libraries that you have previously installed in your machine (step 3)
import os
import sys
import pandas as pd

def main(list_csv, evalue): 
    
    ## Provide the Hit Table List information and give a header name to the 12 columns.  
    input_list = pd.read_csv(list_csv, names=["Query_Num", "Subject_Num", "Identity(%)", "Align_Length", "Mismatches", "Gap_Opens", "Query_Start", "Query_End", "Subject_Start", "Subject_End", "E-Value", "Bit_Score", "Positives(%)"]) 
    evalue = float(ev)
    for index, row in input_list.iterrows(): output_list = input_list.loc[lambda input_list: input_list['E-Value'] <= evalue]
    print("\n There are ",len(input_list)," entries in the input file. The filtered output result is: ", len(output_list),".Please check the generated output file for more information.")
    output_list.to_csv(fname[:-4]+"_output.csv", index=False)   

if __name__ == '__main__':
    input_list = sys.argv[1]
    ev = sys.argv[2]
    fname = os.path.basename(input_list)
    main(input_list, ev)

    



