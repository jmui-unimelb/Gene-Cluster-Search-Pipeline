#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

##Imported modules
import os
from os import system
import sys
import time
import multiprocessing
from multiprocessing import Process, freeze_support
import random
import fileinput

global GUI
global OUTBOX
global FRAME
global CURRENTDIR
global MGBPATH
global APPDATA
global TEMP
global DBPATH


#Find path to mgb files if run from another directory
pathfolders = os.environ['PATH'].split(os.pathsep)
pathfolders.reverse()
pathfolders.append(os.getcwd())
pathfolders.reverse()
CURRENTDIR = os.getcwd()
MGBPATH = ""
for folder in pathfolders:
  try:
    if "read_input_gui.py" in os.listdir(folder) and "guilib.py" in os.listdir(folder) and "empty.xhtml" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and "mgb_gui.py" in os.listdir(folder):
      MGBPATH = folder
      break
  except:
    pass
try:
  if  MGBPATH == "" and os.sep in sys.argv[0] and "read_input_gui.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]) and "guilib.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
    MGBPATH = sys.argv[0].rpartition(os.sep)[0]
    os.chdir(MGBPATH)
except:
  pass
if MGBPATH == "":
  print "Error: Please add the MultiGeneBlast installation directory to your $PATH environment variable before running the executable from another folder."
  sys.exit(1)
#Find path to Application Data
if sys.platform == ('win32'):
  APPDATA = os.environ['ALLUSERSPROFILE'] + os.sep + 'Application Data'
elif sys.platform == ('darwin'):
  APPDATA = os.path.expanduser("~") + "/Library/Application Support"
else:
  try:
    if os.path.exists(os.getcwd() + os.sep + "multigeneblast_data"):
      APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
    else:
      os.mkdir(os.getcwd() + os.sep + "multigeneblast_data")
      APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
  except:
    try:
      if os.path.exists(os.environ['HOME'] + os.sep + "multigeneblast_data"):
        APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
      else:
        os.mkdir(os.environ['HOME'] + os.sep + "multigeneblast_data")
        APPDATA = os.environ['HOME'] + os.sep + "multigeneblast_data"
    except:
      print "No permission to write to installation folder. Please change user or save somewhere else."
      sys.exit()
if sys.platform == ('darwin') or sys.platform == ('win32'):
  try:
    os.mkdir(APPDATA + os.sep + 'MultiGeneBlast')
    APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
  except:
    if os.path.exists(APPDATA + os.sep + 'MultiGeneBlast'):
      APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
#Find path to temporary files
if sys.platform == ('win32'):
  TEMP = os.environ['TEMP']
elif sys.platform == ('darwin'):
  TEMP = os.environ['TMPDIR']
else:
  try:
    os.mkdir(os.environ['HOME'] + os.sep + ".mgbtemp")
    TEMP = os.environ['HOME'] + os.sep + ".mgbtemp"
  except:
    TEMP = APPDATA
#Set other environment variables
os.environ['EXEC'] = MGBPATH + os.sep + "exec"
os.environ['PATH'] = os.environ['EXEC'] + os.pathsep + os.environ['PATH']

from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *
from string import ascii_letters
import urllib2
from urllib2 import Request,urlopen,URLError,HTTPError
import httplib
from httplib import BadStatusLine,HTTPException
import urllib
import tarfile
import cPickle as pickle
from Tkinter import *
from tkMessageBox import askyesno, showerror
import shutil

class Options(dict):
    """Simple options access class, first step to use Optparse"""
    def __init__(self, indict=None):
        if indict is None:
            indict = {}
        dict.__init__(self, indict)
        self.__initialized = True

    def __getattr__(self, attr):
        try:
            return self.__getitem__(attr)
        except KeyError:
            raise AttributeError(attr)

    def __setattr__(self, attr, value):
        if not self.__dict__.has_key('_Options__initialized'):
            return dict.__setattr__(self, attr, value)
        elif attr in self:
            dict.__setattr__(self, attr, value)
        else:
            self.__setitem__(attr, value)

##Functions necessary for this script
def get_sequence(fasta):
    """get the description and trimmed dna sequence"""
    in_file = open(fasta, 'r')
    content = in_file.readlines()
    in_file.close()
    content2 = []
    for i in content:
      if i != "":
        content2.append(i)
    content = content2
    while content[0] == "" or content[0] == "\n":
      content = content[1:]
    header = content[0]
    content = content[1:]
    content = [x.rstrip() for x in content]
    seq = "".join(content)
    if ">" not in header or ">" in seq:
      print >> sys.stderr, "FASTA file not properly formatted; should be single sequence starting with '>' and sequence name."
      sys.exit(1)
    return seq

def complement(seq):  
  complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}  
  complseq = []
  for base in seq:
    if base in complement.keys():
      complbase = complement[str(base)]
      complseq.append(complbase)
    else:
      complbase = 'n'
      complseq.append(complbase)
  return complseq 

def reverse_complement(seq):  
  seq = list(seq)  
  seq.reverse()  
  revcompl = complement(seq)
  revcomplstr = str()
  for i in revcompl:
    revcomplstr = revcomplstr + str(i)
  return  revcomplstr

def fastaseqlengths(proteins):
  names = proteins[0]
  seqs = proteins[1]
  seqlengths = {}
  a = 0
  for i in names:
    seq = seqs[a]
    seqlength = len(seq)
    seqlengths[i] = seqlength
    a += 1
  return seqlengths
  
def parsegenes(genes):
  genedict = {}
  genelist = []
  joinlist = []
  joindict = {}
  accessiondict = {}
  locustagdict = {}
  genenr = 0
  for i in genes:
    i = i.split("     gene            ")[0]
    join = "no"
    genenr += 1
    #Find gene location info for each gene
    if "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("   /")[0]
      while ")" not in location.replace(" ","")[-3:]:
        locationlist = location.split("\n")
        locationlist = locationlist[:-1]
        location = ""
        for i in locationlist:
          location = location + "i"
      location = location.replace("\n","")
      location = location.replace(" ","")
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("/")[0]
      while ")" not in location.replace(" ","")[-3:]:
        locationlist = location.split("\n")
        locationlist = locationlist[:-1]
        location = ""
        for i in locationlist:
          location = location + "i"
      location = location.replace("\n","")
      location = location.replace(" ","")
    else:
      location = i.split("\n")[0]
    #location info found in embl file, now extract start and end positions
    if "complement" in location.lower():
      location = location.lower()
      location = location.split("complement(")[1][:-1]
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.split("join(")[1][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","")
        end = location.split("..")[1]
        end = end.replace(">","")
      strand = "-"
    else:
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.split("join(")[1][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","")
        end = location.split("..")[1]
        end = end.replace(">","")
      strand = "+"
    if int(start) > int(end):
      start2 = end
      end2 = start
      start = start2
      end = end2
    #Correct for alternative codon start positions
    if "codon_start=" in i.lower():
      codonstart = i.lower().split("codon_start=")[1][0]
      if strand == "+":
        start = str(int(start) +  (int(codonstart) - 1))
      elif strand == "-":
        end = str(int(end) - (int(codonstart) - 1))
    #Find gene name for each gene, preferably locus_tag, than gene, than protein_ID
    a = 0
    b = 0
    genename = ""
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      if "protein_id=" in line:
        genename = (line.split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genename = genename.replace('"','')
        b += 1
      elif "protein_id=" in line.lower():
        genename = (line.lower().split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genename = genename.replace('"','')
        b += 1
      elif a == (nrlines - 1):
        genename = ""
        b += 1
      else:
        a += 1
    if len(genename) > 1:
      accnr = genename
    else:
      accnr = "no_accession_number_found"
    #Find gene name or locus tag
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      locustag = ""
      if "locus_tag=" in line:
        locustag = (line.split("locus_tag=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        locustag = locustag.replace('"','')
        b += 1
      elif "locus_tag=" in line.lower():
        locustag = (line.lower().split("locus_tag=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        locustag = locustag.replace('"','')
        b += 1
      elif a == (nrlines - 1):
        if locustag == "":
          locustag = "none"
        b += 1
      else:
        a += 1
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      if "gene=" in line:
        genename = (line.split("gene=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genename = genename.replace('"','')
        b += 1
      elif "gene=" in line.lower():
        genename = (line.lower().split("gene=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genename = genename.replace('"','')
        b += 1
      elif a == (nrlines - 1):
        if genename == "":
          genename = "none"
        b += 1
      else:
        a += 1
    if locustag != "none":
      locustagdict[accnr.rpartition(".")[0]] = locustag
    if accnr == "no_accession_number_found" and locustag != "none":
      accnr = locustag
      genename = locustag
    #Find sequence for each gene
    a = 0                                             ###Not all gbks contain protein sequences as translations, therefore sequences from gene clusters are now extracted from the database at a later stage if sequence is not in gbk
    b = 0
    sequence = ""
    while b < 2:
      line = i.split("\n")[a]
      if "translation=" in line:
        sequence = line.split("translation=")[1][1:]
        b += 1
        a += 1
        if line.count('"') > 1:
          sequence = line.split("translation=")[1][1:-1]
          b = 2
      elif "translation=" in line.lower():
        sequence = line.lower().split("translation=")[1][1:]
        b += 1
        a += 1
        if line.count('"') > 1:
          sequence = line.lower().split("translation=")[1][1:-1]
          b = 2
      elif a == (nrlines - 2) or a == (nrlines - 1):
        sequence = ""
        b = 2
      elif b == 1:
        if '"' in line:
          seqline = line.replace(" ","")
          seqline = seqline.split('"')[0]
          sequence = sequence + seqline
          b += 1
        else:
          seqline = line.replace(" ","")
          sequence = sequence + seqline
        a += 1
      else:
        a += 1
    sequence = sequence.upper()
    #Quality-check sequence
    forbiddencharacters = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for z in forbiddencharacters:
      if z in sequence:
        sequence = ""
    #Find annotation for each gene
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      if "product=" in line:
        annotation = line.split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif "product=" in line.lower():
        annotation = line.lower().split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif a == (nrlines - 1):
        annotation = "not_annotated"
        b += 1
      else:
        a += 1
    accessiondict[genename] = accnr
    if join == "yes":
      joinlist.append(genename)
      joindict[genename] = joinedparts2
    #Remove illegal chars
    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|} '''
    genename = "".join([char for char in genename if char not in illegal_chars])
    if len(genename) < 2:
      genename = "orf" + "_" + str(genenr)
    #Save data to dictionary
    if len(genename) > 1:
      genedict[genename] = [start,end,strand,annotation,sequence,accnr,genename]
    genelist.append(genename)
  return [genelist, genedict, joinlist, joindict, accessiondict, locustagdict]

def cleandnaseq(dnaseq):
  dnaseq = dnaseq.replace(" ","")
  dnaseq = dnaseq.replace("\t","")
  dnaseq = dnaseq.replace("\n","")
  dnaseq = dnaseq.replace("0","")
  dnaseq = dnaseq.replace("1","")
  dnaseq = dnaseq.replace("2","")
  dnaseq = dnaseq.replace("3","")
  dnaseq = dnaseq.replace("4","")
  dnaseq = dnaseq.replace("5","")
  dnaseq = dnaseq.replace("6","")
  dnaseq = dnaseq.replace("7","")
  dnaseq = dnaseq.replace("8","")
  dnaseq = dnaseq.replace("9","")
  dnaseq = dnaseq.replace("/","")
  return dnaseq

def extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict, locustagdict):
  names = []
  seqs = []
  for i in genelist:
    genename = i
    if locustagdict.has_key(genename):
      locustag = locustagdict[genename]
    elif locustagdict.has_key(genename.partition(".")[0]):
      locustag = locustagdict[genename.partition(".")[0]]
    elif accessiondict.has_key(genename.partition(".")[0]) and locustagdict.has_key(accessiondict[genename].partition(".")[0]):
      locustag = locustagdict[accessiondict[genename].partition(".")[0]]
    elif accessiondict.has_key(genename) and locustagdict.has_key(accessiondict[genename]):
      locustag = locustagdict[accessiondict[genename]]
    else:
      locustag = "no_locus_tag"
    #If suitable translation found in gbk, use that
    if len(genedict[i][4]) > 5:
      protseq = genedict[i][4]
      i = genedict[i]
    #If no suitable translation found in gbk, extract from DNA sequence
    else:
      i = genedict[i]
      y = int(i[0])
      z = int(i[1])
      if i[2] == "+":
        if genename in joinlist:
          geneseq = ""
          for j in joindict[genename]:
            partstart = int(j.split("..")[0])
            if ".." in j:
              partend = int(j.split("..")[1])
            else:
              partend = int(j)
            geneseqpart = dnaseq[(partstart - 1):partend]
            geneseq = geneseq + geneseqpart
        else:
          geneseq = dnaseq[(y - 1):z]
        protseq = translate(geneseq)
      elif i[2] == "-":
        if genename in joinlist:
          geneseq = ""
          joinlistrev = joindict[genename]
          joinlistrev.reverse()
          for j in joinlistrev:
            partstart = int(j.split("..")[0])
            if ".." in j:
              partend = int(j.split("..")[1])
            else:
              partend = int(j)
            geneseqpart = rc_dnaseq[(len(rc_dnaseq) - partend):(len(rc_dnaseq) - partstart + 1)]
            geneseq = geneseq + geneseqpart
        else:
          geneseq = rc_dnaseq[(len(rc_dnaseq) - z):(len(rc_dnaseq) - y + 1)]
        protseq = translate(geneseq)
    genedict[genename] = i[:-1] + [locustag]
    name = "input" + "|" + "c1" + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename + "|" + i[3] + "|" + i[5] + "|" + locustag
    seqs.append(protseq)
    names.append(name)
  proteins = [names,seqs,genelist,genedict,accessiondict]
  return proteins

def gbk2proteins(gbkfile):
  try:
    file = open(gbkfile,"r")
  except:
    print "Error: no or invalid input file: " + gbkfile
    sys.exit(1)
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  if "     CDS             " not in filetext or "\nORIGIN" not in filetext:
    print >> sys.stderr, "Exit: GBK file not properly formatted, no sequence found"
    sys.exit(1)
  cdspart = filetext.split("\nORIGIN")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = filetext.split("\nORIGIN")[1]
  dnaseq = cleandnaseq(dnaseq)
  dnaseqlength = len(dnaseq)
  rc_dnaseq = reverse_complement(dnaseq)
  #Extract genes
  genes = cdspart.split("     CDS             ")
  genes = genes[1:]
  genesdetails = parsegenes(genes)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  locustagdict = genesdetails[5]
  #Locate all genes on DNA sequence and translate to protein sequence
  proteins = extractprotfasta(genelist, genedict, dnaseq, rc_dnaseq, joinlist, joindict, accessiondict, locustagdict)
  textlines = filetext.split("\n//")[0]
  textlines = textlines.split("\n")
  accession = ""
  definition = ""
  definitionfound = "n"
  for i in textlines:
    if accession == "":
      if "LOCUS       " in i:
        j = i.split("LOCUS       ")[1]
        accession = j.split(" ")[0]
        if len(accession) < 4:
          accession = ""
    if definition == "":
      if "DEFINITION  " in i:
        j = i.split("DEFINITION  ")[1]
        definition = j
        definitionfound = "y"
    if definitionfound == "y":
      if "            " in i:
        definitionfound = "n"
        definition = definition + i.split("           ")[1]
      else:
        definitionfound = "n"
  #Test if accession number is probably real GenBank/RefSeq acc nr
  if testaccession(accession) == "n":
    accession = ""
  return [proteins, accession, dnaseqlength, definition]

def parse_dna_from_embl(embl_string):
    "Parse DNA sequence from EMBL input"
    seq_array = []
    lines = embl_string.split('\n')
    for line in lines:
        if line.lower().find('sequence') > -1:
            continue
        line = line.strip()
        line = line.rstrip('0123456789')
        line = line.rstrip('/')
        line = line.strip()
        seq_array.append(line)

    return "".join(seq_array)

def embl2proteins(emblfile):
  try:
    file = open(emblfile,"r")
  except:
    print "Error: no or invalid input file: " + emblfile
    sys.exit(1)
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  if "FT   CDS " not in filetext or ("\nSQ" not in filetext):
      log("Exit: EMBL file not properly formatted, no sequence found or no " \
          "CDS annotation found.", exit=True)
  cdspart = filetext.split("\nSQ  ")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = parse_dna_from_embl(filetext.split("\nSQ  ")[1])
  dnaseq = cleandnaseq(dnaseq)
  sequence = dnaseq
  if (sequence.count('N') + sequence.count('n') + sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
      log("Protein GBK/EMBL file provided. Please provide nucleotide " \
          "GBK/EMBL file.", exit=True)
  dnaseqlength = len(dnaseq)
  rc_dnaseq = reverse_complement(dnaseq)
  if dnaseqlength < 1:
      log("No sequence found in GBK/EMBL file. Please provide an annotated " \
          "nucleotide GBK/EMBL file with a DNA sequence.", exit=True)
  #Extract genes
  genes = cdspart.split("FT   CDS             ")
  genes = genes[1:]
  genesdetails = parsegenes(genes)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  locustagdict = genesdetails[5]
  #Locate all genes on DNA sequence and translate to protein sequence
  proteins = extractprotfasta(genelist, genedict, dnaseq, rc_dnaseq, joinlist, joindict, accessiondict, locustagdict)
  textlines = filetext.split("SQ   ")[0]
  textlines = textlines.split("\n")
  accession = ""
  definition = ""
  for i in textlines:
    if accession == "":
      if "AC   " in i:
        j = i.split("AC   ")[1]
        j = j.replace(" ","")
        accession = j.split(";")[0]
        if len(accession) < 4:
          accession = ""
      if definition == "":
        if "DE   " in i:
          j = i.split("DE   ")[1]
          definition = j
  #Test if accession number is probably real GenBank/RefSeq acc nr
  if testaccession(accession) == "n":
    accession = ""
  return [proteins, accession, dnaseqlength, definition]

def translate(sequence):
  #Translation table standard genetic code; according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
  transldict = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 
                 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', 
                 'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 
                 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 
                 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 
                 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 
                 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 
                 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
                 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 
                 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
                 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
                 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 
                 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 
                 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
                 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
                 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                 'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                 'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
                 'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
                 'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
                 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
                 'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
                 'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
                 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
                 'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
                 'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
                 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
                 'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
                 'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
                 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
                 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
  triplets = []
  triplet = ""
  a = 0
  for i in sequence:
    if a < 2:
      a += 1
      triplet = triplet + i
    elif a == 2:
      triplet = triplet + i
      triplets.append(triplet)
      triplet = ""
      a = 0
  protseq = ""
  aanr = 0
  for i in triplets:
    aanr += 1
    if aanr == 1:
      protseq = protseq + "M"
    else:
      if "n" in i or "N" in i or i not in transldict.keys():
        protseq = protseq + "X"
      else:
        protseq = protseq + transldict[i]
  if  len(protseq) > 0 and protseq[-1] == "*":
    protseq = protseq[:-1]
  return protseq

def writefasta(names,seqs,file):
  e = 0
  f = len(names) - 1
  try:
    out_file = open(file,"w")
    while e <= f:
      out_file.write(">")
      out_file.write(names[e])
      out_file.write("\n")
      out_file.write(seqs[e])
      out_file.write("\n")
      e += 1
    out_file.close()
  except(IOError,OSError,NotImplementedError):
    print >> sys.stderr, "FASTA file not created."

def testaccession(accession):
  #Test if accession number is probably real GenBank/RefSeq acc nr
  numbers = range(0,10)
  letters = []
  for i in ascii_letters:
    letters.append(i)
  nrnumbers = 0
  nrletters = 0
  for i in accession:
    if i in letters:
      nrletters += 1
    try:
      j = int(i)
      if j in numbers:
        nrnumbers += 1
    except:
      pass
  test = "y"
  if nrnumbers < 3 or nrletters < 1:
    test = "n"
  return test

def sortdictkeysbyvalues(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    return [key for value, key in items]

def sortdictvaluesbykeys(dict):
    items = [(key, value) for key, value in dict.items()]
    items.sort()
    return [value for key, value in items]

def sortdictkeysbyvaluesrev(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]
    
def sortdictkeysbyvaluesrevv(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    items.reverse()
    return [value for value, key in items]

def blastparse(blasttext, minseqcoverage, minpercidentity, seqlengths, seqdict, dbname, dbtype):
  blastdict = {}
  querylist = []
  blastlines = blasttext.split("\n")[:-1]
  #Filter for best blast hits (of one query on each subject)
  query_subject_combinations = []
  blastlines2 = []
  for i in blastlines:
    tabs = i.split("\t")
    query = tabs[0]
    subject = tabs[1]
    query_subject_combination = query + "_" + subject
    if query_subject_combination in query_subject_combinations:
      pass
    else:
      query_subject_combinations.append(query_subject_combination)
      blastlines2.append(i)
  blastlines = blastlines2
  frame_update()
  #Filters blastlines to get rid of hits that do not meet criteria
  blastlines2 = []
  for i in blastlines:
    tabs = i.split("\t")
    query = tabs[0]
    subject = tabs[1]
    perc_ident = int(tabs[2].split(".")[0])
    alignmentlength = float(tabs[3])
    evalue = str(tabs[10])
    blastscore = int(tabs[11].split(".")[0])
    if seqlengths.has_key(query):
      perc_coverage = (float(tabs[3]) / seqlengths[query]) * 100
    else:
      perc_coverage = 0
      print seqlengths
      print "Error: no sequence length found for", query
      sys.exit()
    if perc_ident > minpercidentity and (perc_coverage > minseqcoverage or alignmentlength > 40):
      blastlines2.append(i)
  blastlines = blastlines2
  #Goes through the blastlines. For each query, creates a querydict and hitlist, and adds these to the blastdict when finding the next query
  firstquery = "y"
  hitnr = 1
  for i in blastlines:
    frame_update()
    tabs = i.split("\t")
    query = tabs[0]
    if dbtype == "prot":
        subject = tabs[1].split("|")[3].split(".")[0]
    else:
        subject = tabs[1] + "_" + str(hitnr)
    internalblast = "n"
    if subject == "+" or subject == "-":
      internalblast = "y"
      subject = tabs[1].split("|")[4].split(".")[0]
    perc_ident = int(tabs[2].split(".")[0])
    alignmentlength = float(tabs[3])
    hit_start = str(tabs[8])
    hit_end = str(tabs[9])
    evalue = str(tabs[10])
    blastscore = int(tabs[11].split(".")[0])
    if seqlengths.has_key(query):
      perc_coverage = (float(tabs[3]) / seqlengths[query]) * 100
    else:
      seqlength = len(seqdict[query.split("|")[4]])
      perc_coverage = (float(tabs[3]) / seqlength) * 100
    if firstquery == "y": #Only until the first blastline with good hit
      if dbtype == "nucl" or testaccession(subject) == "y" or internalblast == "y":
        firstquery = "n"
        querylist.append(query)
        subjectlist = []
        querydict = {}
        subjectlist.append(subject)
        querydict[subject] = [perc_ident,blastscore,perc_coverage,evalue,hit_start,hit_end]
        last_query = query
    elif i == blastlines[-1]: #Only for the last blastline
      if query not in querylist:
        if dbtype == "nucl" or testaccession(subject) == "y" or internalblast == "y":
          blastdict[last_query] = [subjectlist,querydict]
          querylist.append(query)
          subjectlist = []
          querydict = {}
          subjectlist.append(subject)
          querydict[subject] = [perc_ident,blastscore,perc_coverage,evalue,hit_start,hit_end]
          blastdict[query] = [subjectlist,querydict]
          querylist.append(query)
      else:
        if dbtype == "nucl" or testaccession(subject) == "y" or internalblast == "y":
          subjectlist.append(subject)
          querydict[subject] = [perc_ident,blastscore,perc_coverage,evalue,hit_start,hit_end]
          blastdict[query] = [subjectlist,querydict]
    else: #For all but the first and last blastlines
      if query not in querylist:
        if dbtype == "nucl" or testaccession(subject) == "y" or internalblast == "y" or "genbank" not in dbname:
          blastdict[last_query] = [subjectlist,querydict]
          querylist.append(query)
          subjectlist = []
          querydict = {}
          subjectlist.append(subject)
          querydict[subject] = [perc_ident,blastscore,perc_coverage,evalue,hit_start,hit_end]
          last_query = query
      else:
        if dbtype == "nucl" or testaccession(subject) == "y" or internalblast == "y" or "genbank" not in dbname:
          subjectlist.append(subject)
          querydict[subject] = [perc_ident,blastscore,perc_coverage,evalue,hit_start,hit_end]
    hitnr += 1
  return [blastdict,querylist]

def generate_rgbscheme(nr):
  usablenumbers = [1,2,4,8,12,18,24,32,48,64,10000]
  lengthsdict = {1:[1,1,1],2:[1,1,2],4:[1,2,2],8:[2,2,2],12:[2,2,3],18:[2,3,3],24:[3,3,3],32:[3,3,4],48:[3,4,4],64:[4,4,4]}
  shortestdistance = 10000
  for i in usablenumbers:
    distance = i - nr
    if distance >= 0:
      if distance < shortestdistance:
        shortestdistance = distance
        closestnr = i
  toohigh = "n"
  if closestnr == 10000:
    toohigh = "y"
    closestnr = 64
  xyznumbers = lengthsdict[closestnr]
  x = xyznumbers[0]
  y = xyznumbers[1]
  z = xyznumbers[2]
  xpoints = []
  xpoint = (255/z)/2
  for i in range(x):
    xpoints.append(xpoint)
    xpoint += (255/x)
  ypoints = []
  ypoint = (255/z)/2
  for i in range(y):
    ypoints.append(ypoint)
    ypoint += (255/y)
  zpoints = []
  zpoint = (255/z)/2
  for i in range(z):
    zpoints.append(zpoint)
    zpoint += (255/z)
  colorlist = []
  for i in xpoints:
    for j in ypoints:
      for k in zpoints:
        rgb = "rgb(" + str(i) + "," + str(j) + "," + str(k) + ")"
        colorlist.append(rgb)
  if toohigh == "y":
    colorlist = colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist
  if closestnr == 24:
    colorlist = colorlist[:15] + colorlist[18:]
  if closestnr == 32:
    colorlist = colorlist[:21] + colorlist[24:]
  colorlist2 = []
  if closestnr == 1:
    colorlist2.append("red")
  if closestnr == 2:
    colorlist2.append("red")
    colorlist2.append("green")
  if closestnr == 4:
    colorlist2.append("red")
    colorlist2.append("green")
    colorlist2.append("blue")
    colorlist2.append("yellow")
  if closestnr == 8:
    neworder=[4,1,2,5,6,7,3,0]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 12:
    neworder=[6,3,5,9,7,2,11,4,8,1,10,0]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 18:
    neworder=[9,6,2,14,15,8,12,10,3,5,7,11,4,1,16,13,0]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 24:
    neworder=[15,12,9,6,5,0,21,1,16,14,8,17,2,23,22,3,13,7,10,4,18,20,19,11]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr == 32:
    neworder = [21,19,27,6,8,1,14,7,20,13,9,30,4,23,18,12,5,29,24,17,11,31,2,28,22,15,26,3,20,16,10,25]
    colorlist2 = [colorlist[i] for i in neworder]
  if closestnr > 32:
    random.shuffle(colorlist)
    colorlist2 = colorlist
  colorlist = colorlist2
  return colorlist

def _gene_arrow(start,end,strand,color,base,height):
    halfheight = height/2
    if start > end:
      start2 = end
      end2 = start
      start = start2
      end = end2
    oh = ShapeBuilder()
    if (end - start) < halfheight:
        if (strand == "+"):
            pointsAsTuples=[(start,base),
                            (end,base - halfheight),
                            (start,base - height),
                            (start,base)
    ]
        if (strand == "-"):
            pointsAsTuples=[(start,base - halfheight),
                            (end,base - height),
                            (end,base),
                            (start,base - halfheight)
                            ]
    else:
        if (strand == "+"):
            arrowstart = end-halfheight
            pointsAsTuples=[(start,base),
                            (arrowstart,base),
                            (end,base-halfheight),
                            (arrowstart,base - height),
                            (start,base - height),
                            (start,base)
                            ]
        if (strand == "-"):
            arrowstart = start + halfheight
            pointsAsTuples=[(start,base - halfheight),
                            (arrowstart,base - height),
                            (end,base - height),
                            (end,base),
                            (arrowstart,base),
                            (start,base - halfheight)
                            ]
    pg=oh.createPolygon(points=oh.convertTupleArrayToPoints(pointsAsTuples),strokewidth=1, stroke='black', fill=color)
    return pg        

def relativepositions(starts, ends, largestclustersize, screenwidth):
  rel_starts = []
  rel_ends = []
  #Assign relative start and end sites for visualization
  lowest_start = int(starts[0])
  leftboundary = lowest_start
  for i in starts:
    i = float(float(int(i) - int(leftboundary)) / largestclustersize) * float(screenwidth * 0.75)
    i = int(i)
    rel_starts.append(i)
  for i in ends:
    i = float(float(int(i) - int(leftboundary)) / largestclustersize) * float(screenwidth * 0.75)
    i = int(i)
    rel_ends.append(i)
  return [rel_starts,rel_ends]

def startendsitescheck(starts,ends):
  #Check whether start sites are always lower than end sites, reverse if necessary
  starts2 = []
  ends2 = []
  a = 0
  for i in starts:
    if int(i) > int(ends[a]):
      starts2.append(ends[a])
      ends2.append(i)
    else:
      starts2.append(i)
      ends2.append(ends[a])
    a += 1
  ends = ends2
  starts = starts2
  return [starts,ends]

def calculate_colorgroups(queryclusternumber,hitclusternumbers,queryclusterdata,internalhomologygroupsdict):
  #Extract data and generate color scheme
  hitclusterdata = queryclusterdata[queryclusternumber][1]
  queryclustergenes = hitclusterdata[hitclusterdata.keys()[0]][3]
  colorgroupsdict = {}
  colorgroupslengthlist = []
  colorgroupslist = []
  for hitclusternumber in hitclusternumbers:
    colorgroups = hitclusterdata[hitclusternumber][0][hitclusternumber]
    colorgroupsdict[hitclusternumber] = colorgroups
    colorgroupslengthlist.append(len(colorgroups))
    colorgroupslist.append(colorgroups)
  metacolorgroups = []
  internalgroups = internalhomologygroupsdict[queryclusternumber]
  for i in internalgroups:
    metagroup = []
    for j in i:
      for m in colorgroupslist:
        for l in m:
          if j in l:
            for k in l:
              if k not in metagroup:
                metagroup.append(k)
    if len(metagroup) > 1 and metagroup not in metacolorgroups:
      metacolorgroups.append(metagroup)
  #Generate RGB scheme
  rgbcolorscheme = generate_rgbscheme(len(metacolorgroups))
  rgbcolorscheme.append("#FFFFFF")
  #Create colorschemedict in which all genes that are hits of the same query gene get the same color
  colorschemedict = {}
  z = 0
  for i in queryclustergenes:
    for j in metacolorgroups:
      if i in j:
        for l in j:
          if colorschemedict.has_key(l):
            pass
          else:
            colorschemedict[l] = z
    if z in colorschemedict.values():
      z += 1
  return colorschemedict,rgbcolorscheme

def clusterblastresults(queryclusternumber,hitclusternumbers,queryclusterdata,colorschemedict,rgbcolorscheme, screenwidth, arch_search, allhits="n"):
  #print "Generating svg for cluster",queryclusternumber
  #Extract data and generate color scheme
  nrhitclusters = queryclusterdata[1][0]
  hitclusterdata = queryclusterdata[1][1]
  if nrhitclusters == 0:
    s = svg(x = 0, y = 0, width = (screenwidth * 0.75), height = (2770))
    viewbox = "0 0 " + str(screenwidth * 0.8) + " " + str(2950)
    s.set_viewBox(viewbox)
    s.set_preserveAspectRatio("none")
    return [s,[{},{},{}]]
  queryclustergenes = hitclusterdata[hitclusterdata.keys()[0]][3]
  queryclustergenesdetails = hitclusterdata[hitclusterdata.keys()[0]][4]
  colorgroupsdict = {}
  colorgroupslengthlist = []
  colorgroupslist = []
  for hitclusternumber in hitclusternumbers:
    colorgroups = hitclusterdata[hitclusternumber][0][hitclusternumber]
    colorgroupsdict[hitclusternumber] = colorgroups
    colorgroupslengthlist.append(len(colorgroups))
    colorgroupslist.append(colorgroups)
  #Find out whether hit gene cluster needs to be inverted compared to query gene cluster
  strandsbalancedict = {}
  for m in hitclusternumbers:
    hitclustergenesdetails = hitclusterdata[m][2]
    strandsbalance = 0
    for i in queryclustergenes:
      refstrand = queryclustergenesdetails[i][2]
      for j in colorgroupsdict[m]:
        if i in j:
          for k in j:
            if k in hitclusterdata[m][1] and hitclustergenesdetails[k][2] == refstrand:
              strandsbalance += 1
            elif k in hitclusterdata[m][1] and hitclusterdata[m][2][k][2] != refstrand:
              strandsbalance = strandsbalance - 1
    strandsbalancedict[m] = strandsbalance
  #Generate coordinates for SVG figure
  qnrgenes = len(queryclustergenes)
  qstarts =[]
  qends = []
  qstrands =[]
  qcolors = []
  for i in queryclustergenes:
    qgenedata = queryclustergenesdetails[i]
    if qgenedata[0] > qgenedata[1]:
      qstarts.append(qgenedata[0])
      qends.append(qgenedata[1])
    else:
      qstarts.append(qgenedata[1])
      qends.append(qgenedata[0])
    qstrands.append(qgenedata[2])
    if colorschemedict.has_key(i):
      qcolors.append(colorschemedict[i])
    else:
      qcolors.append("white")
  qstarts_ends = startendsitescheck(qstarts,qends)
  qstarts = qstarts_ends[0]
  qends = qstarts_ends[1]
  hdata = {}
  for m in hitclusternumbers:
    hitclustergenes = hitclusterdata[m][1]
    hitclustergenesdetails = hitclusterdata[m][2]
    hnrgenes = len(hitclustergenes)
    hstarts =[]
    hends = []
    hstrands =[]
    hcolors = []
    for i in hitclustergenes:
      hgenedata = hitclustergenesdetails[i]
      if int(hgenedata[0]) > int(hgenedata[1]):
        hstarts.append(hgenedata[0])
        hends.append(hgenedata[1])
      else:
        hstarts.append(hgenedata[1])
        hends.append(hgenedata[0])
      hstrands.append(hgenedata[2])
      if colorschemedict.has_key(i):
        hcolors.append(colorschemedict[i])
      else:
        hcolors.append("white")
    #Invert gene cluster if needed
    if strandsbalancedict[m] < 0:
      hstarts2 = []
      hends2 = []
      hstrands2 = []
      for i in hstarts:
        hstarts2.append(str(100000000 - int(i)))
      hstarts = hstarts2
      hstarts.reverse()
      for i in hends:
        hends2.append(str(100000000 - int(i)))
      hends = hends2
      hends.reverse()
      hstarts, hends = hends, hstarts
      for i in hstrands:
        if i == "+":
          hstrands2.append("-")
        elif i == "-":
          hstrands2.append("+")
      hstrands = hstrands2
      hstrands.reverse()
      hcolors.reverse()
    #Sort genes properly and remove duplicates
    stranddict = {}
    colorsdict = {}
    y = 0
    sortstarts = []
    for n in hstarts:
      while n in sortstarts:
        n = str(int(n) + 1)
      sortstarts.append(n)
    for n in sortstarts:
      stranddict[int(n)] = hstrands[y]
      w = y + 1
      try:
        nextstart = sortstarts[w]
      except:
        nextstart = 0
      color = hcolors[y]
      if color == "white":
        while int(nextstart) == int(n):
          if len(hcolors) > w and hcolors[w] != 'white':
            color = hcolors[w]
            break
          w += 1
          try:
            nextstart = sortstarts[w]
          except:
            break
      if not colorsdict.has_key(int(n)) or colorsdict[int(n)] == "white":
        colorsdict[int(n)] = color
      y += 1
    hstarts = [int(l) for l in hstarts]
    #hstarts = dict.fromkeys(hstarts).keys()
    hstarts.sort()
    hstarts = [str(n) for n in hstarts]
    hends = [int(l) for l in hends]
    #hends = dict.fromkeys(hends).keys()
    hends.sort()
    hends = [str(n) for n in hends]
    hstrands = sortdictvaluesbykeys(stranddict)
    hcolors = sortdictvaluesbykeys(colorsdict)
    hstarts_ends = startendsitescheck(hstarts,hends)
    hstarts = hstarts_ends[0]
    hends = hstarts_ends[1]
    hdata[m] = [hstarts,hends,hstrands,hcolors]
  #Resize all gene clusters to normal sizes
  #Find largest hit cluster
  approxclustersizes = []
  for m in hitclusternumbers:
    hstarts,hends,hstrands,hcolors = hdata[m]
    x = 0
    first = -1
    last = int(hends[-1])
    for n in hcolors:
      if n != "white" and first == -1:
        first = int(hstarts[x])
        last = int(hends[x])
      elif n != "white" and first != -1:
        last = int(hends[x])
      x += 1
    approxclustersizes.append((int(last)-int(first)))
  approxclustersizes.append(int(qends[-1]) - int(qstarts[0]))
  largestsize = int(int(max(approxclustersizes)) + 5000)
  #Resize all clusters
  hdata2 = {}
  savedpositions = {}
  for m in hitclusternumbers:
    hstarts,hends,hstrands,hcolors = hdata[m]
    x = 0
    first = -1
    last = 0
    for n in hcolors:
      if n != "white" and first == -1:
        first = min(int(hstarts[x]), int(hends[x]))
        if max(int(hstarts[x]), int(hends[x])) > last:
          last = max(int(hstarts[x]), int(hends[x]))
      elif n != "white" and first != -1:
        if min(int(hstarts[x]), int(hends[x])) < first:
          first = min(int(hstarts[x]), int(hends[x]))
        if max(int(hstarts[x]), int(hends[x])) > last:
          last = max(int(hstarts[x]), int(hends[x]))
      x += 1
    approxclustersize = (int(last)-int(first))
    piecetobeadded = (largestsize - approxclustersize) / 2
    #if min([(first - int(hstarts[0])),(int(hends[-1]) - last)]) < piecetobeadded - 1:
    #  piecetobeadded = min([(first - int(hstarts[0])),(int(hends[-1]) - last)])
    if piecetobeadded < 0:
      piecetobeadded = 0
    newcstart = int(first) - piecetobeadded
    newcend = int(last) + piecetobeadded
    firstentry = 1000000000
    lastentry = -1
    x = 0
    for i in hstarts:
      hstart = int(i)
      hend = int(hends[x])
      if firstentry == 1000000000 and hend >= newcstart:
        firstentry = x
        lastentry = x + 1
      elif hstart <= newcend:
        lastentry = x + 1
      x += 1
    #print str(cstart) + " " + str(cend) + " " + str(newcstart) + " " + str(newcend)
    hstarts = hstarts[firstentry:lastentry]
    hends = hends[firstentry:lastentry]
    hstrands = hstrands[firstentry:lastentry]
    hcolors = hcolors[firstentry:lastentry]
    hdata2[m] = [hstarts,hends,hstrands,hcolors]
    savedpositions[m] = [hstarts,hends]
  hdata = hdata2
  #Find cluster size of largest cluster of query & all hit clusters assessed
  clustersizes = []
  for m in hitclusternumbers:
    pstarts = [int(n) for n in hdata[m][1]]
    pends = [int(n) for n in hdata[m][0]]
    locations = pstarts + pends
    hclustersize = abs(max(locations) - min(locations))
    clustersizes.append(hclustersize)
  qpositions = [int(q) for q in qends] + [int(q) for q in qstarts]
  qclustersize = abs(max(qpositions) - min(qpositions))
  clustersizes.append(qclustersize)
  largestclustersize = max(clustersizes)
  #Find relative positions
  qrelpositions = relativepositions(qstarts,qends,largestclustersize, screenwidth)
  qrel_starts = qrelpositions[0]
  qrel_ends = qrelpositions[1]
  qdata = [qrel_starts,qrel_ends,qstrands,qcolors]
  hdata2 = {}
  qdata2 = []
  q_adjusted = False
  for m in hitclusternumbers:
    hclustersize = int(hdata[m][1][-1]) - int(hdata[m][0][0])
    hrelpositions = relativepositions(hdata[m][0],hdata[m][1],largestclustersize, screenwidth)
    hrel_starts = hrelpositions[0]
    hrel_ends = hrelpositions[1]
    #Center-align smallest gene cluster
    if largestclustersize == hclustersize:
      if q_adjusted == False:
        q_adjusted = True
        qrel_ends2 = []
        qrel_starts2 = []
        for i in qrel_starts:
          qrel_starts2.append(int(i) + int(float(float((largestclustersize - qclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
        for i in qrel_ends:
          qrel_ends2.append(int(i) + int(float(float((largestclustersize - qclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
        qrel_ends = qrel_ends2
        qrel_starts = qrel_starts2
    else:
      hrel_ends2 = []
      hrel_starts2 = []
      for i in hrel_starts:
        hrel_starts2.append(int(i) + int(float(float((largestclustersize - hclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
      for i in hrel_ends:
        hrel_ends2.append(int(i) + int(float(float((largestclustersize - hclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
      hrel_ends = hrel_ends2
      hrel_starts = hrel_starts2
    hdata2[m] = [hrel_starts,hrel_ends,hdata[m][2],hdata[m][3]]
    qdata2 = [qrel_starts,qrel_ends,qdata[2],qdata[3]]
  hdata = hdata2
  qdata = qdata2
  s = svg(x = 0, y = 0, width = (screenwidth * 0.75), height = 2770)
  viewbox = "0 0 " + str(screenwidth * 0.8) + " " + str(2680)
  s.set_viewBox(viewbox)
  s.set_preserveAspectRatio("none")
  #Add line behind query gene cluster gene arrows, except for architecture searches
  oh = ShapeBuilder()
  if arch_search == "n":
    group = g()
    group.addElement(oh.createLine(10,35,10 + (screenwidth * 0.75),35, strokewidth = 1, stroke = "grey"))
    s.addElement(group)
  #Add query gene cluster gene arrows
  a = 0
  y = 0
  for x in range(qnrgenes):
    group = g()
    #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
    if qcolors[a] == "white":
      group.addElement(_gene_arrow(10 + qrel_starts[a],10 + qrel_ends[a],qstrands[a],rgbcolorscheme[-1],40,10))
    else:
      group.addElement(_gene_arrow(10 + qrel_starts[a],10 + qrel_ends[a],qstrands[a],rgbcolorscheme[qcolors[a]],40,10))
    #Can be used for domains
    #group.addElement(oh.createRect(rel_starts[a],45,(rel_ends[a]-rel_starts[a]),10, strokewidth = 2, stroke = "black", fill="#237845"))
    if allhits == "n":
      group.set_id("q" + str(queryclusternumber) + "_" + str(hitclusternumbers[0]) + "_" + "%s"%x)
    else:
      group.set_id("all_" + str(queryclusternumber) + "_0_" + "%s"%x)
    s.addElement(group)
    if y == 0:
      y = 1
    elif y == 1:
      y = 0
    a += 1
  for m in hitclusternumbers:
    group = g()
    group.addElement(oh.createLine(10,35 + 50 * (hitclusternumbers.index(m) + 1),10 + (screenwidth * 0.75),35 + 50 * (hitclusternumbers.index(m) + 1), strokewidth = 1, stroke = "grey"))
    s.addElement(group)
    #Add hit gene cluster gene arrows
    hitclustergenes = hitclusterdata[m][1]
    hrel_starts = hdata[m][0]
    hnrgenes = len(hrel_starts)
    hrel_ends = hdata[m][1]
    hstrands = hdata[m][2]
    hcolors = hdata[m][3]
    a = 0
    y = 0
    for x in range(hnrgenes):
      group = g()       
      #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
      if hcolors[a] == "white":
        group.addElement(_gene_arrow(10 + hrel_starts[a],10 + hrel_ends[a],hstrands[a],rgbcolorscheme[-1],40 + 50 * (hitclusternumbers.index(m) + 1),10))
      else:
        group.addElement(_gene_arrow(10 + hrel_starts[a],10 + hrel_ends[a],hstrands[a],rgbcolorscheme[hcolors[a]],40 + 50 * (hitclusternumbers.index(m) + 1),10))
      #Can be used for domains
      #   group.addElement(oh.createRect(rel_starts[a],45,(rel_ends[a]-rel_starts[a]),10, strokewidth = 2, stroke = "black", fill="#237845"))
      if allhits == "n":
        group.set_id("h" + str(queryclusternumber) + "_" + str(m) + "_" + "%s"%x)
      else:
        group.set_id("all_" + str(queryclusternumber) + "_" + str(m) + "_" + "%s"%x)
      s.addElement(group)
      if y == 0:
        y = 1
      elif y == 1:
        y = 0
      a += 1
  return [s,[qdata,hdata,strandsbalancedict,savedpositions]]

def log(message, exit=False, retcode=1, stdout=False):
    if GUI == "y":
      OUTBOX.text_insert(message + "\n")
      FRAME.update()
      if exit:
        sys.exit(retcode)
    else:
      "Log to stderr and logfile and optionally exit"
      if stdout:
        print message
      else:
        print >> sys.stderr, message
      logfile = open('multigeneblast.log', 'a', 1)
      logfile.write(message + '\n')
      logfile.close()
      if exit:
        sys.exit(retcode)

def inputinstructions():
    return """MultiGeneBlast 1.1.0 arguments:

Usage: multigeneblast [options]
Options (x is an integer number)

-in <file name>       : Query file name: GBK/EMBL file for homology search,
                        FASTA file with multiple protein sequences for
                        architecture search
-from <x>             : Start position of query region
-to <x>               : End position of query region
-genes <acc,acc;...>  : Accession codes of genes constituting query
                        multigene module
-out <folder name>    : Output folder in which results will be stored
-db <db name>         : Blast database to be queried (default: genbank_mf)
-cores <x>            : Number of parallel CPUs to use for threading
                        (default: all)
-hitspergene          : Number of Blast hits per query gene to be taken
                        into account (default: 250).
-minseqcov <x>        : Minimal % coverage of a Blast hit on hit protein
                        sequence to be taken into account (default: 25)
-minpercid <x>        : Minimal % identity of a Blast hit on hit protein
                        sequence to be taken into account (default: 30)
-distancekb <x>       : Maximum kb distance between two blast hits to be
                        counted as belonging to the same locus (default: 10)
-syntenyweight <x>    : Weight of synteny conservation in hit sorting score
                      : (default: 0.5)
-muscle <y/n>         : generate Muscle multiple sequence alignments of
                        all hits of each input gene  (default: n)
-outpages <x>         : Maximum number of output pages (with 50 hits each)
                        to be generated  (default: 5)"""

def default_options(opts):
    #Implement defaults
    opts.db = "genbank_mf"
    opts.cores = "all"
    opts.minseqcov = 25
    opts.minpercid = 30
    opts.screenwidth = 1024
    opts.hitspergene = 250
    opts.distancekb = 20000
    opts.muscle = "n"
    opts.startpos = "N/A"
    opts.endpos = "N/A"
    opts.ingenes = "N/A"
    opts.pages = 5
    opts.gui = "n"
    opts.syntenyweight = 0.5

def invalidoptions(argument):
    print "From the command line, input multigeneblast -help for more information."
    if len(argument) > 0:
        log("Invalid options input: %s" % argument, exit=True)

def collect_identifiers(options):
    #identify option identifiers
    identifiers = []
    for i in options:
        if i[0] == "-":
            if i not in identifiers:
                identifiers.append(i)
            else:
                invalidoptions("No '-' in given options or option given twice.")
    return identifiers

def determine_cpu_nr(cores):
    #Determine number of CPUs used
    if cores == "all":
        try:
            nrcpus = multiprocessing.cpu_count()
        except(IOError,OSError,NotImplementedError):
            nrcpus = 1
    else:
        try:
            nrcpus = multiprocessing.cpu_count()
        except(IOError,OSError,NotImplementedError):
            nrcpus = 1
        if cores < nrcpus:
            nrcpus = cores
    return nrcpus

def process_identifiers(identifiers, opts, options):
    infile, startpos, endpos, ingenes, outfolder = "n","n","n","n","n"
    genes_tag_used = "n"
    fromto_tag_used = "n"
    fastafile = "n"
    global CURRENTDIR
    for i in identifiers:
        if i == "-help" or i == "--help" or i == "-h":
          print inputinstructions()
          sys.exit(1)
        else:
          value = options[options.index(i) + 1].strip()
          if i == "-from":
              fromto_tag_used = "y"
              if genes_tag_used == "y":
                print "Please either the -from and -to tags, or the -genes tag to select genes."
                invalidoptions(i)
              if value.isdigit():
                  opts.startpos = int(value)
                  startpos = "y"
              else:
                  invalidoptions(i)
          elif i == "-to":
              fromto_tag_used = "y"
              if genes_tag_used == "y":
                print "Please either the -from and -to tags, or the -genes tag to select genes."
                invalidoptions(i)
              if value.isdigit():
                  opts.endpos = int(value)
                  endpos = "y"
              else:
                  invalidoptions(i)
          elif i == "-genes":
              genes_tag_used = "y"
              if fromto_tag_used == "y":
                print "Please either the -from and -to tags, or the -genes tag to select genes."
                invalidoptions(i)
              if "," in value:
                  opts.ingenes = [gene for gene in value.split(",") if gene != ""]
                  ingenes = "y"
              else:
                  invalidoptions(i)
          elif i == "-in":
              if sys.platform == ('win32') and os.sep not in value:
                value = CURRENTDIR + os.sep + value
              elif os.sep not in value or value[0] != os.sep:
                value = CURRENTDIR + os.sep + value
              root, ext = os.path.splitext(value)
              if ext.lower() not in [".gbk",".gb",".genbank",".embl",".emb",".fasta",".fas",".fa",".fna"]:
                  print "Please supply input file with valid GBK / EMBL extension (homology search) or FASTA extension (architecture search)."
                  invalidoptions(i)
              if value in os.listdir(".") or (os.sep in value and os.path.exists(value.rpartition(os.sep)[0]) and value.rpartition(os.sep)[2] in os.listdir(value.rpartition(os.sep)[0])):
                  opts.infile = value
                  infile = "y"
                  if ext.lower() in [".fasta",".fas",".fa",".fna"]:
                    fastafile = "y"
              else:
                  print "Specified input file not found..."
                  invalidoptions(i)
          elif i == "-out":
	      value = value.replace(".","").replace("_","")
              if not value.replace("_","").replace("_","").replace("..","").replace(os.sep,"").isalnum():
                  print "Not a valid output folder name. Please use alpha-numerical characters only"
                  invalidoptions(i)
              if sys.platform == ('win32') and value.count(os.sep) == 1 and value[0] == os.sep:
                invalidoptions(i)
              opts.outputfolder = value
              if opts.outputfolder[0] == os.sep and ".." in opts.outputfolder:
                invalidoptions(i)
              elif os.sep in value[0] and not os.path.exists(value.rpartition(os.sep)[0]):
                invalidoptions(i)
              elif os.sep in opts.outputfolder and ".." in opts.outputfolder:
                startdir = CURRENTDIR
                while ".." in opts.outputfolder:
                  if ".." not in opts.outputfolder.partition(os.sep)[0]:
                    invalidoptions(i)
                  opts.outputfolder = opts.outputfolder.partition(os.sep)[2]
                  startdir = startdir.rpartition(os.sep)[0]
                  if len(startdir) < 1:
                    invalidoptions(i)
                if opts.outputfolder[0] == os.sep:
                  opts.outputfolder = startdir + opts.outputfolder
                else:
                  opts.outputfolder = startdir + os.sep + opts.outputfolder
              elif os.sep not in opts.outputfolder:
                opts.outputfolder = CURRENTDIR + os.sep + opts.outputfolder
              elif opts.outputfolder[0] == os.sep:
                opts.outputfolder = opts.outputfolder
              elif os.sep in opts.outputfolder and os.sep not in opts.outputfolder[0]:
                opts.outputfolder = CURRENTDIR + os.sep + opts.outputfolder
              else:
                invalidoptions(i)
              if os.path.exists(opts.outputfolder):
                  print "Warning: Overwriting existing folder"
                  for xhtmlfile in [filename for filename in os.listdir(opts.outputfolder) if "xhtml" in filename]:
                    os.remove(opts.outputfolder + os.sep + xhtmlfile)
                  outfolder = "y"
              else:
                  try:
                    os.mkdir(opts.outputfolder)
                  except:
                    invalidoptions(i)
                  outfolder = "y"
          elif i == "-db":
              global MGBPATH
              global DBPATH
              value = value.partition(".pal")[0].partition(".nal")[0]
              if sys.platform == ('win32') and os.sep not in value:
                value = CURRENTDIR + os.sep + value
              elif os.sep not in value or value[0] != os.sep:
                value = CURRENTDIR + os.sep + value
              if not value + ".pal" in os.listdir(MGBPATH) and not value + ".pal" in os.listdir(".") and not value + ".pal" in os.listdir(CURRENTDIR) and not (os.sep in value and os.path.exists(value.rpartition(os.sep)[0]) and value.rpartition(os.sep)[2] + ".pal" in os.listdir(value.rpartition(os.sep)[0])):
                if not value + ".phr" in os.listdir(MGBPATH) and not value + ".phr" in os.listdir(".") and not (os.sep in value and os.path.exists(value.rpartition(os.sep)[0]) and value.rpartition(os.sep)[2] + ".phr" in os.listdir(value.rpartition(os.sep)[0])):
                  if not value + ".nal" in os.listdir(MGBPATH) and not value + ".nal" in os.listdir(".") and not value + ".nal" in os.listdir(CURRENTDIR) and not (os.sep in value and os.path.exists(value.rpartition(os.sep)[0]) and value.rpartition(os.sep)[2] + ".nal" in os.listdir(value.rpartition(os.sep)[0])):
                    if not value + ".nhr" in os.listdir(MGBPATH) and not value + ".nhr" in os.listdir(".") and not (os.sep in value and os.path.exists(value.rpartition(os.sep)[0]) and value.rpartition(os.sep)[2] + ".nhr" in os.listdir(value.rpartition(os.sep)[0])):
                      print "Error: Database not found; database should have accompanying .phr, .psq and .pin files."
                      invalidoptions(i)
              opts.db = value
              if value + ".pal" in os.listdir(MGBPATH) or value + ".nal" in os.listdir(MGBPATH):
                DBPATH = MGBPATH
              elif value + ".pal" in os.listdir(CURRENTDIR) or value + ".nal" in os.listdir(CURRENTDIR):
                DBPATH = CURRENTDIR
              elif value + ".pal" in os.listdir(".") or value + ".nal" in os.listdir("."):
                DBPATH = os.getcwd()
              elif os.sep in value and value[0] != os.sep and os.path.exists(os.getcwd() + os.sep + value.rpartition(os.sep)[0]) and ((value.rpartition(os.sep)[2] + ".pal" in os.listdir(os.getcwd() + os.sep + value.rpartition(os.sep)[0])) or (value.rpartition(os.sep)[2] + ".nal" in os.listdir(os.getcwd() + os.sep + value.rpartition(os.sep)[0]))):
                DBPATH = os.getcwd() + os.sep + value.rpartition(os.sep)[0]
                opts.db = value.rpartition(os.sep)[2]
              elif os.sep in value and os.path.exists(value.rpartition(os.sep)[0]) and (value.rpartition(os.sep)[2] + ".pal" in os.listdir(value.rpartition(os.sep)[0]) or value.rpartition(os.sep)[2] + ".nal" in os.listdir(value.rpartition(os.sep)[0])):
                DBPATH = value.rpartition(os.sep)[0]
                opts.db = value.rpartition(os.sep)[2]
              else:
                print "Error: Database not found; database should have accompanying .phr, .psq and .pin or .nhr, .nsq and .nin files."
                invalidoptions(i)
              os.environ['BLASTDB'] = DBPATH
              if opts.db + ".pal" in os.listdir(DBPATH):
                opts.dbtype = "prot"
              elif opts.db + ".nal" in os.listdir(DBPATH):
                opts.dbtype = "nucl"
          elif i == "-cores":
              if value.isdigit() and int(value) in range(1,1000):
                  opts.cores = int(value)
              else:
                  invalidoptions(i)
          elif i == "-minseqcov":
              if value.isdigit() and int(value) in range(0,100):
                  opts.minseqcov = int(value)
              else:
                  invalidoptions(i)
          elif i == "-minpercid":
              if value.isdigit() and int(value) in range(0,100):
                  opts.minpercid = int(value)
              else:
                  invalidoptions(i)
          elif i == "-distancekb":
              if value.isdigit() and int(value) in range(1,100):
                  opts.distancekb = int(value) * 1000
              else:
                  print "Error: please select a number between 1-100."
                  invalidoptions(i)
          elif i == "-syntenyweight":
              if (value.isdigit() or (value.count(".") == 1 and value.partition(".")[0].isdigit() and value.partition(".")[2].isdigit())) and float(value) <= 2.0 and float(value) >= 0.0:
                  opts.syntenyweight = float(value)
              else:
                  print "Error: please select a number between 0.0 and 2.0."
                  invalidoptions(i)
          elif i == "-hitspergene":
              if value.isdigit() and int(value) in range(50,10001):
                  opts.hitspergene = int(value)
              else:
                  print "Error: please select a number between 50-10000."
                  invalidoptions(i)
          elif i == "-muscle":
              if value == "y" or value == "n":
                  opts.muscle = value
              else:
                  invalidoptions(i)
          elif i == "-outpages":
              if value.isdigit() and int(value) in range(0,41):
                  opts.pages = int(value)
              else:
                  print "Error: please select a number between 0-40."
                  invalidoptions(i)
          else:
              invalidoptions(i)
    #Stop process if options are provided with a FASTA file that do not belong with it
    if fastafile == "y" and (genes_tag_used == "y" or fromto_tag_used == "y"):
      print "Error: -from, -to and -genes tags are incompatible with architecture search (FASTA input)"
      sys.exit(1)
    #Stop process if inadequate options are supplied.
    if infile == "n" or ("n" in [startpos, endpos] and ingenes == "n" and ".fa" not in opts.infile.partition(".")[1] + opts.infile.partition(".")[2]) or outfolder == "n":
        print "Input error. An input file, an outputfolder and a query region (for EMBL/GBK inputs) must be supplied."
        invalidoptions(" ".join(options))

def parse_options(args, opts):
    #Run GUI if no arguments supplied
    if len(args) < 2:
      args = "-h"
    default_options(opts)
    #Read user-specified options which may override defaults
    if len(args) >= 2:
        options = args
        if "-" in options[-1] and (args[1] != "-help" and args[1] != "--help" and args[1] != "-h"):
            invalidoptions(options[-1])
        identifiers = collect_identifiers(options)
        process_identifiers(identifiers, opts, options)
    nrcpus = determine_cpu_nr(opts.cores)
    opts.nrcpus = nrcpus

def generate_architecture_data(fastafile):
  try:
    file = open(fastafile,"r")
  except:
    log("Error: no or invalid input file: " + fastafile, exit=True)
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  fasta_entries = [">" + entry.replace("\n\n","\n").replace("\n\n","\n") for entry in filetext.split(">")][1:]
  querytags = []
  seqs = []
  seqlengths = {}
  seqdict = {}
  for entry in fasta_entries:
    if "\n" not in entry or len(entry.partition("\n")[2]) < 2:
      log("FASTA file wrongly formatted at. Please check your input file.")
      log("Wrong entry: " + entry.replace("\n",""), exit=True)
    fname = entry.partition("\n")[0][1:]
    #Generate name without forbidden characters
    forbiddencharacters = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    fname_censored = ""
    for z in fname:
      if z not in forbiddencharacters:
        fname_censored = fname
    fname = fname_censored.replace(" ","_")[:20].replace('|','_')
    seq = entry.partition("\n")[2].replace(" ","").replace("\n","")
    if fname in querytags:
      log("Non-unique sequence name in input FASTA file. Please reformat and try again.", exit=True)
    querytags.append(fname)
    seqs.append(seq)
    seqlengths[fname] = len(seq)
    seqdict[fname] = seq
  #Determine alignment distribution / lengths for display in SVG to put in names as pseudo-locations
  names = []
  genedict = {}
  accessiondict = {}
  totalnt = 0
  for tag in querytags:
    startsite = str(totalnt) 
    fname = "input|c1|" + startsite
    totalnt += (seqlengths[tag] * 3)
    endsite = str(totalnt)
    fname = fname + "-" + endsite + "|+|" + tag + "|" + tag + "|" + tag + "|" + tag
    seqlengths[fname] = seqlengths[tag]
    totalnt += 100
    names.append(fname)
    accessiondict[tag] = tag
    genedict[tag] = [startsite, endsite, "+", tag, seqdict[tag], tag, tag]
  genelist = names
  proteins = [names, seqs, genelist, genedict, accessiondict]
  return proteins, querytags, seqdict, names, seqs

def parse_absolute_paths(infile):
    #Parse absolute paths if found
    originalfilename = infile
    if "/" in infile or "\\" in infile:
        lastpos = max([infile.rfind("\\"),infile.rfind("/")])
        originpath = infile[:(lastpos + 1)]
        infile = infile[(lastpos + 1):]
        #if os.getcwd() != originalfilename[:lastpos] and os.getcwd() != originalfilename[:lastpos].replace("/","\\"):
        #  shutil.copyfile(originalfilename, infile)
    return infile

def read_input_file(infile, startpos, endpos, ingenes, gui, outbox=None, frame=None):
  global GUI
  GUI = gui
  global OUTBOX
  OUTBOX = outbox
  global FRAME
  FRAME = frame
  log("Reading and parsing input GenBank file...")
  #infile = parse_absolute_paths(infile)
  ext = infile.rpartition(".")[2]
  if ext.lower() in ["gbk","gb","genbank"]:
    proteins = gbk2proteins(infile)
  elif ext.lower() in ["embl","emb"]:
    proteins = embl2proteins(infile)
  elif ext.lower() in ["fasta","fas","fa","fna"]:
    nucname = "Architecture Search FASTA input"
    genomic_accnr = ""
    dnaseqlength = 0
    proteins, querytags, seqdict, names, seqs = generate_architecture_data(infile)
    writefasta(names,seqs,"query.fasta")
    arch_search = "y"
    return proteins, genomic_accnr, dnaseqlength, nucname, querytags, names, seqs, seqdict, arch_search
  arch_search = "n"
  genomic_accnr = proteins[1]
  dnaseqlength = proteins[2]
  nucname = proteins[3]
  proteins = proteins[0]
  querytags = []
  z = 0
  names = []
  seqs = []
  seqdict = {}
  if startpos != "N/A" and endpos != "N/A":
    for i in proteins[0]:
      seq = proteins[1][z]
      pstartpos = int(i.split("|")[2].split("-")[0])
      pendpos = int(i.split("|")[2].split("-")[1])
      if (pstartpos > startpos and pstartpos < endpos) or (pendpos > startpos and pendpos < endpos):
        names.append(i)
        seqs.append(seq)
        seqdict[i] = seq
        querytags.append(i.split("|")[4])
      z += 1
    if len(names) == 0:
      log("Error: no genes found within the specified region of the input file.", exit=True)
  elif ingenes != "N/A":
    for i in proteins[0]:
      seq = proteins[1][z]
      if i.split("|")[4] in ingenes or i.split("|")[6] in ingenes or i.split("|")[6].partition(".")[0] in ingenes or i.split("|")[6].partition(".")[0] in [gene.partition(".")[0] for gene in ingenes] or i.split("|")[7] in ingenes:
        names.append(i)
        seqs.append(seq)
        seqdict[i] = seq
        querytags.append(i.split("|")[4])
      z += 1
    if len(names) == 0:
      log("Error: no genes found with these names in the input file.", exit=True)
  writefasta(names,seqs,"query.fasta")
  return proteins, genomic_accnr, dnaseqlength, nucname, querytags, names, seqs, seqdict, arch_search

def internal_blast(minseqcoverage, minpercidentity, names, proteins, seqdict, nrcpus):
  #Run BLAST on gene cluster proteins of each cluster on itself to find internal homologs, store groups of homologs - including singles - in a dictionary as a list of lists accordingly
  log("Finding internal homologs..")
  internalhomologygroupsdict = {}
  clusternumber = 1
  #Make Blast db for internal search
  makeblastdbcommand = "makeblastdb -in query.fasta -out query.fasta -dbtype prot"
  makeblastdb_stdout = os.popen4(makeblastdbcommand)
  makeblastdb_stdout = makeblastdb_stdout[1].read()
  z = 0
  while "error" in makeblastdb_stdout.lower():
    log(makeblastdb_stdout)
    log("Error running BLAST. Retrying...")
    makeblastdb_stdout = os.popen4(makeblastdbcommand)
    makeblastdb_stdout = makeblastdb_stdout[1].read()
    if z > 2:
      log("Error generating internal Blast database, exiting. Please check your system.", exit=True)
    z += 1
  #Run and parse BLAST search
  blastsearch = "blastp  -db query.fasta -query query.fasta -outfmt 6 -max_target_seqs 1000 -evalue 1e-05 -out internal_input.out -num_threads " + str(nrcpus)
  blast_stdout = os.popen4(blastsearch)
  blast_stdout = blast_stdout[1].read()
  z = 0
  while "error" in blast_stdout.lower():
    log(blast_stdout)
    log("Error running BLAST. Retrying...")
    blast_stdout = os.popen4(blastsearch)
    blast_stdout = blast_stdout[1].read()
    if z > 2:
      log("Error running Blast, exiting. Please check your system.", exit=True)
    z += 1
  blastoutput = open("internal_input.out","r").read()
  seqlengths = fastaseqlengths(proteins)
  iblastinfo = blastparse(blastoutput, minseqcoverage, minpercidentity, seqlengths, seqdict, "internal", "prot")
  iblastdict = iblastinfo[0]
  #find and store internal homologs
  groups = []
  for j in names:
    frame_update()
    if iblastdict.has_key(j):
      hits = iblastdict[j][0]
      group = []
      for k in hits:
        if k[:2] == "h_":
          group.append(k[2:])
        elif k.count("|") > 4:
          group.append(k.split("|")[4])
        else:
          group.append(k)
      if j.split("|")[4] not in group:
        group.append(j.split("|")[4])
      x = 0
      for l in groups:
        for m in group:
          if m in l:
            del groups[x]
            for n in l:
              if n not in group:
                group.append(n)
            break
        x += 1
      group.sort()
      groups.append(group)
    else:
      groups.append([j.split("|")[4]])
  internalhomologygroupsdict[clusternumber] = groups
  return internalhomologygroupsdict, seqlengths

def runblast(args, dbtype):
  #return ##CAN BE UNCOMMENTED TEMPORARILY FOR FAST TESTING
  if dbtype == "prot":
    blastsearch = "blastp  " + args
  else:
    blastsearch = "tblastn  " + args
  blast_stdout = os.popen4(blastsearch)
  blast_stdout = blast_stdout[1].read()
  z = 0
  while "error" in blast_stdout.lower():
    print blast_stdout
    print "Error running BLAST. Retrying..."
    blast_stdout = os.popen4(blastsearch)
    blast_stdout = blast_stdout[1].read()
    if z > 2:
      print "Error running Blast, exiting. Please check your system."
      sys.exit(1)
    z += 1

def db_blast(names, seqs, db, nrcpus, hitspergene, dbtype="prot"):
  ##Run BLAST on genbank_mf database
  log("Running NCBI BLAST+ searches on GenBank database..")
  queryclusternames = names
  queryclusterseqs = seqs
  writefasta(queryclusternames,queryclusterseqs,"input.fasta")
  #blastsearch = "blastp  -db " + db + " -query input.fasta -outfmt 6 -max_target_seqs 1000 -num_descriptions 1000 -num_alignments 500 -evalue 1e-05 -out input.out -num_threads " + str(nrcpus)
  #blast_stdout = os.popen4(blastsearch)  ##CAN BE COMMENTED OUT TEMPORARILY FOR FAST TESTING
  #blast_stdout = blast_stdout[1].read()
  #z = 0
  #while "error" in blast_stdout.lower():
  #  log(blast_stdout)
  #  log("Error running BLAST. Retrying...")
  #  blast_stdout = os.popen4(blastsearch)  ##CAN BE COMMENTED OUT TEMPORARILY FOR FAST TESTING
  #  blast_stdout = blast_stdout[1].read()
  #  if z > 2:
  #    log("Error running Blast, exiting. Please check your system.", exit=True)
  #  z += 1
  args = "-db " + db + " -query input.fasta -outfmt 6 -max_target_seqs " + str(hitspergene) + " -evalue 1e-05 -out input.out -num_threads " + str(nrcpus)
  mgbprocess = Process(target=runblast, args=[args, dbtype])
  mgbprocess.start()
  while True:
      processrunning = "n"
      if mgbprocess.is_alive():
          processrunning = "y"
      if processrunning == "y" and GUI == "y":
          FRAME.update()
      elif processrunning == "y":
          pass
      else:
          break
  try:
      blastoutputfile = open("input.out","r")
  except:
      log("Error while execution NCBI Blast: no output file generated", exit=True)
  blastoutput = blastoutputfile.read()
  blastoutputfile.close()
  return blastoutput
    
def parse_blast(blastoutput, minseqcoverage, minpercidentity, seqlengths, seqdict, dbname, dbtype):
  #Read BLAST output and parse
  log("Blast search finished. Parsing results...")
  blastinfo = blastparse(blastoutput, minseqcoverage, minpercidentity, seqlengths, seqdict, dbname, dbtype)
  blastdict = blastinfo[0]
  if len(blastdict.keys()) == 0:
    log("No BLAST hits found above significance tresholds. Exiting MultiGeneBlast...", exit=True)
  querylist = blastinfo[1]
  return blastdict, querylist

def frame_update():
  if GUI == "y":
    global FRAME
    FRAME.update()

def load_genecluster_info(dbname, allgenomes):
  #Load gene cluster info to memory
  DBPATH = os.environ['BLASTDB']
  clusters = {}
  allgenomes_tags = [genomename[:6] for genomename in allgenomes]
  for i in fileinput.input(DBPATH + os.sep + dbname + "_all_descrs.txt"):
    tabs = i.split("\t")
    if len(tabs) > 0 and (tabs[0] in allgenomes or tabs[0] in allgenomes_tags):
      accession = tabs[0]
      clusterdescription = tabs[1]
      clusters[accession] = clusterdescription
  nucdescriptions = clusters
  frame_update()
  return nucdescriptions, clusters

class ProteinInfo(object):
    __slots__ = ['genome', 'pstart', 'pend', 'strand', 'annotation', 'locustag']
    def __init__(self, genome, pstart, pend, strand, annotation, locustag):
        self.genome = genome
        if pstart.isdigit():
            self.pstart = int(pstart)
        if pend.isdigit():
            self.pend = int(pend)
        if strand == "+":
            self.strand = 1
        else:
            self.strand = -1
        self.annotation = annotation
        self.locustag = locustag

def load_dbproteins_info(querylist, blastdict, dbname):
  ##Load needed gene cluster database proteins info into memory
  global DBPATH
  DBPATH = os.environ['BLASTDB']
  allhitprots = []
  nucdict = {}
  temp_proteininfo = {}
  proteininfo = {}
  for i in querylist:
    if blastdict.has_key(i):
      subjects = blastdict[i][0]
      for j in subjects:
        if j not in allhitprots:
          allhitprots.append(j)
  allhitprots.sort()
  allgenomes = []
  proteininfo_archive = tarfile.open(DBPATH + os.sep + dbname + ".pinfo.tar")
  same = "n"
  next_idx = 0
  for j in allhitprots:
    next_idx += 1
    frame_update()
    if same == "n":
      try:
        infofile = proteininfo_archive.extractfile("proteininfo/" + j[:4].upper() + ".pickle")
        proteininfodict = pickle.load(infofile)
      except:
        log(j[:4].upper() + ".pickle" + " missing in proteininfo tar")
        continue        ##########THIS NEEDS TO BE FIXED IN THE DB TO PREVENT THIS ERROR / MAKE THIS UNNECESSARY
    if not proteininfodict.has_key(j):
      log(j + " missing in proteininfodict")
      continue        ##########THIS NEEDS TO BE FIXED IN THE DB TO PREVENT THIS ERROR / MAKE THIS UNNECESSARY
    pinfo = str(proteininfodict[j])
    #Fix for 'stuttering' of nucl accession number
    if pinfo.count(pinfo.split("|")[0] + "|") > 1 and pinfo.split("|")[0] != "":
      log("Correcting" + j + ":\n" + str(pinfo) + "\n" + pinfo.split("|")[0] + "|" + pinfo.rpartition(pinfo.split("|")[0])[2])
      pinfo = pinfo.split("|")[0] + "|" + pinfo.rpartition(pinfo.split("|")[0])[2]
    #Fix for faulty coordinates
    if not pinfo.split("|")[1].replace("-","").isdigit():
      if pinfo.split("|")[1].replace("-","").replace(")","").isdigit():
        pinfo = pinfo.split("|")[0] + "|" + pinfo.split("|")[1].replace(")","") + "|" + "|".join(pinfo.split("|")[2:])
      if "," in pinfo.split("|")[1].replace("-",""):
        pinfo = pinfo.split("|")[0] + "|" + pinfo.split("|")[1].partition(",")[0] + "-" + pinfo.split("|")[1].rpartition("-")[2] + "|" + "|".join(pinfo.split("|")[2:])
    if "|" not in str(pinfo) or pinfo.split("|")[0] == "":
      log("DNA accession number missing for " + j + "\n" + str(pinfo))
      continue        ##########THIS NEEDS TO BE FIXED IN THE DB TO PREVENT THIS ERROR / MAKE THIS UNNECESSARY
    tabs = pinfo.split("|")
    if len(tabs) < 4:
      log("Faulty info for " + j + ":\n" + str(pinfo))
      continue
    if "-" not in tabs[1] and "-" in tabs[2]:
      del tabs[1]
    protein = tabs[3]
    genome = tabs[0]
    location = tabs[1]
    strand = tabs[2]
    annotation = tabs[4]
    locustag = tabs[5]
    pstart = location.partition("-")[0]
    pend = location.partition("-")[2]
    if not pend.isdigit():
      pend = str(int(pstart) + 100)
    if genome not in allgenomes:
      allgenomes.append(genome)
      temp_proteininfo[genome] = [j, ProteinInfo(genome,pstart,pend,strand,annotation,locustag)]
    else:
      if temp_proteininfo.has_key(genome):
        old_entry = temp_proteininfo[genome]
        proteininfo[old_entry[0]] = old_entry[1]
        nucdict[old_entry[0]]  = old_entry[1].genome
        del temp_proteininfo[genome]
      proteininfo[j] = ProteinInfo(genome,pstart,pend,strand,annotation,locustag)
      nucdict[j]  = genome
    if not (len(allhitprots) > next_idx and j[:4] == allhitprots[next_idx][:4]):
      infofile.close()
      same = "n"
    else:
      same = "y"
    lasthitprot = j
  proteininfo_archive.close()
  allgenomes.sort()
  return allgenomes, nucdict, proteininfo

def load_ndb_info(querylist, blastdict, dbname):
  ##Load needed gene cluster database proteins info into memory
  global DBPATH
  DBPATH = os.environ['BLASTDB']
  allhitprots = []
  nucdict = {}
  proteininfo = {}
  for i in querylist:
    if blastdict.has_key(i):
      subjects = blastdict[i][0]
      for j in subjects:
        if j not in allhitprots:
          allhitprots.append(j)
        genome = j.rpartition("_")[0]
        pstart = min([blastdict[i][1][j][4],blastdict[i][1][j][5]])
        pend = max([blastdict[i][1][j][4],blastdict[i][1][j][5]])
        if int(blastdict[i][1][j][5]) > int(blastdict[i][1][j][4]):
          strand = "+"
        else:
          strand = "-"
        annotation = "tblastn hit"
        locustag = "tblastn_hit_" + j.rpartition("_")[2]
        proteininfo[j] = ProteinInfo(genome,pstart,pend,strand,annotation,locustag)
  allhitprots.sort()
  allgenomes = []
  frame_update()
  for j in allhitprots:
    genome = j.rpartition("_")[0]
    if genome not in allgenomes:
        allgenomes.append(genome)
    nucdict[j] = genome
  allgenomes.sort()
  return allgenomes, nucdict, proteininfo

def load_other_genes(allgenomes, proteininfo, dbname, blastdict):
  #Add all other genes from same genomes to proteininfo
  global DBPATH
  DBPATH = os.environ['BLASTDB']
  genecords_archive = tarfile.open(DBPATH + os.sep + dbname + ".cords.tar")
  same = "n"
  for i in allgenomes:
    frame_update()
    if same == "n":
      infofile = genecords_archive.extractfile("genecords/" + i[:5].upper() + ".pickle")
      genecordsdict = pickle.load(infofile)
    genecords = genecordsdict[i]
    genomeprotpositions = [proteininfo[prot].pstart for prot in proteininfo.keys() if i == proteininfo[prot].genome]
    genomeprotpositions.extend([proteininfo[prot].pend for prot in proteininfo.keys() if i == proteininfo[prot].genome])
    blastsubjects = [item for sublist in [blastdict[key][0] for key in blastdict.keys()] for item in sublist]
    for j in genecords:
      if j.count(j.split("|")[0]) > 1:
        j = j.rpartition(j.split("|")[1]) + j.rpartition(j.split("|")[2])
      if "|" in j:
        tabs = j.split("|")
        if "-" not in tabs[1] and "-" in tabs[2]:
          del tabs[1]
        protein = tabs[3]
        genome = tabs[0]
        location = tabs[1]
        strand = tabs[2]
        annotation = tabs[4]
        locustag = tabs[5]
        pstart = location.partition("-")[0]
        pend = location.partition("-")[2]
        if not pstart.isdigit():
          pstart = 0
          pend = 0
        elif not pend.isdigit():
          pend = str(int(pstart) + 100)
        gene_is_close_to_hit = False
        for position in genomeprotpositions:
          if abs(int(pstart) - position) < 20000:
            gene_is_close_to_hit = True
        if not gene_is_close_to_hit and protein not in blastsubjects:
          if not (len(allgenomes) > (allgenomes.index(i) + 1) and i[:5] == allgenomes[allgenomes.index(i) + 1][:5]):
            infofile.close()
            same = "n"
          else:
            same = "y"
          continue
        if not proteininfo.has_key(protein):
          proteininfo[protein] = ProteinInfo(genome,pstart,pend,strand,annotation,locustag)
        correct = "n"
        z = 0
        while correct == "n" and z < 2:
          try:
            number = 1 + int(proteininfo[protein].pstart)
            number = 1 + int(proteininfo[protein].pend)
            correct = "y"
          except:
            j = j.rpartition(genome)[1] + j.rpartition(genome)[2]
            tabs = j.split("|")
            protein = tabs[3]
            genome = tabs[0]
            location = tabs[1]
            strand = tabs[2]
            annotation = tabs[4]
            locustag = tabs[5]
            proteininfo[protein] = ProteinInfo(genome,location.partition("-")[0],location.partition("-")[2],strand,annotation,locustag)
            z += 1
    if not (len(allgenomes) > (allgenomes.index(i) + 1) and i[:5] == allgenomes[allgenomes.index(i) + 1][:5]):
      infofile.close()
      same = "n"
    else:
      same = "y"
  genecords_archive.close()
  return proteininfo

def load_databases(querylist, blastdict, processnr, dbname, dbtype):
  #Load GenBank positional info into memory
  log("Loading GenBank positional info into memory...")
  if dbtype == "prot":
    allgenomes, nucdict, proteininfo = load_dbproteins_info(querylist, blastdict, dbname)
    proteininfo = load_other_genes(allgenomes, proteininfo, dbname, blastdict)
  else:
    allgenomes, nucdict, proteininfo = load_ndb_info(querylist, blastdict, dbname)
  nucdescriptions, clusters = load_genecluster_info(dbname, allgenomes)
  return nucdescriptions, nucdict, proteininfo

def find_hits_positions(blastdict, proteininfo, querylist):
  #Find db xrefs, start positions, end positions, strand info for each hit and add to blast dict
  log("Finding gene info of all hit genes...")
  frame_update()
  blastdict2 = {}
  for i in querylist:
    if blastdict.has_key(i):
      subjects = blastdict[i][0]
      subjects2 = []
      querydict = blastdict[i][1]
      querydict2 = {}
      for j in subjects:
        #genome_acc = nucdict[j]
        #geneinfo = [start,end,strand,annotation,sequence,accnr,genername]
        #blastparse = [perc_ident,blastscore,perc_coverage,evalue]
        #goal = [subject_genecluster,subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
        #goal now = [subject_start,subject_end,subject_strand,subject_annotation,perc_ident,blastscore,perc_coverage,evalue,locustag]
        if proteininfo.has_key(j):
          oldblastdictinfo = querydict[j]
          if proteininfo[j].strand == 1:
            strand = "+"
          else:
            strand = "-"
          newblastdictinfo = [proteininfo[j].genome, proteininfo[j].pstart, proteininfo[j].pend, strand, proteininfo[j].annotation] + oldblastdictinfo + [j]
          querydict2[j] = newblastdictinfo
          subjects2.append(j)
        else:
          print "WARNING:", j, "accession number not taken into account; data entry invalid."
      blastdict2[i] = [subjects2,querydict2]
  blastdict = blastdict2
  return blastdict

def sort_genes_per_nucleotide(querylist, blastdict, nucdict):
  log("Locating gene clusters on nucleotide scaffolds...")
  frame_update()
  universalquerydict = {}
  for i in querylist:
    if blastdict.has_key(i):
      querydict = blastdict[i][1]
      universalquerydict.update(querydict)
  #Make dictionary that sorts all hit genes per nucleotide scaffold
  sourcedict = {}
  source2dict = {}
  multiplelist = []
  multiplehitlist = []
  nucdictkeys = nucdict.keys()
  nucdictkeys.sort()
  for j in nucdict.keys():
    nucsource = nucdict[j]
    if source2dict.has_key(nucsource):
      if nucsource not in multiplelist:
        multiplehitlist.append(source2dict[nucsource][0])
        multiplelist.append(nucsource)
        if j not in multiplehitlist:
          multiplehitlist.append(j)
        sourcedict[nucsource] = source2dict[nucsource] + [j]
      else:
        sourcedict[nucsource].append(j)
        if j not in multiplehitlist:
          multiplehitlist.append(j)
    else:
      source2dict[nucsource] = [j]
  return sourcedict, multiplehitlist, universalquerydict

def find_geneclusters(sourcedict, universalquerydict, allowedgenedistance, nucdescriptions, proteininfo, dnaseqlength):
  #For every genomic scaffold, find gene clusters based on positions and save gene-cluster links in dictionary
  clusters = {}
  clusterdict = {}
  geneposdict = {}
  hitclusters = []
  extraspace = 20000
  if int(dnaseqlength) / 2 > 20000:
    extraspace = int(int(dnaseqlength) / 2)
  if extraspace > int(allowedgenedistance):
    extraspace = int(allowedgenedistance)
  nuccodes = dict.fromkeys(sourcedict.keys()).keys()
  for i in nuccodes:
    frame_update()
    extragenes = []
    for j in proteininfo.keys():
      if i == proteininfo[j].genome:
        extragenes.append(j)
    scaffoldgenes = sourcedict[i]
    nuc_acc = i
    protstartlocations = []
    protendlocations = []
    for j in scaffoldgenes:
      startpos = int(universalquerydict[j][1])
      endpos = int(universalquerydict[j][2])
      if startpos > endpos:
        startpos, endpos = endpos, startpos
      protstartlocations.append(startpos)
      protendlocations.append(endpos)
      geneposdict[j] = [startpos,endpos]
    protstartlocations.sort()
    protendlocations.sort()
    nrlocations = len(protstartlocations)
    a = 0
    clusterstarts = []
    clusterends = []
    for j in protstartlocations:
      if a == 0:
        cstart = str(int(j) - extraspace)
        if int(cstart) < 0:
          cstart = "0"
        clusterstarts.append(cstart)
        if len(protendlocations) == 1:
          clusterends.append(str(int(protendlocations[a]) + extraspace))
      elif a == nrlocations - 1:
        if j < ((protendlocations[a - 1]) + allowedgenedistance):
          clusterends.append(str(int(protendlocations[a]) + extraspace))
        else:
          cend = str(int(protendlocations[a - 1]) + extraspace)
          clusterends.append(cend)
          cstart = str(int(j) - extraspace)
          if int(cstart) < 0:
            cstart = "0"
          clusterstarts.append(cstart)
          clusterends.append(str(protendlocations[a]))
      else:
        if j > ((protendlocations[a - 1]) + allowedgenedistance):
          clusterends.append(str(int(protendlocations[a - 1]) + extraspace))
          cstart = str(j - extraspace)
          if int(cstart) < 0:
            cstart = "0"
          clusterstarts.append(cstart)
        else:
          pass
      a += 1
    geneclusternumber = 0
    for j in clusterstarts:
      geneclustername = nuc_acc + "_" + str(geneclusternumber)
      if geneclustername not in hitclusters:
        hitclusters.append(geneclustername)
      cstart = int(j)
      cend = int(clusterends[geneclusternumber])
      clustergenes = []
      clustergenesdict = {}
      geneclusternumber += 1
      for k in scaffoldgenes:
        startpos = int(geneposdict[k][0])
        endpos = int(geneposdict[k][1])
        if (startpos >= cstart and startpos <= cend) or (endpos >= cstart and endpos <= cend):
          clusterdict[k] = geneclustername
          clustergenes.append(k)
          clustergenesdict[k] = startpos
      for j in extragenes:
        if (int(proteininfo[j].pstart) >= cstart and int(proteininfo[j].pstart) <= cend) or (int(proteininfo[j].pend) >= cstart and int(proteininfo[j].pend) <= cend):
          if j not in clustergenes:
            clustergenes.append(j)
            clustergenesdict[j] = int(proteininfo[j].pstart)
      nucdescription = ""
      for i in nucdescriptions.keys():
        if i in nuc_acc.split(".")[0]:
          nucdescription = nucdescriptions[i]
          break
      clustergenes = sortdictkeysbyvalues(clustergenesdict)
      clusters[geneclustername] = [clustergenes,nucdescription]
  return clusterdict, geneposdict, hitclusters, clusters

def update_blastdict(blastdict, querylist, clusterdict, multiplehitlist):
  #Update blastdict
  blastdict2 = {}
  frame_update()
  for i in querylist:
    if blastdict.has_key(i):
      subjects = blastdict[i][0]
      querydict = blastdict[i][1]
      querydict2 = {}
      multiplehitlistsubjects = "n"
      for j in subjects:
        if j in multiplehitlist:
          multiplehitlistsubjects = "y"
          geneclustername = clusterdict[j]
          oldblastdictinfo = querydict[j][1:]
          newblastdictinfo = [geneclustername] + oldblastdictinfo
          querydict2[j] = newblastdictinfo
      if multiplehitlistsubjects == "y":
        blastdict2[i] = [subjects,querydict2]
    else:
      blastdict2[i] = [[],{}]
  blastdict = blastdict2
  return blastdict

def find_genomic_loci(blastdict, nucdict, proteininfo, allowedgenedistance, querylist, nucdescriptions, dnaseqlength):
  ##In Genbank files, locate clusters of hits within 20 kb distance using new info from blastdict
  ##Save clusters and sort first based on nr query genes having homologs in them, second based on cumulative BLAST bit score; update blastdict with cluster ID
  #Write dictionary with genes and positions for all nucleotide scaffolds
  blastdict = find_hits_positions(blastdict, proteininfo, querylist)
  sourcedict, multiplehitlist, universalquerydict = sort_genes_per_nucleotide(querylist, blastdict, nucdict)
  clusterdict, geneposdict, hitclusters, clusters = find_geneclusters(sourcedict, universalquerydict, allowedgenedistance, nucdescriptions, proteininfo, dnaseqlength)
  blastdict = update_blastdict(blastdict, querylist, clusterdict, multiplehitlist)
  return blastdict, geneposdict, hitclusters, clusters, multiplehitlist

def score_blast(hitclusters, querylist, blastdict, clusters, multiplehitlist, arch_search, syntenyweight):
  #Score BLAST output on all gene clusters
  #Rank gene cluster hits based on 1) number of protein hits covering >25% sequence length or at least 100aa alignment, with >30% identity and 2) cumulative blast score
  #Find number of protein hits and cumulative blast score for each gene cluster
  log("   Scoring Blast outputs...")
  hitclusterdict = {}
  hitclusterdata = {}
  for i in hitclusters:
    frame_update()
    hitclusterdatalist = []
    nrhits = float(0)
    cumblastscore = float(0)
    hitpositions = []
    hitposcorelist = []
    for j in querylist:
      querynrhits = 0
      querycumblastscore = float(0)
      nrhitsplus = "n"
      if blastdict.has_key(j):
        for k in blastdict[j][0]:
          if k in multiplehitlist and i == blastdict[j][1][k][0]:
            if [querylist.index(j),clusters[i][0].index(blastdict[j][1][k][11])] not in hitpositions:
              nrhitsplus = "y"
              querynrhits += 1
              blastscore = float(blastdict[j][1][k][6]) / 1000000
              querycumblastscore = querycumblastscore + blastscore
              hitclusterdatalist.append([j,k,blastdict[j][1][k][5],blastdict[j][1][k][6],blastdict[j][1][k][7],blastdict[j][1][k][8]])
              hitclusterdata[i] = hitclusterdatalist
              hitpositions.append([querylist.index(j),clusters[i][0].index(blastdict[j][1][k][11])])
        if nrhitsplus == "y":
          nrhits += 1
          for hit in range(querynrhits):
            hitposcorelist.append(0)
          cumblastscore = cumblastscore + float(querycumblastscore)
    query_givenscores_querydict = {}
    query_givenscores_hitdict = {}
    #Find groups of hits
    hitgroupsdict = {}
    for p in hitpositions:
      if not hitgroupsdict.has_key(p[0]):
        hitgroupsdict[p[0]] = [p[1]]
      else:
        hitgroupsdict[p[0]].append(p[1])
    #Calculate synteny score; give score only if more than one hits (otherwise no synteny possible), and only once for every query gene and every hit gene
    if arch_search == "n":
      synteny_score = 0
      z = 1
      if nrhits > 1:
        for p in hitpositions[:-1]:
          tandem = "n"
          #Check if a gene homologous to this gene has already been scored for synteny in the previous entry
          if p[1] in hitgroupsdict[hitpositions[z][0]]:
            tandem = "y"
          #Score entry
          if ((not query_givenscores_querydict.has_key(p[0])) or query_givenscores_querydict[p[0]] == 0) and ((not query_givenscores_hitdict.has_key(p[1])) or query_givenscores_hitdict[p[1]] == 0) and tandem == "n":
            q = hitpositions[z]
            if (abs(p[0] - q[0]) < 2) and abs(p[0]-q[0]) == abs(p[1]-q[1]):
              synteny_score += 1
              if hitposcorelist[z - 1] == 1 or hitposcorelist[z] == 1:
                synteny_score += 1
              query_givenscores_querydict[p[0]] = 1
              query_givenscores_hitdict[p[1]] = 1
            else:
              query_givenscores_querydict[p[0]] = 0
              query_givenscores_hitdict[p[1]] = 0
          z += 1
      #Weigh synteny score by factor
      synteny_score = float(synteny_score) * syntenyweight
      #sorting score is based on number of hits (discrete values) & cumulative blast score (behind comma values)
      sortingscore = nrhits + synteny_score + cumblastscore
    else:
      sortingscore = nrhits + cumblastscore
    hitclusterdict[i] = float(sortingscore)
  #Sort gene clusters
  rankedclusters = sortdictkeysbyvaluesrev(hitclusterdict)
  rankedclustervalues = sortdictkeysbyvaluesrevv(hitclusterdict)
  return rankedclusters, rankedclustervalues, hitclusterdata

def write_txt_output(rankedclusters, rankedclustervalues, hitclusterdata, proteins, proteininfo, querytags, infile, clusters, nucdescriptions, pages):
  global dbname
  #Check if the number of set pages is not too large for the number of results
  possible_pages = len(rankedclusters) / 50
  if len(rankedclusters) % 50 > 1:
    possible_pages += 1
  if possible_pages < int(pages):
    pages = possible_pages
  if pages == 0:
    pages = 1
  #Output for each hit: table of genes and locations of input cluster, table of genes and locations of hit cluster, table of hits between the clusters
  log("   Writing TXT output file...")
  out_file = open(dbname + "clusterblast_output.txt","w")
  out_file.write("ClusterBlast scores for " + infile)
  out_file.write("\n")
  out_file.write("\n")
  out_file.write("Table of genes, locations, strands and annotations of query cluster:")
  out_file.write("\n")
  maxpos = 0
  minpos = 1000000000000
  for i in querytags:
    out_file.write(i)
    out_file.write("\t")
    out_file.write(proteins[3][i][0])
    out_file.write("\t")
    out_file.write(proteins[3][i][1])
    out_file.write("\t")
    out_file.write(proteins[3][i][2])
    out_file.write("\t")
    out_file.write(proteins[3][i][3])
    out_file.write("\t")
    out_file.write(proteins[3][i][-1])
    out_file.write("\t")
    out_file.write("\n")
    if int(proteins[3][i][0]) > maxpos:
      maxpos = int(proteins[3][i][0])
    if int(proteins[3][i][1]) > maxpos:
      maxpos = int(proteins[3][i][1])
    if int(proteins[3][i][0]) < minpos:
      minpos = int(proteins[3][i][0])
    if int(proteins[3][i][1]) < minpos:
      minpos = int(proteins[3][i][1])
  frame_update()
  #Add extra genes from query file that are in between selected genes
  for i in proteins[0]:
    if int(i.split("|")[2].partition("-")[0]) > minpos and int(i.split("|")[2].partition("-")[2]) > minpos and int(i.split("|")[2].partition("-")[0]) < maxpos and int(i.split("|")[2].partition("-")[2]) < maxpos:
      if i.split("|")[4] in querytags or i.split("|")[-1] in querytags:
        continue
      try:
        out_file.write(i.split("|")[4])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[4]][0])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[4]][1])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[4]][2])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[4]][3])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[4]][-1])
        out_file.write("\t")
      except:
        out_file.write(i.split("|")[-1])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[-1]][0])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[-1]][1])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[-1]][2])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[-1]][3])
        out_file.write("\t")
        out_file.write(proteins[3][i.split("|")[-1]][-1])
        out_file.write("\t")
      out_file.write("\n")
  out_file.write("\n")
  out_file.write("\n")
  out_file.write("Significant hits: ")
  out_file.write("\n")
  z = 0
  for i in rankedclusters[:(pages * 50)]:
    out_file.write(str(z+1) + ". " + i + "\t" + clusters[i][1])
    out_file.write("\n")
    z += 1
  out_file.write("\n")
  out_file.write("\n")
  z = 0
  out_file.write("Details:")
  for i in rankedclusters[:(pages * 50)]:
    frame_update()
    value = str(rankedclustervalues[z])
    nrhits = value.split(".")[0]
    if nrhits > 0:
      out_file.write("\n\n")
      out_file.write(">>")
      out_file.write("\n")
      mgbscore = cumblastscore = str(float(value.split(".")[0] + "." + value.split(".")[1][0]))
      cumblastscore = str(int(float("0." + value.split(".")[1][1:]) * 100000))
      out_file.write("\n")
      out_file.write(str(z+1) + ". " + i)
      out_file.write("\n")
      nucleotidename = ""
      for j in i.split("_")[:-1]:
        nucleotidename = nucleotidename + j + "_"
      nucleotidename = nucleotidename[:-1].split(".")[0]
      nucdescription = ""
      for j in nucdescriptions.keys():
        if j in nucleotidename:
          nucdescription = nucdescriptions[j]
          break
      out_file.write("Source: " + nucdescription)
      out_file.write("\n")
      out_file.write("Number of proteins with BLAST hits to this cluster: " + nrhits)
      out_file.write("\n")
      out_file.write("MultiGeneBlast score: " + mgbscore)
      out_file.write("\n")
      out_file.write("Cumulative Blast bit score: " + cumblastscore)
      out_file.write("\n")
      out_file.write("\n")
      out_file.write("Table of genes, locations, strands and annotations of subject cluster:")
      out_file.write("\n")
      clusterproteins = clusters[i][0]
      for j in clusterproteins:
        #if proteinlocations.has_key(j) and proteinannotations.has_key(j) and proteinstrands.has_key(j):
        out_file.write(j)
        out_file.write("\t")
        out_file.write(str(proteininfo[j].pstart))
        out_file.write("\t")
        out_file.write(str(proteininfo[j].pend))
        out_file.write("\t")
        if proteininfo[j].strand == 1:
            out_file.write("+")
        else:
            out_file.write("-")
        out_file.write("\t")
        out_file.write(str(proteininfo[j].annotation))
        out_file.write("\t")
        out_file.write(str(proteininfo[j].locustag.replace("\n","")))
        out_file.write("\n")
      out_file.write("\n")
      out_file.write("Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):")
      out_file.write("\n")
      if i in hitclusterdata.keys():
        tabledata = hitclusterdata[i]
        for x in tabledata:
          w = 0
          for y in x:
            if w == 0:
              out_file.write(str(y).split("|")[4])
              out_file.write("\t")
              w += 1
            else:
              out_file.write(str(y))
              out_file.write("\t")
          out_file.write("\n") 
      else:
        "data not found"
        out_file.write("\n") 
      out_file.write("\n")
      z += 1
  out_file.close()
  return pages

def score_blast_output(hitclusters, querylist, blastdict, multiplehitlist, proteins, proteininfo, querytags, infile, clusters, nucdescriptions, pages, arch_search, syntenyweight):
  rankedclusters, rankedclustervalues, hitclusterdata = score_blast(hitclusters, querylist, blastdict, clusters, multiplehitlist, arch_search, syntenyweight)
  pages = write_txt_output(rankedclusters, rankedclustervalues, hitclusterdata, proteins, proteininfo, querytags, infile, clusters, nucdescriptions, pages)
  return pages

def read_multigeneblast_data(page):
  queryclusterdata = {}
  nrhitgeneclusters = {}
  clusterblastfile = open(dbname + "clusterblast_output.txt","r")
  clusterblastfile = clusterblastfile.read()
  clusterblastfile = clusterblastfile.replace("\r","\n")
  tophitclusters = []
  #Identify top 50 hits for visualization
  hitlines = [i for i in ((clusterblastfile.split("Significant hits: \n")[1]).split("\nDetails:")[0]).split("\n") if i != ""]
  a = 0
  cb_accessiondict = {}
  b = 1
  for i in hitlines:
    if " " in i:
      cb_accessiondict[b] = (i.split("\t")[0]).split(" ")[1]
    b += 1
    if a < page * 50 and a >= (page - 1) * 50:
      if len(i) < 140:
        tophitclusters.append(i)  
      elif len(i) >= 140:
        j = i[0:137] + "..."
        tophitclusters.append(j)
    a += 1
  details = (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:]
  nrhitclusters = len(tophitclusters)
  frame_update()
  #Save query gene cluster data
  querylines = ((clusterblastfile.split("Table of genes, locations, strands and annotations of query cluster:\n")[1]).split("\n\n\nSignificant hits:")[0]).split("\n")
  queryclustergenes = []
  queryclustergenesdetails = {}
  for i in querylines:
    tabs = i.split("\t")
    queryclustergenes.append(tabs[0])
    queryclustergenesdetails[tabs[0]] = [tabs[1],tabs[2],tabs[3],tabs[4], tabs[5]]
  #Sort query cluster genes by start position
  starts = [max([int(queryclustergenesdetails[gene][0]), int(queryclustergenesdetails[gene][1])]) for gene in queryclustergenes]
  genesAndStarts = zip(starts, queryclustergenes)
  genesAndStarts.sort()
  starts, queryclustergenes = zip(*genesAndStarts)
  return queryclusterdata, nrhitgeneclusters, nrhitclusters, cb_accessiondict, queryclustergenes, queryclustergenesdetails, tophitclusters, details

def process_multigeneblast_data(nrhitgeneclusters, nrhitclusters, cb_accessiondict, queryclustergenes, queryclustergenesdetails, internalhomologygroupsdict, tophitclusters, details, page):
  #For every gene cluster, store hit genes and details
  colorgroupsdict = {}
  hitclusterdata = {}
  blastdetails = {}
  mgb_scores = {}
  hitclusternr = 1
  for i in details:
    frame_update()
    hitclustergenes = []
    hitclustergenesdetails = {}
    #Only calculate for specified hit gene clusters
    if not (hitclusternr <= page * 50 and hitclusternr >= (page - 1) * 50):
      hitclusternr += 1
    else:
      nrhitgeneclusters[1] = hitclusternr
      accession = cb_accessiondict[hitclusternr]
      #Store mgb score
      mgbscore = i.partition("MultiGeneBlast score: ")[2].partition("\n")[0]
      cumblastscore = i.partition("Cumulative Blast bit score: ")[2].partition("\n")[0]
      mgb_scores[accession] = [mgbscore, cumblastscore]
      #Store Blast details
      blastdetailslines = [line for line in (i.split("%coverage, e-value):\n")[1]).split("\n") if line != ""]
      for line in blastdetailslines:
        tabs = line.split("\t")
        if not blastdetails.has_key(accession):
          blastdetails[accession] = {}
        if not blastdetails[accession].has_key(tabs[1]):
          blastdetails[accession][tabs[1]] = [[tabs[0], tabs[2], tabs[3], tabs[4], tabs[5]]]
        else:
          blastdetails[accession][tabs[1]].append([tabs[0], tabs[2], tabs[3], tabs[4], tabs[5]])
      #Store Basic Blast output data
      hitclustergeneslines = [line for line in ((i.split("Table of genes, locations, strands and annotations of subject cluster:\n")[1]).split("\n\nTable of Blast hits ")[0]).split("\n") if line != ""]
      locations = []
      for j in hitclustergeneslines:
        tabs = j.split("\t")
        hitclustergenes.append(tabs[0])
        hitclustergenesdetails[tabs[0]] = [tabs[1],tabs[2],tabs[3],tabs[4], tabs[5]]
        locations = locations + [int(tabs[1]),int(tabs[2])]
      #cstart = min(locations)
      #cend = max(locations)
      #print [proteininfo[j][0] for j in proteininfo.keys()]
      #z = 0
      #for j in proteininfo.keys():
      #  if accession.rpartition("_")[0] == proteininfo[j][0]:
      #    z += 1
      #    if z == 200:
      #      z = 0
      #    if int(proteininfo[j][1]) > cstart and int(proteininfo[j][2]) < cend:
      #      #proteininfo[j] = [genome,location.partition("-")[0],location.partition("-")[2],strand,annotation]
      #      if j not in hitclustergenes:
      #        hitclustergenes.append(j)
      #        hitclustergenesdetails[j] = proteininfo[j][1:]
      blasthitslines = [line for line in ((i.split("%coverage, e-value):\n")[1]).split("\n\n")[0]).split("\n") if line != ""]
      if len(blasthitslines) > 0:
        blasthitdict = {}
        blastdetailsdict = {}
        querygenes = []
        revblasthitdict = {}
        hitgenes = []
        for i in blasthitslines:
          tabs = i.split("\t")
          if blasthitdict.has_key(tabs[0]):
            hits = blasthitdict[tabs[0]]
            hits.append(tabs[1])
            blasthitdict[tabs[0]] = hits
            if revblasthitdict.has_key(tabs[1]):
              revhits = revblasthitdict[tabs[1]]
              revhits.append(tabs[0])
              revblasthitdict[tabs[1]] = revhits
            else:
              revblasthitdict[tabs[1]] = [tabs[0]]
            blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[4],tabs[3]]
            if tabs[0] not in querygenes:
              querygenes.append(tabs[0])
            hitgenes.append(tabs[1])
          else:
            blasthitdict[tabs[0]] = [tabs[1]]
            if revblasthitdict.has_key(tabs[1]):
              revhits = revblasthitdict[tabs[1]]
              revhits.append(tabs[0])
              revblasthitdict[tabs[1]] = revhits
            else:
              revblasthitdict[tabs[1]] = [tabs[0]]
            blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[4],tabs[3]]
            if tabs[0] not in querygenes:
              querygenes.append(tabs[0])
            hitgenes.append(tabs[1])
        #Make groups of genes for coloring
        colorgroups = []
        internalgroups = internalhomologygroupsdict[1]
        for i in internalgroups:
          querygenes_and_hits = []
          for j in i:
            #Make list of query gene and its hits
            additionalhits = []
            #For each hit, check if it was also hit by another gene; if so, only add it to the group if this hit had the lowest blast score
            queryscore = 0
            if blasthitdict.has_key(j):
              for k in blasthitdict[j]:
                otherscores = []
                for l in blastdetailsdict.keys():
                  if j == l.partition("_|_")[0] and k == l.rpartition("_|_")[2]:
                    queryscore = blastdetailsdict[l][1]
                  if k in l and j not in l:
                    otherscores.append(blastdetailsdict[l][1])
                allscores = otherscores + [queryscore]
                if int(queryscore) == max([int(m) for m in allscores]):
                  additionalhits.append(k)
               #Add additional hits to the querygenes_and_hits list that will form a colorgroup
              querygenes_and_hits = querygenes_and_hits + additionalhits
              if j not in querygenes_and_hits:
                querygenes_and_hits.append(j)
          if len(querygenes_and_hits) > 0:
            colorgroups.append(querygenes_and_hits)
        colorgroupsdict[hitclusternr] = colorgroups
        hitclusterdata[hitclusternr] = [colorgroupsdict,hitclustergenes,hitclustergenesdetails,queryclustergenes,queryclustergenesdetails,tophitclusters,accession]
        hitclusternr += 1
      else:
        nrhitclusters = nrhitclusters - 1
  if len(details) == 0:
    log("MultiGeneBlast found no significant hits. Exiting...", exit=True)
  return hitclusterdata, nrhitclusters, blastdetails, mgb_scores
    
def write_svg_files(queryclusterdata, hitclusterdata, nrhitclusters, internalhomologygroupsdict, svgfolder, page, screenwidth, arch_search):
  queryclusterdata[1] = [nrhitclusters,hitclusterdata]
  clusterblastpositiondata = {}
  i = page
  #Create alignment svg for each pair of hit&query
  hitclusters = [nr + (page - 1) * 50 for nr in range(queryclusterdata[1][0] + 1)[1:]]
  #Create svgs for pairwise gene cluster alignment
  colorschemedict,rgbcolorscheme = calculate_colorgroups(1,hitclusters,queryclusterdata,internalhomologygroupsdict)
  for k in hitclusters:
    frame_update()
    cresults = clusterblastresults(page,[k],queryclusterdata,colorschemedict,rgbcolorscheme, screenwidth, arch_search)
    s = cresults[0]
    clusterblastpositiondata[str(i) + "_"+str(k)] = cresults[1]
    outfile = open(svgfolder + "clusterblast" + str(i) + "_" + str(k) + ".svg","w")
    outfile.write(s.getXML())
    outfile.close()
  #Create svgs for multiple gene cluster alignment
  cresults = clusterblastresults(page,hitclusters,queryclusterdata,colorschemedict,rgbcolorscheme, screenwidth, arch_search, allhits="y")
  s = cresults[0]
  clusterblastpositiondata[str(i) + "_all"] = cresults[1]
  outfile = open(svgfolder + "clusterblast" + str(i) + "_all.svg","w")
  outfile.write(s.getXML())
  outfile.close()
  return clusterblastpositiondata, colorschemedict

def write_svgs(page, screenwidth, internalhomologygroupsdict, arch_search):
  log("Writing visualization SVGs and XHTML")
  svgfolder = "svg/"
  try:
    os.mkdir(svgfolder)
  except(IOError,OSError):
    pass
  #Read in MultiGeneBlast output data
  queryclusterdata, nrhitgeneclusters, nrhitclusters, cb_accessiondict, queryclustergenes, queryclustergenesdetails, tophitclusters, details = read_multigeneblast_data(page)
  hitclusterdata, nrhitclusters, blastdetails, mgb_scores = process_multigeneblast_data(nrhitgeneclusters, nrhitclusters, cb_accessiondict, queryclustergenes, queryclustergenesdetails, internalhomologygroupsdict, tophitclusters, details, page)
  clusterblastpositiondata, colorschemedict = write_svg_files(queryclusterdata, hitclusterdata, nrhitclusters, internalhomologygroupsdict, svgfolder, page, screenwidth, arch_search)
  return queryclusterdata, colorschemedict, clusterblastpositiondata, blastdetails, mgb_scores

def runmuscle(args):
  os.system("muscle " + args)

def align_muscle(include_muscle, colorschemedict, seqdict):
  #Create Muscle alignments of colourgroups
  musclegroups = []
  if include_muscle == "y":
    log("Aligning homologous sequences with Muscle")
    try:
      os.mkdir("fasta")
    except(IOError,OSError):
      pass
    orthogroupsdup = colorschemedict.values()
    orthogroups = dict.fromkeys(orthogroupsdup).keys()
    for k in orthogroups:
      frame_update()
      accessions = []
      for l in colorschemedict.keys():
        if colorschemedict[l] == k:
          accessions.append(l)
      seqdict2 = {}
      for key in seqdict.keys():
        seqdict2[key.split("|")[-1]] = seqdict[key]
      queryseqs = [">" + acc + "\n" + seqdict2[acc] + "\n" for acc in accessions if seqdict2.has_key(acc)]
      accessions = [acc for acc in accessions if testaccession(acc) == "y"]
      if len(queryseqs) + len(accessions) < 2:
        continue
      musclegroups.append(k)
      url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&tool="multigeneblast"&ID=' + ",".join(accessions) + "&rettype=fasta&retmode=text"
      urltry = "n"
      tries = 0
      while urltry == "n":
        try:
          time.sleep(3)
          req = urllib2.Request(url)
          response = urllib2.urlopen(req)
          output = response.read()
          if len(output) > 5:
            urltry = "y"
          if ">" not in output:
            log("Downloading of FASTA sequences failed")
            break
        except (IOError,httplib.BadStatusLine,URLError,httplib.HTTPException):
          tries += 1
          if tries == 5:
            break
          log("Waiting for connection... (4)")
          time.sleep(60)
      outfile = open("fasta" + os.sep + "orthogroup" + str(k) + ".fasta","w")
      for seq in queryseqs:
        outfile.write(seq)
      outfile.write(output)
      outfile.close()
      args = "-quiet -in fasta" + os.sep + "orthogroup" + str(k) + ".fasta -out fasta" + os.sep + "orthogroup" + str(k) + "_muscle.fasta"
      muscleprocess = Process(target=runmuscle, args=[args])
      muscleprocess.start()
      while True:
          processrunning = "n"
          if muscleprocess.is_alive():
              processrunning = "y"
          if processrunning == "y":
              frame_update()
          elif processrunning == "y":
              pass
          else:
              break
  return musclegroups

def create_xhtml_template(queryclusterdata, page, pages):
  #Create HTML file with gene cluster info in hidden div tags
  htmlfile = open(MGBPATH + os.sep + "empty.xhtml","r")
  html = htmlfile.read()
  html = html.replace("\r","\n")
  htmlparts = html.split("<SPLIT HERE>")
  htmloutfile = open(dbname + "displaypage" + str(page) + ".xhtml","w")
  htmloutfile.write(htmlparts[0] + '  displaycblastresults("' + str(page) + '","all")' + htmlparts[1])
  htmloutfile.write("var list=[" + ",".join([str(nr + (page - 1) * 50 + 1) for nr in range(queryclusterdata[1][0])]) + ",'all'];" + htmlparts[2])
  htmloutfile.write("var list=[" + ",".join([str(nr + (page - 1) * 50 + 1) for nr in range(queryclusterdata[1][0])]) + ",'all'];" + htmlparts[3])
  htmloutfile.write('<a class="bigtext"><br/><br/>&nbsp;Results pages: ')
  for pagenr in [pagenr + 1 for pagenr in range(pages)]:
    htmloutfile.write('<a href="displaypage' + str(pagenr) + '.xhtml" class="bigtext">' + str(pagenr) + '</a>')
    if pagenr != pages:
      htmloutfile.write(", ")
    else:
      htmloutfile.write("</a>")
  htmloutfile.write(htmlparts[4])
  return htmloutfile, htmlparts

def write_xhtml_output(htmloutfile, queryclusterdata, clusters, clusterblastpositiondata, nucname, page, screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, dbtype):
  #Write ClusterBlast divs with pictures and description pop-up tags
  frame_update()
  htmloutfile.write('<div id="clusterblastview" class="clusterdescr">\n\n')
  #Add menu bar 3
  htmloutfile.write('<div id="bartext3" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:3px; left:20px;"><b>MultiGeneBlast hits</b></div>')
  htmloutfile.write('<div id="descrbar3" style="position:absolute; z-index:1; top:0px;"><img src="images/bar.png" height="25" width="' + str(int(0.75*screenwidth)) + '"/></div>')
  htmloutfile.write('<div class="help" id="help3" style="position:absolute; z-index:1; top:2px; left:' + str(int(screenwidth * 0.75) - 30) + 'px;"><a href="http://multigeneblast.sourceforge.net/usage.html" target="_blank"><img border="0" src="images/help.png"/></a></div>')
  qclusternr = page
  nrhitclusters = queryclusterdata[1][0]
  hitclusterdata = queryclusterdata[1][1]
  htmloutfile.write('<div id="qcluster' + str(qclusternr) + '">\n<br/><br/>\n<div align="left">\n<form name="clusterform' + str(qclusternr) + '">\n<select name="selection' + str(qclusternr) + '" onchange="javascript:navigate(this);">\n')
  htmloutfile.write('<option value="">Select gene cluster alignment</option>\n')
  for i in range(nrhitclusters):
    cdescription = hitclusterdata[i + 1 + (page - 1) * 50][5][i].replace("&","&amp;")
    if len(cdescription) > 80:
      cdescription = cdescription[:77] + "..."
    htmloutfile.write('<option value="javascript:displaycblastresults(' + str(page) + ',' + str(i+1 + (page - 1) * 50) + ')">' + cdescription + '</option>\n')
  htmloutfile.write('</select>\n</form>\n\n</div>')
  htmloutfile.write('<div style="position:absolute; top:33px; left:' + str(screenwidth*0.625) + 'px;"><img src="images/button.gif" name="button' + str(qclusternr) + '" onclick="javascript:displaybutton(' + str(qclusternr) + ');"/></div>')
  for i in range(nrhitclusters):
    frame_update()
    hitclusterdata = queryclusterdata[1][1]
    queryclustergenes = hitclusterdata[hitclusterdata.keys()[0]][3]
    queryclustergenesdetails = hitclusterdata[hitclusterdata.keys()[0]][4]
    hitclusternumber =  i + 1 + (page - 1) * 50
    cluster_acc = hitclusterdata[hitclusternumber][6]
    cluster_blastdetails = blastdetails[cluster_acc]
    mgbscore = mgb_scores[cluster_acc][0]
    cumblastscore = mgb_scores[cluster_acc][1]
    hitclustergenes = clusters[cluster_acc][0]
    hitclustergenesdetails = hitclusterdata[hitclusternumber][2]
    relpositiondata = clusterblastpositiondata[str(page) + "_" + str(i + 1 + (page - 1) * 50)]
    qrel_starts = relpositiondata[0][0]
    qrel_ends = relpositiondata[0][1]
    hrel_starts = relpositiondata[1][hitclusternumber][0]
    hrel_ends = relpositiondata[1][hitclusternumber][1]
    strandsbalance = relpositiondata[2][hitclusternumber]
    hstarts = relpositiondata[3][hitclusternumber][0]
    hends = relpositiondata[3][hitclusternumber][1]
    invertedhstarts = [str(100000000 - int(l)) for l in hstarts]
    invertedhends = [str(100000000 - int(l)) for l in hends]
    if strandsbalance < 0:
      hitclustergenes.reverse()
    htmloutfile.write('<div id="hitcluster' + str(qclusternr) + '_' + str(i + 1 + (page - 1) * 50) + '">\n')
    #Load svg and embed it into XHTML
    svglines = open("svg" + os.sep + "clusterblast" + str(qclusternr) + '_' + str(i + 1 + (page - 1) * 50) + ".svg","r").read().split("\n")
    htmloutfile.write("\n" + svglines[0][:-1] + 'id="svg' + str(qclusternr) + '_' + str(i + 1 + (page - 1) * 50) + '" >' + "\n")
    for svgline in svglines[1:]:
      htmloutfile.write(svgline + "\n")
    #Insert gene cluster descriptions
    cgbkdescription = hitclusterdata[i + 1 + (page - 1) * 50][5][i].replace("&","&amp;").replace("\t"," ").partition(" ")[2].partition(" ")[2].split(", whole")[0].split(", complete")[0].split(", partial")[0]
    if len(cgbkdescription) > 90:
      cgbkdescription = cgbkdescription[:87] + "..."
    if testaccession(cluster_acc.rpartition("_")[0]) == "y":
      cdescription = '<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.rpartition("_")[0] + '" target="_blank"> ' + cluster_acc.rpartition("_")[0] + "</a>" + " : " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp;Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
    else:
      cdescription = cluster_acc.rpartition("_")[0] + " : " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp;Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
    if len(nucname) < 90:
      qdescription = "Query: " + nucname
    else:
      qdescription = "Query: " + nucname[0:87] + "..."
    htmloutfile.write('<div id="descriptionquery" style="text-align:left; position:absolute; top:60px; left:10px; font-size:10px; font-style:italic">' + qdescription + '</div>\n')
    htmloutfile.write('<div id="description' + str(qclusternr) + '" style="text-align:left; position:absolute; top:115px; left:10px; font-size:10px; font-style:italic">' + cdescription + '</div>\n')
    #Insert NCBI links
    htmloutfile.write('<div id="pub_pics" style="position:absolute; top:175px; left:' + str(int(screenwidth * 0.0)) + 'px; font-size:10px"> Hit cluster cross-links: \n')
    htmloutfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.rpartition("_")[0] + '" target="_blank"><img align="absmiddle" border="0" src="images/genbank.gif"/></a>\n')
    htmloutfile.write('</div>\n\n')
    #Create gene pop-ups
    a = 0
    for j in queryclustergenes:
      j_accession = j
      htmloutfile.write('<div id="q' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(100) + 'px; left:' + str(int(float(qrel_starts[a])*0.875)) + 'px;">\n')
      htmloutfile.write(queryclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
      link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j_accession + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
      if j != queryclustergenesdetails[j][4] and testaccession(j) == "y":
        htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
      htmloutfile.write("<br/>Location: " + str(queryclustergenesdetails[j][0]) + "-" + str(queryclustergenesdetails[j][1]) + "\n")
      htmloutfile.write("</div>\n\n")
      htmloutfile.write('<div id="q' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(75) + 'px; left:' + str(int(float((float(qrel_starts[a])+float(qrel_ends[a]))/2)*0.9375)) + 'px;">\n')
      if queryclustergenesdetails[j][4] != "" and queryclustergenesdetails[j][4] != "no_locus_tag":          
        htmloutfile.write(queryclustergenesdetails[j][4])
      else:
        htmloutfile.write(j)
      htmloutfile.write("</div>\n\n")        
      a+= 1
    a = 0
    for j in hitclustergenes:
      if ((hitclustergenesdetails[j][0] in hstarts or hitclustergenesdetails[j][0] in hends) and (hitclustergenesdetails[j][1] in hends or hitclustergenesdetails[j][1] in hstarts)) or ((hitclustergenesdetails[j][1] in invertedhstarts or hitclustergenesdetails[j][1] in invertedhends) and (hitclustergenesdetails[j][0] in invertedhends or hitclustergenesdetails[j][0] in invertedhstarts)):
        j_accession = j
        htmloutfile.write('<div id="h' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(151) + 'px; left:' + str(int(float(hrel_starts[a])*0.875)) + 'px;">\n')
        htmloutfile.write(hitclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
        link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j_accession + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
        if dbtype == "nucl":
          htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/nuccore/' + j.rpartition("_")[0] + '" target="_blank">' + j.rpartition("_")[0] + "</a>\n")
        elif j != hitclustergenesdetails[j][4] and testaccession(j) == "y":
          htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
        htmloutfile.write("<br/>Location: " + str(hitclustergenesdetails[j][0]) + "-" + str(hitclustergenesdetails[j][1]) + "\n")
        if cluster_blastdetails.has_key(j):
          for blasthit in cluster_blastdetails[j]:
            htmloutfile.write("<br/><br/><b>BlastP hit with " + blasthit[0] + "</b>\n<br/>Percentage identity: " + blasthit[1] + " %\n")
            htmloutfile.write("<br/>BlastP bit score: " + blasthit[2] + "\n<br/>Sequence coverage: " + blasthit[3].partition(".")[0] + " %\n")
            htmloutfile.write("<br/>E-value: " + blasthit[4] + "\n<br/>")
        if testaccession(j) == "y" and dbtype != "nucl":
          htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
        if colorschemedict.has_key(j) and colorschemedict[j] in musclegroups:
          htmloutfile.write("<br/><a href=\"fasta" + os.sep + "orthogroup" + str(colorschemedict[j]) + "_muscle.fasta\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n")
        htmloutfile.write("</div>\n\n")
        htmloutfile.write('<div id="h' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(126) + 'px; left:' + str(int(float((float(hrel_starts[a])+float(hrel_ends[a]))/2)*0.9375)) + 'px;">\n')
        if hitclustergenesdetails[j][4] != "" and hitclustergenesdetails[j][4] != "no_locus_tag":          
          htmloutfile.write(hitclustergenesdetails[j][4])
        else:
          htmloutfile.write(j)
        htmloutfile.write("</div>\n\n")    
        a += 1
    htmloutfile.write('</div>\n')
  #Find new relative positions for display of all gene clusters in one picture
  relpositiondata = clusterblastpositiondata[str(page) + "_all"]
  if len(relpositiondata[0]) > 0:
    qrel_starts = relpositiondata[0][0]
    qrel_ends = relpositiondata[0][1]
    htmloutfile.write('<div id="hitcluster' + str(page) + '_all" style="display:none">\n')
    #Load svg and embed it into XHTML
    svglines = open("svg" + os.sep + "clusterblast" + str(qclusternr) + "_all.svg","r").read().split("\n")
    htmloutfile.write("\n" + svglines[0][:-1] + 'id="svg' + str(qclusternr) + '_all" >' + "\n")
    for svgline in svglines[1:]:
      htmloutfile.write(svgline + "\n")
    if len(nucname) < 90:
      qdescription = "Query: " + nucname
    else:
      qdescription = "Query: " + nucname[0:87] + "..."
    htmloutfile.write('<div id="descriptionquery" style="text-align:left; position:absolute; top:60px; left:10px; font-size:10px; font-style:italic">' + qdescription + '</div>\n')
    for i in range(nrhitclusters):
      frame_update()
      hitclusterdata = queryclusterdata[1][1]
      queryclustergenes = hitclusterdata[hitclusterdata.keys()[0]][3]
      queryclustergenesdetails = hitclusterdata[hitclusterdata.keys()[0]][4]
      hitclusternumber =  i + 1 + (page - 1) * 50
      hrel_starts = relpositiondata[1][hitclusternumber][0]
      hrel_ends = relpositiondata[1][hitclusternumber][1]
      cluster_acc = hitclusterdata[hitclusternumber][6]
      cluster_blastdetails = blastdetails[cluster_acc]
      mgbscore = mgb_scores[cluster_acc][0]
      cumblastscore = mgb_scores[cluster_acc][1]
      hitclustergenes = clusters[cluster_acc][0]
      hitclustergenesdetails = hitclusterdata[hitclusternumber][2]
      strandsbalance = relpositiondata[2][hitclusternumber]
      hstarts = relpositiondata[3][hitclusternumber][0]
      hends = relpositiondata[3][hitclusternumber][1]
      invertedhstarts = [str(100000000 - int(l)) for l in hstarts]
      invertedhends = [str(100000000 - int(l)) for l in hends]
      cgbkdescription = hitclusterdata[i + 1 + (page - 1) * 50][5][i].replace("&","&amp;").replace("\t"," ").partition(" ")[2].partition(" ")[2].split(", whole")[0].split(", complete")[0].split(", partial")[0]
      if len(cgbkdescription) > 90:
        cgbkdescription = cgbkdescription[:87] + "..."
      if testaccession(cluster_acc.rpartition("_")[0]) == "y":
        cdescription = str(i+1 + (page - 1) * 50) + ". : " + '<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.rpartition("_")[0] + '" target="_blank"> ' + cluster_acc.rpartition("_")[0] + "</a> " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp; Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
      else:
        cdescription = str(i+1 + (page - 1) * 50) + ". : " + cluster_acc.rpartition("_")[0] + " " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp; Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
      htmloutfile.write('<div id="description' + str(qclusternr) + '" style="text-align:left; position:absolute; top:' + str(int(63 + (51.7 * (hitclusternumber - (page - 1) * 50)))) + 'px; left:10px; font-size:10px; font-style:italic">' + cdescription + '</div>\n')
      if hitclusternumber == 1 + (page - 1) * 50:
        a = 0
        for j in queryclustergenes:
          htmloutfile.write('<div id="all_' + str(qclusternr) + "_0_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(100) + 'px; left:' + str(int(float(qrel_starts[a])*0.875)) + 'px; z-index:2;">\n')
          htmloutfile.write(queryclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
          link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
          if j != queryclustergenesdetails[j][4] and testaccession(j) == "y":
            htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
          htmloutfile.write("<br/>Location: " + str(queryclustergenesdetails[j][0]) + "-" + str(queryclustergenesdetails[j][1]) + "\n")
          if testaccession(j) == "y":
            htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
          if colorschemedict.has_key(j) and colorschemedict[j] in musclegroups:
            htmloutfile.write("<br/><a href=\"fasta" + os.sep + "orthogroup" + str(colorschemedict[j]) + "_muscle.fasta\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n")
          htmloutfile.write("</div>\n\n")
          htmloutfile.write('<div id="all_' + str(qclusternr) + "_0_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(75) + 'px; left:' + str(int(float((float(qrel_starts[a])+float(qrel_ends[a]))/2)*0.9375)) + 'px;">\n')
          if queryclustergenesdetails[j][4] != "" and queryclustergenesdetails[j][4] != "no_locus_tag":          
            htmloutfile.write(queryclustergenesdetails[j][4])
          else:
            htmloutfile.write(j)
          htmloutfile.write("</div>\n\n")        
          a+= 1
      a = 0
      for j in hitclustergenes:
        if ((hitclustergenesdetails[j][0] in hstarts or hitclustergenesdetails[j][0] in hends) and (hitclustergenesdetails[j][1] in hends or hitclustergenesdetails[j][1] in hstarts)) or ((hitclustergenesdetails[j][1] in invertedhstarts or hitclustergenesdetails[j][1] in invertedhends) and (hitclustergenesdetails[j][0] in invertedhends or hitclustergenesdetails[j][0] in invertedhstarts)):
          htmloutfile.write('<div id="all_' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(int(100 + 51.7 * (hitclusternumber - (page - 1) * 50))) + 'px; left:' + str(int(float(hrel_starts[a])*0.875)) + 'px; z-index:2;">\n')
          htmloutfile.write(hitclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
          link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
          if dbtype == "nucl":
            htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/nuccore/' + j.rpartition("_")[2] + '" target="_blank">' + j.rpartition("_")[2] + "</a>\n")
          elif j != hitclustergenesdetails[j][4] and testaccession(j) == "y":
            htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
          htmloutfile.write("<br/>Location: " + str(hitclustergenesdetails[j][0]) + "-" + str(hitclustergenesdetails[j][1]) + "\n")
          if cluster_blastdetails.has_key(j):
            for blasthit in cluster_blastdetails[j]:
              htmloutfile.write("<br/><br/><b>BlastP hit with " + blasthit[0] + "</b>\n<br/>Percentage identity: " + blasthit[1] + " %\n")
              htmloutfile.write("<br/>BlastP bit score: " + blasthit[2] + "\n<br/>Sequence coverage: " + blasthit[3].partition(".")[0] + " %\n")
              htmloutfile.write("<br/>E-value: " + blasthit[4] + "\n<br/>")
          if testaccession(j) == "y" and dbtype != "nucl":
            htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
          if colorschemedict.has_key(j) and colorschemedict[j] in musclegroups:
            htmloutfile.write("<br/><a href=\"fasta" + os.sep + "orthogroup" + str(colorschemedict[j]) + "_muscle.fasta\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n")
          htmloutfile.write("</div>\n\n")
          htmloutfile.write('<div id="all_' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(int(75 + 51.7 * (hitclusternumber - (page - 1) * 50))) + 'px; left:' + str(int(float((float(hrel_starts[a])+float(hrel_ends[a]))/2)*0.9375)) + 'px;">\n')
          if hitclustergenesdetails[j][4] != "" and hitclustergenesdetails[j][4] != "no_locus_tag":          
            htmloutfile.write(hitclustergenesdetails[j][4])
          else:
            htmloutfile.write(j)
          htmloutfile.write("</div>\n\n")    
          a += 1
    htmloutfile.write('</div>\n')
    htmloutfile.write('</div>\n\n')
  else:
    htmloutfile.write('<br/>No homologous gene clusters found.</div>\n')
  htmloutfile.write('</div>\n')
  htmloutfile.write('<div id="creditsbar' + str(i) + '" class="banner" style="position:absolute; width:' + str(int(0.98 * screenwidth)) +'px; align:\'left\'; height:75; top:2750px; left:0px; color:#000066; z-index:-1;">')
  htmloutfile.write('<div style="float:center; font-size:0.9em;">\n<div style="position:absolute; top:0px; left:30px;">\n<img src="images/ruglogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n<img src="images/gbblogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n</div>\n<div style="position:absolute; top:10px; left:340px;">\nDetecting sequence homology at the gene cluster level with MultiGeneBlast.\n<br/>Marnix H. Medema, Rainer Breitling &amp; Eriko Takano (2013)\n<br/><i>Molecular Biology and Evolution</i> , 30: 1218-1223.\n</div>\n</div>\n</div>')
  
def finalize_xhtml(htmloutfile, htmlparts):
  #Add final part of HTML file
  htmloutfile.write(htmlparts[-1])
  #Copy accessory files for HTML viewing
  #if sys.platform == ('win32'):
  #  copycommand1 = "copy/y vis\\* " + genomename + " > nul"
  #  copycommand2 = "copy/y vis\\html\\* " + genomename + "\\html > nul"
  #  copycommand3 = "copy/y vis\\images\\* " + genomename + "\\images > nul"
  #elif sys.platform == ('linux2'):
  #  copycommand1 = "cp vis/* " + genomename + " > /dev/null"
  #  copycommand2 = "cp -r vis/html " + genomename + "/html > /dev/null"
  #  copycommand3 = "cp -r vis/images " + genomename + "/images > /dev/null"
  #os.system(copycommand1)
  #os.system(copycommand2)
  #os.system(copycommand3)

  #Close open html file
  htmloutfile.close()

def create_xhtml_file(queryclusterdata, clusters, clusterblastpositiondata, nucname, page, pages, screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, dbtype):
  htmloutfile, htmlparts = create_xhtml_template(queryclusterdata, page, pages)
  write_xhtml_output(htmloutfile, queryclusterdata, clusters, clusterblastpositiondata, nucname, page, screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, dbtype)
  finalize_xhtml(htmloutfile, htmlparts)

def move_outputfiles(foldername, pages):
  global MGBPATH
  #Move output files to specified output folder. Overwriting when files/folders are already present there
  try:
    os.mkdir(foldername)
  except:
    pass
  try:
    shutil.rmtree(foldername + os.sep + "svg")
  except:
    pass
  try:
    shutil.rmtree(foldername + os.sep + "fasta")
  except:
    pass
  for page in range(pages):
    try:
      os.remove(foldername + os.sep + dbname + "displaypage" + str(page + 1) + ".xhtml")
    except:
      pass
    try:
      shutil.move(dbname + "displaypage" + str(page + 1) + ".xhtml", foldername + os.sep + dbname + "displaypage" + str(page + 1) + ".xhtml")
    except:
      pass
  filestomove = [dbname + "clusterblast_output.txt", "svg", "fasta"]
  frame_update()
  for f in filestomove:
    try:
      os.remove(foldername + os.sep + f)
    except:
      try:
        shutil.rmtree(foldername + os.sep + f)
      except:
        pass
    try:
      shutil.move(f, foldername + os.sep + f)
    except:
      pass
  frame_update()
  filestocopy = ["style.css", "jquery.svg.js", "jquery-1.4.2.min.js", "jquery.svgdom.js"]
  for f in filestocopy:
    try:
      os.remove(foldername + os.sep + f)
    except:
      pass
    shutil.copy(MGBPATH + os.sep + f, foldername + os.sep + f)
  folderstocopy = ["images"]
  for f in folderstocopy:
    try:
      shutil.rmtree(foldername + os.sep + f)
    except:
      pass
    shutil.copytree(MGBPATH + os.sep + f, foldername + os.sep + f)


def main():
  global GUI
  global TEMP
  global MGBPATH
  os.environ['BLASTDB'] = MGBPATH
  
  os.chdir(TEMP)
  GUI = "n"
  starttime = time.time()
  
  opts = Options()
  #Step 1: parse options
  parse_options(sys.argv, opts)
  dbname = opts.db
  global dbname
  print "Step 1/11: Time since start: " + str((time.time() - starttime))
  #Step 2: Read GBK / EMBL file, select genes from requested region and output FASTA file
  proteins, genomic_accnr, dnaseqlength, nucname, querytags, names, seqs, seqdict, arch_search = read_input_file(opts.infile, opts.startpos, opts.endpos, opts.ingenes, opts.gui)
  print "Step 2/11: Time since start: " + str((time.time() - starttime))
  #Step 3: Run internal BLAST
  internalhomologygroupsdict, seqlengths = internal_blast(opts.minseqcov, opts.minpercid, names, proteins, seqdict, opts.nrcpus)
  print "Step 3/11: Time since start: " + str((time.time() - starttime))
  #Step 4: Run BLAST on genbank_mf database
  blastoutput = db_blast(names, seqs, opts.db, opts.nrcpus, opts.hitspergene, opts.dbtype)
  print "Step 4/11: Time since start: " + str((time.time() - starttime))
  #Step 5: Parse BLAST output
  blastdict, querylist = parse_blast(blastoutput, opts.minseqcov, opts.minpercid, seqlengths, seqdict, opts.db, opts.dbtype)
  print "Step 5/11: Time since start: " + str((time.time() - starttime))
  #Step 6: Load genomic databases into memory
  nucdescriptions, nucdict, proteininfo = load_databases(querylist, blastdict, opts.nrcpus, opts.db, opts.dbtype)
  print "Step 6/11: Time since start: " + str((time.time() - starttime))
  #Step 7: Locate Blast hits in genomes
  blastdict, geneposdict, hitclusters, clusters, multiplehitlist = find_genomic_loci(blastdict, nucdict, proteininfo, opts.distancekb, querylist, nucdescriptions, dnaseqlength)
  print "Step 7/11: Time since start: " + str((time.time() - starttime))
  #Step 8: Score Blast output on all loci
  opts.pages = score_blast_output(hitclusters, querylist, blastdict, multiplehitlist, proteins, proteininfo, querytags, opts.infile, clusters, nucdescriptions, opts.pages, arch_search, opts.syntenyweight)
  print "Step 8/11: Time since start: " + str((time.time() - starttime))
  #Output. From here, iterate for every page
  for page in [pagenr + 1 for pagenr in range(opts.pages)]:

    #Step 9: Write MultiGeneBlast SVGs
    queryclusterdata, colorschemedict, clusterblastpositiondata, blastdetails, mgb_scores = write_svgs(page, opts.screenwidth, internalhomologygroupsdict, arch_search)
    print "Step 9/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime))

    #Step 10: Create muscle alignments
    musclegroups = align_muscle(opts.muscle, colorschemedict, seqdict)
    print "Step 10/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime))

    #Step 11: Create XHTML output file
    create_xhtml_file(queryclusterdata, clusters, clusterblastpositiondata, nucname, page, opts.pages, opts.screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, opts.dbtype)
    print "Step 11/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime))

  #Move all files to specified output folder
  move_outputfiles(opts.outputfolder, opts.pages)  

  #Close log file
  print "MultiGeneBlast successfully finished in " + str((time.time() - starttime)) + " seconds.\n"

if __name__ == '__main__':
  freeze_support()
  main()
  
