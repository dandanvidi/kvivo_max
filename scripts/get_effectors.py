#!/usr/bin/python

import sys
import re
import argparse
import csv
from brenda_soap import BrendaSoap

parser = argparse.ArgumentParser(description='Download data from BRENDA.')
parser.add_argument('-i', '--ifile', nargs='?', type=argparse.FileType('r'), required=True)
parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'), required=True)
args = parser.parse_args()

brenda = BrendaSoap("dan.davidi@weizmann.ac.il", "300004686")

#allEcNumbers = brenda.GetEcNumbersFromActivatingCompound()


fEcNumbers = args.ifile
ecNumbers = map(lambda k:k.strip(), fEcNumbers.readlines())
fEcNumbers.close()

excl = re.compile("!")

row_dicts = []
for ec in ecNumbers:
    sys.stderr.write('%10s ... ' % ec)
    row_dicts += brenda.GetActivatingCompound(ec)
    row_dicts += brenda.GetInhibitors(ec)
    sys.stderr.write('done\n')

fieldnames = ['ecNumber', 'organism', 'activatingCompound', 'inhibitors']
for d in row_dicts:
    for f in d.keys():
        if f not in fieldnames:
            fieldnames.append(f)

fTurnover = csv.writer(args.outfile, delimiter='\t')
fTurnover.writerow(fieldnames)
for d in row_dicts:
   fTurnover.writerow(map(lambda k: d.get(k, '').encode('utf-8'), fieldnames))
   
args.outfile.close()
sys.stderr.write('All Done !\n')




