#!/usr/bin/python

# use the csv library to process csv files
import csv

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

infile = "2323_9_31_38_K562.csv"
outfile = "2323_9_31_38_K562_aligned.csv"

# open the input file
csvreader = csv.reader(open(infile, "r"))

# open a file in which to dump the results
csvwriter = csv.writer(open(outfile, "w"), dialect='excel')

# get the first line
firstline = csvreader.next()

# headers for three columns. Each column is a list of ID/value pairs.
col1header = firstline[1]
col2header = firstline[4]
col3header = firstline[7]

# to hold the id/value pairs for each column
col1entries = {}
col2entries = {}
col3entries = {}

#breakpoint()
# for remaining lines
for row in csvreader:
        # we could do all these assignments in one line, but this is clearer
        try:
                col1id = row[0]
                col1data = row[1]
        except IndexError:
                print "no entry in first column:"
                print row
        try:
                col2id = row[3]
                col2data = row[4]
        except IndexError:
                print "no entry in second column:"
                print row
        try:
                col3id = row[6]
                col3data = row[7]
        except IndexError:
                print "no entry in third column:"
                print row
        if col1id != "":
                col1entries[col1id] = col1data
        if col2id != "":
                col2entries[col2id] = col2data
        if col3id != "":
                col3entries[col3id] = col3data

# output based on sorted ids
col1ids = col1entries.keys()
col1ids.sort()

# now we're ready to output. Start with the header line. Don't forget the blank spacers!
csvwriter.writerow([ "", col1header, "", "", col2header, "", "", col3header ])

# now align everything with column1 ids
for c1id in col1ids:
        c1data = col1entries[c1id]
        row = [ c1id, c1data, "" ]
        # add column 2 entry, if present
        if c1id in col2entries:
                c2data = col2entries[c1id]
                row.extend([ c1id, c2data, "" ])
                # if we found it, remove so we don't output it again later
                del col2entries[c1id]
        else:
                row.extend([ "", "", "" ])
        # add column 3 entry, if present
        if c1id in col3entries:
                c3data = col3entries[c1id]
                row.extend([ c1id, c3data, "" ])
                # if we found it, remove so we don't output it again later
                del col3entries[c1id]
        else:
                row.extend([ "", "", "" ])
        # write out the row we've built
        csvwriter.writerow(row)
        
# if there are any left in col2...
col2ids = col2entries.keys()
col2ids.sort()

for c2id in col2ids:
        # obviously there was nothing in col1, so add a spacer
        row = [ "", "", "" ]
        c2data = col2entries[c2id]
        row.extend([ c2id, c2data, "" ])
        # add column 3 entry, if present
        if c2id in col3entries:
                c3data = col3entries[c2id]
                row.extend([ c2id, c3data, "" ])
                # if we found it, remove so we don't output it again later
                del col3entries[c2id]
        else:
                row.extend([ "", "", "" ])
        # write out the row we've built
        csvwriter.writerow(row)
        
# if there are any left in col3...
col3ids = col3entries.keys()
col3ids.sort()

for c3id in col3ids:
        # obviously there was nothing in col1 or col2, so add a spacer
        row = [ "", "", "", "", "", "" ]
        c3data = col3entries[c3id]
        row.extend([ c3id, c3data, "" ])
        # write out the row we've built
        csvwriter.writerow(row)
        
