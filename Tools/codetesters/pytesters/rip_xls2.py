#!/usr/bin/python

import xlrd
import sys
import os
import csv
import time

def dumpOneTab(tab, name):
    fh = open(name, "w")
    for rownum in range(tab.nrows):
        line = []
        for colnum in range(tab.ncols):
            cellval = tab.cell_value(rowx=rownum, colx=colnum)
            try:
                cellval = str(cellval)
            except:
                cellval = "--garbled--"
            line += [cellval]
        print >> fh, ",".join(line)
    fh.close()

def dumpOneSheet(sheet, outfname, encoding):
    fh = open(outfname, "wb")
    writer = csv.writer(fh)
    writerow = writer.writerow
    ctype_text = xlrd.XL_CELL_TEXT
    row_types = sheet.row_types
    row_values = sheet.row_values
    for rownum in xrange(sheet.nrows):
        writerow([
            cval.encode(encoding) if ctype == ctype_text else cval
            for (ctype, cval) in zip(row_types(rownum), row_values(rownum))
            ])
    fh.close()


if __name__ == "__main__":

    filename = sys.argv[1]
    filebase = os.path.basename(filename)[:-4]
    which = sys.argv[2]
    assert which in ["old", "new"]
    if which == "new":
        try:
            enc = sys.argv[3]
        except IndexError:
            enc = 'utf8'
    try:
        t0 = time.clock()
        import psyco
        psyco.full()
        t1 = time.clock()
        print "psyco setup time: %.2f seconds" % (t1 - t0)
    except ImportError:
        pass
    print "reading workbook..."
    t0 = time.clock()
    xls = xlrd.open_workbook(filename)
    t1 = time.clock()
    print "open time: %.2f seconds" % (t1 - t0)
    t0 = time.clock()
    if which == "old":
        for sheetname in xls.sheet_names():
            dumpfile = "%s-%s.csv" % (filebase, sheetname)
            print "dumping tab", dumpfile, "..."
            sheet = xls.sheet_by_name(sheetname)
            dumpOneTab(sheet, dumpfile)
    else:
        for sheetname in xls.sheet_names():
            dumpfile = "%s-%s-new.csv" % (filebase, sheetname)
            print "dumping sheet", dumpfile, "..."
            sheet = xls.sheet_by_name(sheetname)
            dumpOneSheet(sheet, dumpfile, enc)
    t1 = time.clock()
    print "dump time: %.2f seconds" % (t1 - t0)
