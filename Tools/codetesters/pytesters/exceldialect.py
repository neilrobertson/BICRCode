#!/usr/bin/env python

import csv

writer = csv.writer(open("tester.txt", "w"), dialect="excel")

writer.writerow(["tester", "another", "a third"])
