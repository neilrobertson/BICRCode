#!/usr/bin/env python

fh = open("tester.txt", "ab")

fh2 = open("fhtester.py", "rb")

output = fh2.read()

fh.write(output)
