#!/usr/bin/env python

mydict = { "a" : 100, "b" : 50, "c" : 150 }

mylist = [ "a", "b", "c" ]

mylist.sort(key=mydict.__getitem__)

print mylist
