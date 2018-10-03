#!/usr/bin/python

mylist = {}

try:
	mylist["hello"]
except KeyError:
	print "key error"
except:
	print "other"
