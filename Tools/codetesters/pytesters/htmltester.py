#!/usr/bin/python

import htmllib, formatter
class x(htmllib.HTMLParser):
	def dump(self, tag, attrs):
		print tag,
		for a, v in attrs:
			if a in ['action', 'src', 'href']:
				print a, v,
		print
	def do_img(self, attrs):
		self.dump('img', attrs)
	def start_a(self, attrs):
		self.dump('a', attrs)
	def start_form(self, attrs):
		self.dump('form', attrs)

y = x(formatter.NullFormatter())
y.feed(open('HISTONE.html').read())
y.close() 
