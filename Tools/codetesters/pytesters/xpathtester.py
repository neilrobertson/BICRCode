#!/usr/bin/env python

from lxml import etree
import pxdom


filename = "base.xml"
dom = pxdom.parse(filename)
tree = etree.parse(filename)

root = dom._childNodes[0]

path = "/base/element1"


nodes = tree.xpath(path)

root.appendChild(nodes[0])

print dom.pxdomContent
