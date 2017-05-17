#!/usr/bin/env python

"""
This script is used to compute the RMS between two vectors in two differnt files
"""

import sys
import math

def getMag(vec1,vec2):

	try:
		vec1txt=open(vec1,'r')
	except IOError as e:
		return 'Could not open vector file 1'
		
	try:
		vec2txt=open(vec2,'r')
	except IOError as e:
		return 'Could not open vector file 2'
		
	vec1=map(float,vec1txt.read().split(' '))
	vec2=map(float,vec2txt.read().split(' '))
	
	return(sum([(vec1[i]-vec2[i])**float(2.0) for i in range(len(vec1))]))
	

if len(sys.argv) == 3:
	lastout=getMag(sys.argv[1],sys.argv[2])
	print lastout
if len(sys.argv) == 4:
	lastout=getMag(sys.argv[1],sys.argv[2])
	if float(sys.argv[3]) > lastout:
		sys.exit(1)
	else:
		sys.exit(0)
	