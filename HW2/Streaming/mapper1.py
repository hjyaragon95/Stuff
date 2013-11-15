#!/usr/bin/env python
def upperNum(num):
	ans = int(num*10)*0.1
	if (ans < num):
	 ans += 0.1
	return ans
def num2Bin(num):
	upper = upperNum(num)
	lower = upper - 0.1
	return "{:.1f},{:.1f}".format(lower, upper)
import sys

#input comes from STDIN (standard input)
for line in sys.stdin:
	#remove leading and trailing whitespace
	line = line.strip()
	#split the line into words
	x,y = line.split()
	#convert x,y (currently a string) to float
	x = float(x)
	y = float(y)
	binKey = num2Bin(x) + "," + num2Bin(y)
	print '%s\t%s' % (binKey, 1)
