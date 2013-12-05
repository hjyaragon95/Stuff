#!/usr/bin/env python
import sys

# input comes from STDIN (standard input)
for line in sys.stdin:
    # Need to remove trailing '\n'
    line = line.strip()
    # split the line into (x,y)
    x,y = line.split('\t')
    x = float(x)
    y = float(y)
    if x>0:
        x_lo=float(int(x*10))/10
        x_hi=x_lo+0.1
    elif x<0:
        x_hi=float(int(x*10))/10
        x_lo=x_hi-0.1
    if y>0:
        y_lo=float(int(y*10))/10
        y_hi=y_lo+0.1
    elif y<0:
        y_hi=float(int(y*10))/10
        y_lo=y_hi-0.1
        # write the results to STDOUT (standard output);
        # what we output here will be the input for the
        # Reduce step, i.e. the input for reducer.py
        # tab-delimited; the trivial word count is 1
    print x_lo,x_hi, y_lo, y_hi, 1


#def upperNum(num):
#	ans = int(num*10)*0.1
#	if (ans < num):
#		ans += 0.1
#	return ans
#def num2Bin(num):
#	upper = upperNum(num)
#	lower = upper - 0.1
#	return "{:.1f},{:.1f}".format(lower, upper)


#import sys

# input comes from STDIN (standard input)
#for line in sys.stdin:
    # remove leading and trailing whitespace
 #   line = line.strip()
    # split the line into words
  #  x,y = line.split()
    # convert x,y (currently a string) to float
  #  x = float(x)
  #  y = float(y)
   # binKey = num2Bin(x) + "," + num2Bin(y)
    #print '%s\t%s' % (binKey, 1)