#!/usr/bin/env python

from operator import itemgetter
import sys

current_x_lo = None
current_x_hi = None
current_y_lo = None
current_y_hi = None
current_count = 0


# input comes from STDIN and it is the output of mapper.py
for line in sys.stdin:
        # remove leading and trailing whitespace
        line = line.strip()

        # parse the input we got from mapper.py
        x_lo,x_hi, y_lo, y_hi,count = line.split('\t')

          # convert count (currently a string) to int
        try:
             count = int(count)
        except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
            continue

        # this IF-switch only works because Hadoop sorts map output
        # by key (here: word) before it is passed to the reducer
        if current_x_lo == x_lo and current_x_hi == x_hi and current_y_lo == y_lo and current_y_hi == y_hi:
                current_count += count
        else:
                if current_x_lo and current_x_hi and current_y_lo and current_y_hi:
                        # write result to STDOUT
                        print current_x_lo ， current_x_hi， current_y_lo， current_y_hi，current_count
                current_count = count
                current_x_lo = x_lo
                current_x_hi = x_hi
                current_y_lo = y_lo
                current_y_hi = y_hi

# do not forget to output the last word if needed!
if current_x_lo == x_lo and current_x_hi == x_hi and current_y_lo == y_lo and current_y_hi == y_hi:
        print current_x_lo ， current_x_hi， current_y_lo， current_y_hi，current_count







#!/usr/bin/env python

from operator import itemgetter
import sys

current_binKey = None
current_count = 0
binKey = None

# input comes from STDIN
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()

    # parse the input we got from mapper.py
    binKey, count = line.split('\t', 1)

    # convert count (currently a string) to int
    try:
        count = int(count)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        continue

    # this IF-switch only works because Hadoop sorts map output
    # by key (here: binKey) before it is passed to the reducer
    if current_binKey == binKey:
        current_count += count
    else:
        if current_binKey:
            # write result to STDOUT
            print '%s\t%s' % (current_binKey, current_count)
        current_count = count
        current_binKey = binKey

# do not forget to output the last word if needed!
if current_binKey == binKey:
    print '%s\t%s' % (current_binKey, current_count)


