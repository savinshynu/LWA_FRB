from __future__ import print_function

import sys

success = True
is_multi = False
with open(sys.argv[1], 'r') as fh:
    for line in fh:
        if line.find('Proposed File Time Alignment') != -1:
            is_multi = True
            print("a 'multi' run")
        if '\r' in line:
            parts = line.split('\r')
            if parts[-1].find('[===') != -1:
                line = parts[-1]
            else:
                line = parts[-2]
                 
        if line.find('[===') != -1:
            print(line.strip().rstrip())
            #print(line)
            fields = line.split()
            #print(fields)
            percent = float(fields[-2].strip('%'))
            if percent < 99.0:
                success &= False
                
if not success:
    sys.exit(1)
