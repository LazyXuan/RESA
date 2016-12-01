###
# Copyright 2016 Miler T. Lee, University of Pittburgh
# This file is part of the RESA Suite
#
# RESA Suite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# RESA Suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RESA Suite.  If not, see <http://www.gnu.org/licenses/>.
#
#
# bam_seq_convert.py:
# Given a bamfile that has read IDs of the form <orig_read_id>#<orig_seq>
# replaces the bamfile sequence with the one encoded in orig_seq,
# preserving strand of the bam entry
#
##

import sys
import subprocess
from resa_util import *

def bam_seq_sub(bamfile):
    """
    Iterate through the bamfile, which should have read IDs
    that contain the original read sequence.
    
    This can also be done using a python module called pysam.
    """
    
    cmd = ['samtools', 'view', '-h', bamfile]
    p = subprocess.Popen(cmd, bufsize = 1, stdout = subprocess.PIPE)

    for line in p.stdout:
        line = line.decode('utf-8')
        if line.startswith('@'):
            print(line.strip())
            continue
            
        fields = line.strip().split()
        #In a BAM file, the sequence ID is the first column and
        #the sequence is in the tenth column. The samflag is the
        #second column
        id_fields = fields[0].split('#')
        if len(id_fields) > 1:
            new_id = '#'.join(id_fields[:-1])
            new_seq = id_fields[-1]

            if bam_entry_is_reverse(int(fields[1])):
                new_seq = rc(new_seq)
            
            fields[0] = new_id
            fields[9] = new_seq.upper()
        
        print('\t'.join(fields))

        
if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('Usage: python bam_seq_convert.py <bamfile>\n')
        print('Output is uncompressed SAM. To obtain a bamfile:')
        print('python bam_seq_convert.py <bamfile> | samtools view -hSb -o <out.bam> -\n')
        sys.exit()
        
    bam_seq_sub(sys.argv[1])

