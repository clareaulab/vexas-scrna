import pysam
import sys

bamfile = sys.argv[1]

samfile = pysam.Samfile(bamfile, "rb")
for pileupcolumn in samfile.pileup('chrX', 20534934, 20534934+1):
    if(pileupcolumn.pos == 20534934):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ex = pileupread.alignment
                try:
                    print(base + "," + ex.get_tag("CB"))
                except:
                    pass
samfile.close()
