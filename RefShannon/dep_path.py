extAssembler_paths = {
  "stringtie": "/home/shunfu/software/stringtie-1.3.4d.Linux_x86_64/stringtie",
  "cufflinks": "/home/shunfu/software/cufflinks-2.2.1.Linux_x86_64/cufflinks",
  "TransComb": "/home/shunfu/software/TransComb_v.1.0/TransComb",
  # gunzip < CLASS-2.1.7.tar.gz | tar -xvf -
  "CLASS2": "/home/shunfu/software/CLASS-2.1.7/run_class.pl",
  "ryuto": "/home/shunfu/software/ryuto28/ryuto",
  "strawberry" : "/home/shunfu/software/strawberry1.0.2/strawberry",
  "scallop": "/home/shunfu/software/scallop-0.10.2_linux_x86_64/scallop",
}

"""
samtools:
  - tar xjf samtools-1.7.tar.bz2
  - ./configure --without-curses --disable-bz2 --disable-lzma
  - make
  - make install
"""
tool_paths = {  
  "samtools": "samtools",
  "blat": "/home/shunfu/software/blat/blat/blat",
  "gffread": "/home/shunfu/software/gffread-0.9.12.Linux_x86_64/gffread",
  "gffcompare": "/home/shunfu/software/gffcompare-0.10.4.Linux_x86_64/gffcompare",
}