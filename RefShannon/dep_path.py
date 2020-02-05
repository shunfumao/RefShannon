from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

extAssembler_paths = {
  "stringtie": "/home/shunfu/software/stringtie-1.3.4d.Linux_x86_64/stringtie",
  "cufflinks": "/home/shunfu/software/cufflinks-2.2.1.Linux_x86_64/cufflinks",
  # need to install boost (./b2 install --prefix=/home/shunfu/local/boost_1_72_0), Bamtools, samtools, 
  "TransComb": "/home/shunfu/software/TransComb_v.1.0/TransComb",
  # gunzip < CLASS-2.1.7.tar.gz | tar -xvf -
  # need to use build.sh for dependency issues
  "CLASS2": "/home/shunfu/software/CLASS-2.1.7/run_class.pl",
  "ryuto": "/home/shunfu/software/ryuto28/ryuto",
  "strawberry" : "/home/shunfu/software/strawberry1.0.2/strawberry",
  "scallop": "/home/shunfu/software/scallop-0.10.2_linux_x86_64/scallop",
  # need java 1.8, 
  "trinity": "Trinity", # export in bashrc: /home/shunfu/software/trinity/trinityrnaseq-v2.9.1/
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
  # trinity depends on bowtie2, jellyfish and salmon
  "bowtie2": "bowtie2", # export in bashrc: /home/shunfu/software/bowtie2/bowtie2-2.3.5.1-linux-x86_64/
  "bowtie2-build": "bowtie2-build",
  "jellyfish": "jellyfish", # export in bashrc: /home/shunfu/software/jellyfish/bin/
  "salmon": "salmon", # export in bashrc: /home/shunfu/software/salmon/salmon-latest_linux_x86_64/bin/
  # TOADD
  "tophat2": "tophat2",
  "hisat2": "hisat2", # export in bashrc: /home/shunfu/software/hisat2/hisat2-2.1.0/
  "hisat2-build": "hisat2-build",
}