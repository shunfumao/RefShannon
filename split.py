import sys, os, pdb

'''
hits.sam --> chr1/hits.sam, chr2/hits.sam ...

usage:

python split.py -O chrs_dir -i sam_file

'''

def parse_args_split():

    args = sys.argv

    args_error = False 

    if '-O' in args:
        chrs_dir = args[args.index('-O')+1]
    else:
        args_error = True

    if '-i' in args:
        sam_file = args[args.index('-i')+1]
    else:
        args_error = True

    return [args_error, sam_file, chrs_dir]

def split():

    [args_error, filename, loc]=parse_args_split()
    if args_error==True:
        print('split arguments error')
        return

    out_files = {}
    with open(filename) as f:
        for line in f:
            if line[0] == '@': continue
            fields = line.split()
            flags, genome = int(fields[1]), fields[2]
            if (flags >> 2) & 1:
                genome = 'chrUnmapped'
            if genome not in out_files:
                os.mkdir(loc+genome)
                out_files[genome] = open("{}/hits.sam".format(loc+genome), 'w', 0)
            out_files[genome].write(line)

    for _, out_file in out_files.items():
        out_file.close()

    return

if __name__ == '__main__':
    split()
