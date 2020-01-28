import os,pdb, numpy
from filter_FP_batch import calcFP, calcSens

def run_cmd(cmd):
    print(cmd)
    os.system(cmd)




def waterfall(file_list, choice):
    pdb.set_trace()
    thresh_list = [0.9,0.92,0.94,0.96,0.98]
    for filename in file_list:
        rec = [0]*len(thresh_list)

        if choice==0: #sens
            for (i,thresh) in enumerate(thresh_list):
                rec[i] = calcSens(filename, mtl=200, recThreshold=thresh)[0]

            print('\n'+filename)
            for idx in range(len(thresh_list)):
                print('%f\t%d'%(thresh_list[idx], rec[idx]))

        elif choice==1: #false positive
            for (i,thresh) in enumerate(thresh_list):
                rec[i] = calcFP(filename, mtl=200, tr_labels=None, fpThreshold=thresh)[1]

            print('\n'+filename)
            for idx in range(len(thresh_list)):
                print('%f\t%f'%(thresh_list[idx], rec[idx]))

        else:
            print('unexpected choice')
            pdb.set_trace()

if __name__ == '__main__':
    import sys
    args = sys.argv
    if '-0' in args:
        choice = 0 # 'sens'
    elif '-1' in args:
        choice = 1 # 'false positive'
    else:
        print('unexpected choice')
        pdb.set_trace()

    f = sys.argv[2:]
    waterfall(f, choice)




