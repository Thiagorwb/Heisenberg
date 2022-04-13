import numpy as np

def gray_code(num):
    shift = (num >> 1)
    return (shift ^ num)

def gray_index(num):
    
    return int(np.log2((gray_code(num) ^ gray_code(num + 1))))

def setBitNumber(n):
    if (n == 0):
        return 0;
 
    msb = 0;
    n = int(n / 2);
 
    while (n > 0):
        n = int(n / 2);
        msb += 1;
 
    return (1 << msb);

def gray(i):
    return int(np.log2(setBitNumber(GrayToBinary(gray_index(i)))))

def GrayToBinary(gray_num):

    mask = gray_num
    while (mask):           
        mask >>= 1;
        gray_num ^= mask;
    
    return gray_num

def binary(num, num_bits):
    # returns array of bools corresponding to num
    mask = num
    bins = []
    while mask > 0:
        rem = mask % 2
        bins.append(rem)
        mask = mask//2
        
    nbits = len(bins)
    bins = np.array(bins)
    if nbits <= num_bits:
        pad = np.zeros((num_bits - nbits))
        bins = np.concatenate((bins, pad))
    else:
        bins = bins[:num_bits]

    return bins