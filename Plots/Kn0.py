import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("Kn0.out", "r")

k    = []
knr1 = []
kni1 = []
knr2 = []
kni2 = []
knr3 = []
kni3 = []
knr4 = []
kni4 = []
knr5 = []
kni5 = []
knr6 = []
kni6 = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    c6      = float(numbers[5])
    c7      = float(numbers[6])
    c8      = float(numbers[7])
    c9      = float(numbers[8])
    c10     = float(numbers[9])
    c11     = float(numbers[10])
    c12     = float(numbers[11])
    c13     = float(numbers[12])
    k.append(c1)
    knr1.append(c2)
    kni1.append(c3)
    knr2.append(c4)
    kni2.append(c5)
    knr3.append(c6)
    kni3.append(c7)
    knr4.append(c8)
    kni4.append(c9)
    knr5.append(c10)
    kni5.append(c11)
    knr6.append(c12)
    kni6.append(c13)
                      
fig = plt.figure (figsize=(8.0, 8.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.xlim (0., k[-1])

plt.plot    (k, knr1, color='red',    linewidth = 2, linestyle = 'solid')
plt.plot    (k, kni1, color='red',    linewidth = 2, linestyle = 'dashed')
plt.plot    (k, knr2, color='green',  linewidth = 2, linestyle = 'solid')
plt.plot    (k, kni2, color='green',  linewidth = 2, linestyle = 'dashed')
plt.plot    (k, knr3, color='blue',   linewidth = 2, linestyle = 'solid')
plt.plot    (k, kni3, color='blue',   linewidth = 2, linestyle = 'dashed')
plt.plot    (k, knr4, color='yellow', linewidth = 2, linestyle = 'solid')
plt.plot    (k, kni4, color='yellow', linewidth = 2, linestyle = 'dashed')
plt.plot    (k, knr5, color='cyan',   linewidth = 2, linestyle = 'solid')
plt.plot    (k, kni5, color='cyan',   linewidth = 2, linestyle = 'dashed')
plt.plot    (k, knr6, color='orange', linewidth = 2, linestyle = 'solid')
plt.plot    (k, kni6, color='orange', linewidth = 2, linestyle = 'dashed')
plt.axhline (0.,      color='black',  linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{k}$', fontsize="20")
plt.ylabel(r'$K_{n,0}$', fontsize="20")
#plt.legend(fontsize="20")

plt.show ()
