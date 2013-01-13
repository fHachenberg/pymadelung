#!/usr/bin/python
# -*- coding: utf-8 -*-

# by Fabian Hachenberg, 2010

#this script calculates the Madelung constant for a supercell of size (x,y,z)
#the Madelung constant is normed to the identity distance L = V^(1/3)

import numpy as np
import math
import madelung

modules ={}

def main():
    fl = ["NaCl_POSCAR"]
    for src in fl:
        print src
        m,L = madelungfromposcar(src)
        print "m=", m
	print "L=", L
    
def madelungfromposcar( poscarfn ):
    try:
        f = open(poscarfn)
        
        comment = f.readline()
        alpha = float(f.readline())
        a = tuple( float(v) for v in f.readline().split() if v != "" )
        b = tuple( float(v) for v in f.readline().split() if v != "" )
        c = tuple( float(v) for v in f.readline().split() if v != "" )
        #print alpha,a,b,c
        return sc_mad( (a,b,c), alpha)
    except IOError:
        print "Error reading the POSCAR", poscarfn
        return
    
def sc_mad( aorig, alpha, (x,y,z)=(1,1,1)):    
    
    V = alpha**3 * np.dot(np.cross(x*aorig[0],y*aorig[1]), z*aorig[2])
    L = math.pow(V,0.3333333333)
    a = tuple((x*o[0]/L, y*o[1]/L, z*o[2]/L) for o in aorig)
    n = 10

    m = madelung.madelung(a, (0,0,0), n, 0, False)
    return (m,L)

if __name__ == '__main__':
    main()
