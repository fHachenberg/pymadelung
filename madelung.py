#!/usr/bin/python
# -*- coding: utf-8 -*-

#madelung.py, version 1.0 by Fabian Hachenberg, 2010
#Calculates the Madelung constant for an infinite lattice of point charges, neuralized by a counter-charged jellium background
#Requirements: SciPy, Python 2.6x
#TODO: atm, the vector routines are hardcoded. Using NumPy woud speed things up!

import sys
import math
import itertools #we need this for the creation of the lattice sums
import scipy.special #we need this for the erfc function!
from optparse import OptionParser

#Vector math. Can be replaced by numpy routines for example!
def crossproduct(a,b):
	return ((a[1]*b[2]-a[2]*b[1]),(a[2]*b[0]-a[0]*b[2]),(a[0]*b[1]-a[1]*b[0]))

def product(a,s):
	return (a[0]*s,a[1]*s,a[2]*s)

def dotproduct(a,b):
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(a):
	return math.sqrt(a[0]**2+a[1]**2+a[2]**2)

def transform(x,A):
	return (A[0][0]*x[0]+A[1][0]*x[1]+A[2][0]*x[2], A[0][1]*x[0]+A[1][1]*x[1]+A[2][1]*x[2], A[0][2]*x[0]+A[1][2]*x[1]+A[2][2]*x[2])

def minus(a,b):
	return tuple( x-y for x,y in itertools.izip(a,b) )

#does the actual calculations. input:
# -a an array of 3 lattice vectors in real space. They have to be specified in units of the Madelung constant identity distance. For example if the identity distance is the lattice constant in a cubic lattice, the vectors have to have a length of 1. A common approach compatible for all lattices is to use the third root of the volume of the unit cell V^(1/3) as the identity distance. Then the volume of the lattice vectors specified here has to be 1
# -r the offset in relative coordinates (relative to a). Therefore they have to be contained in a cube of length 1, r is element of (0,1)^3 if the self term is not omitted!
# -iterations Number of iterations, determining the size of the biggest cube to use for summation. For NaCl, convergence up to all digits is obtained already for less than 10 iterations
# -eta Used for numerical purposes, the result of a _CONVERGED_ summation has to be independent of the value of eta, but eta can determine how quickly the sum converges. None or 0 trigger an automatic choice of eta, based on a statement in L Martin "Electronic Structure"
# -selfterm include the selfterm? This is the T=0 term in the lattice sum "sum_T(1/norm(T-r))". You should only omit the self term, if you're using r = 0 
def madelung( a, r, iterations, eta = None, selfterm = True):

	V = dotproduct(crossproduct(a[0],a[1]),a[2]) #calculate the volume of the cell specified by a
	b = (product(crossproduct(a[1],a[2]),2.0*math.pi/V), product(crossproduct(a[2],a[0]),2.0*math.pi/V), product(crossproduct(a[0],a[1]),2.0*math.pi/V)) #calculate the reciprocal lattice vectors
	x,y,z = transform(r,a) #transform r into cartesian coordinates

	if not eta or eta == 0:
		eta = min([norm(v) for v in b]) #l martin writes in "electronic structure", best convergence is obtained for eta = <formula> using the SMALLEST reciprocal lattice vector!

	if selfterm:
		madelung = scipy.special.erfc(eta*norm((x,y,z))) / norm((x,y,z)) #this is the zero term of the lattice sum, infinite for r = 0 
	else:
		madelung = - 2.0*eta/math.sqrt(math.pi) #this is a correction for the reciprocal lattice sum, removing a self-interaction term. It's insufficient to simply omit the T = 0 term in the real space lattice sum!
	#the sum. we're summing up over expanding cubes
	for L in range(1,iterations): #L is the size of the cube
		#ab contains x,y-pairs for lattice points on the xz,yz-face of the cube. first the xz plane, then the -xz plane, then the yz plane, then the -yz plane
		ab =  itertools.chain( itertools.izip(itertools.repeat(L,2*L+1),range(-L,L+1)),itertools.izip(itertools.repeat(-L,2*L+1),range(-L,L+1)), itertools.izip(range(-L+1,L), itertools.repeat(L,2*L-1)), itertools.izip(range(-L+1,L), itertools.repeat(-L,2*L-1)))
		#in the next step, we're uniting the set ab containing only 2d coordinates with the z coordinate. now we're creating planes for real. first we're mutliplying the ab set with the range of z coordinates. this creates the "mantle" of the cube. next we're creating the top plane of the cube, after this we're creating the bottom plane
		abc = itertools.chain( itertools.product(ab,range(-L,L+1)), itertools.izip(itertools.product(range(-L+1,L),range(-L+1,L)), itertools.repeat(L,(2*L-1)**2)), itertools.izip(itertools.product(range(-L+1,L),range(-L+1,L)),itertools.repeat(-L,(2*L-1)**2)) )
		for (u,v),w in abc:#we have to use this strange notations (u,v),w ; now we're going through all lattice points on the surface of the current cube
			G = transform((u,v,w),b) #the reciprocal lattice vector obtained by using the u,v,w integer coordinates
			T = minus(transform((u,v,w),a),(x,y,z)) #the real space lattice vector obtained by using the u,v,w integer coordinates
			nT = norm(T)
			nG = norm(G)
			#this is the ewald sum - the combination of a real space sum (first term) and a reciprocal space sum (second term)
			madelung+=scipy.special.erfc(eta*nT) / nT
			madelung+= 4*math.pi/V * 1.0/(nG**2) * math.exp(-(nG**2)/(4*eta**2)) * math.cos(G[0]*x+G[1]*y+G[2]*z)
	madelung += - math.pi/(eta**2*V) #important. the independence on eta is finally obtained by adding this term
	return madelung

#test function calculating the Madelung constant for a NaCl ionic crystal, using the NN-distance as identity distance
def NaClTest():
	ast = (( 2, 0, 0), (0, 2, 0), (0, 0, 2)) #because in NaCl the next equivalent atom is always the atom next to the next neighbour, we have to use a 2.0 
	n = 10 #number of iterations, result should even be exact up to all digits for less than 10 iterations
	print "Calculating Madelung Constant for NaCl using", n, "iterations"
	#we have to put the actual Madelung constant together from the contributions of the individual non-equivalent ions
	#the first term is the interaction of the ion with it's own images, the second terms are the interaction with the 3 other ions of the same species within the cell, sitting on the faces of the cell. The third term is the interaction with the 3 ions of opposite species sitting on the edges of the cell and finally, the last term, is the interaction with the ion of opposite species sitting in the center of the cell. 
	m = madelung( ast, (0,0,0), n, 0, False) + 3.0 * madelung( ast, (0.5,0.5,0), n, 0, True) - 3.0 * madelung( ast, (0.5,0,0), n, 0, True) - madelung( ast, (0.5,0.5,0.5), n, 0, True)
	print "M =", m

if __name__ == "__test__":
	NaClTest()

if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option("-a", "--a", dest="a", action="store", type="string", help="first lattice vector, format: x,y,z", metavar="A")
	parser.add_option("-b", "--b", dest="b", action="store", type="string", help="second lattice vector, format: x,y,z", metavar="B")
	parser.add_option("-c", "--c", dest="c", action="store", type="string", help="third lattice vector, format: x,y,z", metavar="C")
	parser.add_option("-r", "--r", dest="r", action="store", type="string", help="offset vector in relative lattice coordinates, format: a,b,c", metavar="R")
	parser.add_option("-i", "--iterations", dest="iterations", action="store", default=10, type="int", help="number of iterations", metavar="ITERATIONS")
	parser.add_option("-e", "--eta", dest="eta", action="store", default=0, type="float", help="(optional) specify value for eta - can determine the required number of iterations for convergence", metavar="ETA")
	parser.add_option("-o", "--omitselfterm", dest="omitselfterm", action="store_true", help="(optional) omit the self interaction term")
	(options, args) = parser.parse_args()
	
	if not options.a or not options.b or not options.c:
		print "Error, please specify all lattice vectors a b c!"
		raise SystemExit
	a = ( tuple( float(v) for v in options.a.split(",") if v != "" ), tuple( float(v) for v in options.b.split(",") if v != "" ), tuple( float(v) for v in options.c.split(",") if v != "" ))
	if len(a[0]) != 3 or len(a[1]) != 3 or len(a[2]) != 3:
		print "Error while parsing the provided lattice vectors! Length of one of them is < 3"
		raise SystemExit

	if options.r:
		r = tuple( float(v) for v in options.r.split(",") if v != "" )
	else:
		r = (0,0,0)
		options.omitselfterm = True

	if options.omitselfterm == None:
		selfterm = True
	else:
		selfterm = False
	print madelung( a, r, options.iterations, options.eta, selfterm )
