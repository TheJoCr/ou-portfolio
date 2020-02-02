# First, define the set of global params
N = 11
q = 23
d = 3

# Set up the rings 
# Quotient ring over the integers
R = ZZ[x].quotient( x^N -1 );
# Quotient ring over a finite base field
Rq = GF(q)[x].quotient( x^N -1 )
# Quotient real ring over the reals
Rr = RR[x].quotient( x^N - 1 );

x = ZZ[x].gen()

# Now, define a few helper methods.

from random import shuffle
# To make a random polynomial for f and g. 
# This is a very bad way to do it (not even using a crypto random number gen)
# but it's easy so I'm going to stick with it
def gen_poly( degree, n ):
    l = [0] * (degree - 2*n - 1)
    for i in range (n):
        l.append(1)
        l.append(-1)
    l.append(1);
    shuffle( l )
    return R(l)


def centerlift( poly, m ):
    coeff = poly.list()
    lifted_coeff = map( lambda a: Mod(a, m).lift_centered(), coeff )
    return R(lifted_coeff)

def r( poly ):
    return R( map( round, poly.list() ))

# Returns a pair of keys (private key, verification key)
def gen_key():
    # Key gen is the most complicated part of this process
    f = gen_poly(N, d)
    g = gen_poly(N, d)
    
    f_inv = Rq(f) ^ -1
    # This computed h is going to be the verification key
    h = centerlift(f_inv * g, q)

    # Now we have to find the complementary half basis
    Rf, f1, f2 = f.lift().xgcd(x^N - 1)
    Rg, g1, g2 = g.lift().xgcd(x^N - 1)
    gcdRfRg, Sf, Sg = Rf.xgcd(Rg)

    A = q * Sf * f1
    B = q * Sg * g1

    fi = Rr(f) ^ -1
    gi = Rr(g) ^ -1

    C = .5 * (A * gi - B * fi)
    # Round to the nearest integer
    C = r(C)

    # This is the computed basis - our private key then consists of
    # the basis (f, F) (g, G)
    F = -B - C * f
    G = A - C * g
    return ((f, F, g, G), h)

# Computes the signature of D=(D1,D2) given the private key
def sign( D1, D2, priv_key ):
    f, F, g, G = priv_key
    v1 = r(( D1 * G - D2 * F ) / q)
    v2 = r((-D1 * g + D2 * f) / q)
    s = v1 * f + v2 * F
    return s

def verify( D1, D2, s, h ):
    # And the verification of that signature is then as simple as:
    t = Rq(h * s)
    # Centerlifting the difference between D2 and t gives the smallest
    # possible distance between the two. 
    t_diff = centerlift( t - r(D2), q)
    s_diff = s - r(D1)
    # And compute the actual distance:
    t_sos = sum(map(lambda a: a**2, t_diff))
    s_sos = sum(map(lambda a: a**2, s_diff))
    # Then the total distance is:
    return RR(sqrt( s_sos + t_sos ))

## Example program:
# Samantha generates a key set:
priv, verf = gen_key()
print( "Samantha computes a good basis and publishes h = " + str(verf.list()))

# She then computes the signature to the document:
D1 = Rr([4,7,10,-6,11,11,6,-6,9,-9,3])
D2 = Rr([-7,11,11,1,-3,-7,-2,0,-3,9,5])

s = sign( D1, D2, priv )
print( "Samantha signs with " + str( s.list() ) )

# Victor then verifies that this is sufficently good:
print( "Victor verifies a distance of: " + str(verify( D1, D2, s, verf )) )


