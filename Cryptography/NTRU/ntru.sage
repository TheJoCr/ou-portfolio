# First, define the set of global params
p = 3
N = 7
q = 41
d = 2

# Set up the rings 
R.<x> = ZZ[]
R = R.quotient( x^N -1 );
Rp = GF(p)[x].quotient( x^N -1 )
Rq = GF(q)[x].quotient( x^N -1 )

x = R.gen()

# Now, define a few helper methods.

from random import shuffle
# To make a random polynomial for f and g. 
# This is a very bad way to do it, but it's easy
def gen_g( degree, n ):
    l = [0] * (degree - 2*n)
    for i in range (n):
        l.append(1)
        l.append(-1)
    shuffle( l )
    print(l)
    return R(l)


def gen_f( degree, n ):
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

# Generates a (private key, public key pair)
def gen_key():
    f = gen_f(N, d)
    g = gen_g(N, d)

    Fp = Rp(f) ^ -1
    Fq = Rq(f) ^ -1
    
    h = Fq * g
    return (f, h)

# Then h is the public key

# Then we need to encrypt some data
# The message and an ephemeral key are generated
# this returns the encrypted message given a public key
def encrypt(m, pub):
    r = gen_f(N, d)
    e = Rq(p * pub * r + m)
    return e

def decrypt(e, priv):
    # Since we are lazy with storage, recompute Fp
    f = priv
    Fp = Rp(f) ^ -1
    # And do the actual computation
    # This is done with two centerlifts:
    a = centerlift( Rq( f * e ), q )
    clear_text = centerlift( Rp( Fp * a ), p )
    return clear_text


# Little test program:
(priv, pub) = gen_key()
print("Private key: " + str(priv))
print("Public key: " + str(pub))
# Bob sends Alice the message:
m_bob = -x^5 + x^3 + x^2 - x +1
# Which encrypts to 
e = encrypt( m_bob, pub )
print("Encrypted message: " + str(e))
# Alice decrypts
m_alice = decrypt( e, priv )

print("Bobs message: " + str(m_bob));
print("Alice message: " + str(m_alice));
