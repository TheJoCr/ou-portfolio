
#include <iostream>
#include <random>
using namespace std;

typedef unsigned long long u64;

class RSA {
private:
	u64 p;
	u64 q;
	u64 n;
	u64 d;
	u64 e;
public:
	RSA();
	u64 getN();
	u64 getE();
	u64 decrypt(u64 m);
};

class random {
private:
	mt19937 gen;
	uniform_int_distribution<u64> dis;
public:
	random(u64 min, u64 max);
	u64 rand();
};

class mathUtils {
private:
	const static int tocheck[]; 
public:
	static u64 modexp(u64 base, u64 power, u64 mod);
	static u64 modmult(u64 a, u64 b, u64 mod);
	static u64 modinv(u64 a, u64 mod);
	static u64 randomPrime(u64 min, u64 max);
	static bool isPrime(u64 p, random randomGen);
};


u64 RSAEncode(u64 m, RSA r)
{
	return mathUtils::modexp(m, r.getE(), r.getN());
}

void main() {
	for (int i = 0; i < 100; i++) {
		u64 message = 696969;
		RSA rsa = RSA();
		u64 e = RSAEncode(message, rsa);
		cout << "Encrypted message:" << e << endl;
		cout << "Decrypted message:" << rsa.decrypt(e) << endl;
		cout << "Yay it worked! or not..." << endl;
	}
	cout << "over";
	return;
}


RSA::RSA() 
{
	p = mathUtils::randomPrime(0x000000000000ffffu, 0x00000000ffffffffu); //32 bit prime
	q = mathUtils::randomPrime(0x000000000000ffffu, 0x00000000ffffffffu); //32 bit prime
	n = p*q; //64 bit max b/c the primes are small;
	d = mathUtils::randomPrime(1, (p - 1)*(q - 1));
	e = mathUtils::modinv(d, (p - 1)*(q - 1));

	cout << "d: " << d << ", e: " << e << ", d*e mod(p-1)(q-1)" << mathUtils::modmult(d, e, (p - 1)*(q - 1)) << endl;
}

u64 RSA::getN()
{
	return n;
}

u64 RSA::getE()
{
	return e;
}

u64 RSA::decrypt(u64 m)
{
	return mathUtils::modexp(m, d, n);
}

u64 mathUtils::modexp(u64 base, u64 power, u64 mod)
{
	u64 out = 1;
	base %= mod;
	while (power > 0) {
		out %= mod;
		if (power & 1) {
			out = modmult(out, base, mod);
		}
		base = modmult(base, base, mod);
		power >>= 1;
	}

	return out;
}

u64 mathUtils::modmult(u64 a, u64 b, u64 mod)
{
	u64 res = 0;
	u64 temp;

	if (b >= mod) b %= mod;
	
	while (a != 0) {
		if (a & 1) {
			if (b >= mod - res)
				res -= mod;
			res += b;
		}
		a >>= 1;

		temp = b;
		if (b >= mod - b) 
			temp -= mod;
		b += temp;
	}
	return res;
}

u64 mathUtils::modinv(u64 a, u64 mod)
{
	u64 t = 0;
	u64 newT = 1;
	u64 r = mod;
	u64 newR = a % mod;

	u64 temp;

	while (newR != 0) {
		temp = newT;
		newT = t - (r / newR) * newT;
		t = temp;

		temp = newR;
		newR = r - (r / newR) * newR;
		r = temp;
	}
	return t;
}

u64 mathUtils::randomPrime(u64 min, u64 max) 
{
	u64 out;
	random r = random(min, max);
	do {
		out = r.rand();
	} while (!isPrime(out, r));
	return out;
}

const int mathUtils::tocheck[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

bool mathUtils::isPrime(u64 p, random randomGen) 
{
	if (p % 2 == 0)
		return false; 
	u64 d = p - 1;
	u64 r = 0;
	while (d % 2 == 0) {
		r++;
		d >>= 1;
	}
	for (int i = 0; i < 12; i++) {
		u64 a = tocheck[i];
		u64 x = modexp(a, d, p);
		if (x == 1 || x == p - 1)
			continue;
		for (int j = 0; j < r - 1; j++) {
			x = modexp(x, 2, p);
			if (x == 1)
				return false;
			if (x == p - 1)
				break;
		}
		if (x == p - 1)
			continue;
		return false;
	}
	return true;
}

random::random(u64 min, u64 max) 
{
	random_device rd;
	gen = mt19937(rd());
	dis = uniform_int_distribution<u64>(min, max);
}

u64 random::rand() 
{
	return dis(gen);
}