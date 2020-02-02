# NTRU Encryption Schemes

NTRU is a lattice-based public key encryption scheme that has been
proposed as a quantum secure scheme to replace RSA. 

This directory contains simple proof of concept (which is to say NOT!!!
secure) implementations of several NTRU-adjacent schemes, including:

* The NTRU scheme as originally published, which is vulnerable to certain
  attacks. 

* The NTRU MLS scheme, which addresses some security concerns in the original
  NTRU scheme, and is the scheme currently under consideration by the NIST for
  its post quantum cryptography standard.

* The pqNTRUSign scheme, which uses the same underlying mathematics as NTRU to
  provide a cryptographic signature, similar to RSA. This was eliminated in 
  round 1 of the NIST standardization competition. 

These implementations were produced as part of an independent study I did on 
quantum secure cryptography. 

