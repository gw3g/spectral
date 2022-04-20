
# filename:
diag(d,k,s,mn,mu) = "data/diag.".d."{k=".k."}.(".s.").".mn.".{mu=".mu."}.dat"

# inv. mass
K2(k0,k)=k0*k0-k*k

# only large-frequency
cut(k0,kmax) = k0>kmax ? k0 : NaN
