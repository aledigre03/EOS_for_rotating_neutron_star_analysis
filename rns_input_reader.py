import re

with open('output_kepler.txt', 'r') as file:
    data = file.read()

#the actual script for findall in the '' was found by an AI
m_values = re.findall(r'\n\s*([\d\.e\+-]+)\s+M\s+\(M_sun\)', data)
re_values = re.findall(r'\n\s*([\d\.e\+-]+)\s+R_e\s+\(km\)', data)
rp_re_values = re.findall(r'\n\s*([\d\.e\+-]+)\s+r_p/r_e', data)
ec_values = re.findall(r'\n\s*([\d\.e\+-]+)\s+e_c\s+\(10\^15 gr/cm\^3\)', data)
omega_values = re.findall(r'\n\s*([\d\.e\+-]+)\s+Omega\s+\(10\^4 s\^-1\)', data)

M = [float(m) for m in m_values]
R_e = [float(re) for re in re_values]
r_p_over_r_e = [float(rp_re) for rp_re in rp_re_values]
e_c = [float(ec) for ec in ec_values]
Omega = [float(omega) for omega in omega_values]

assert len(M) == 30, "Expected 30 M values"
assert len(R_e) == 30, "Expected 30 R_e values"
assert len(r_p_over_r_e) == 30, "Expected 30 r_p/r_e values"
assert len(e_c) == 30, "Expected 30 e_c values"
assert len(Omega) == 30, "Expected 30 Omega values"

print("M =", M)
print("R_e =", R_e)
print("r_p_over_r_e =", r_p_over_r_e)
print("e_c =", e_c)
print("Omega =", Omega)

