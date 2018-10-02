function [dm] = calc_chem_RR_Wiehe1(V,m,k,dt)

% chemical kinetic parameters
k1=abs(k(1)); % MA -> maMA+ + (1-ma)D
k2=abs(k(2)); % AS -> aAS+ + bMA+ + (1-a-b)D
k3=abs(k(3)); % MA+ -> AS+

ma=abs(k(4)); % ma 
a=abs(k(5));  % a
b=abs(k(6));  % b

N = numel(m);
dm = zeros(N,1); 

Vdt = V*dt;

C = m/V;

Cd=C(1,1);
Cmap=C(2,1);
Cmal=C(3,1);
Cmah=C(4,1);
Casp=C(5,1);
Cas=C(6,1);
CH2O=C(7,1);

d1_mal = k1*Cmal*V*dt;
d1_mah = k1*Cmah*V*dt;
d2 = k2*Cas*V*dt;
d3 = k3*Cmap*V*dt;

dm(1,1) = (1-ma)*(d1_mal + d1_mah) + (1-a-b)*d2;
dm(2,1) = ma*(d1_mal + d1_mah) + b*d2 - d3;
dm(3,1) = -d1_mal;   
dm(4,1) = -d1_mah;
dm(5,1) = a*d2 + d3;
dm(6,1) = -d2;
dm(7,1) = 0;

