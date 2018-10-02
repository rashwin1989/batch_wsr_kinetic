function [dm] = calc_chem_RR_Wiehe3(V,m,k,dt)

% chemical kinetic parameters
k1=abs(k(1)); % MA -> maMA+ + (1-ma-g)D + gG
              % D -> ma1MA+ + (1-ma1-g1)Dnr + g1G
k2=abs(k(2)); % AS -> aAS+ + bMA+ + (1-a-b-g)D + gG
k3=abs(k(3)); % MA+ -> AS+

ma=abs(k(4)); % ma 
a=abs(k(5));  % a
b=abs(k(6));  % b
g=abs(k(7));  % g
ma1=abs(k(12)); %ma1
g1=abs(k(13)); %g1

N = numel(m);
dm = zeros(N,1); 

Vdt = V*dt;

C = m/V;

Cg=C(1,1);
Cd=C(2,1);
Cmap=C(4,1);
Cmal=C(5,1);
Cmah=C(6,1);
Casp=C(7,1);
Cas=C(8,1);
CH2O=C(9,1);

d1_mad = k1*Cd*V*dt;
d1_mal = k1*Cmal*V*dt;
d1_mah = k1*Cmah*V*dt;
d2 = k2*Cas*V*dt;
d3 = k3*Cmap*V*dt;

dm(1,1) = g*(d1_mal + d1_mah) + g*d2 + g1*d1_mad;
dm(2,1) = (1-ma-g)*(d1_mal + d1_mah) + (1-a-b-g)*d2 - d1_mad;
dm(3,1) = (1-ma1-g1)*d1_mad;
dm(4,1) = ma*(d1_mal + d1_mah) + b*d2 - d3 + ma1*d1_mad;
dm(5,1) = -d1_mal;   
dm(6,1) = -d1_mah;
dm(7,1) = a*d2 + d3;
dm(8,1) = -d2;
dm(9,1) = 0;

