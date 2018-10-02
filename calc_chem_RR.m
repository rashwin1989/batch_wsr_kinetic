function [dm] = calc_chem_RR(V,m,k,dt)

% chemical kinetic parameters
k1=abs(k(1)); % AS -> C
k2=abs(k(2)); % AS -> G
k3=abs(k(3)); % MAvr -> AS; MAd -> AS
k4=abs(k(4)); % MAvr -> MAd+; MAvr -> MAd 
k5=abs(k(5)); % MAvr -> G 

a1=abs(k(7));
a2=abs(k(8));
b1=abs(k(9)); 
b2=abs(k(10));

N = numel(m);
dm = zeros(N,1); 

Vdt = V*dt;

C = m/V;

Cma_dr = C(1,1);
Cma_vr = C(2,1);
Cas = C(3,1);
Cw = C(4,1);
Cma_dnr = C(5,1);
Cg = C(6,1);
Cc = C(7,1);

das = (Cas^a1)*(Cc^b1)*Vdt;
dma_vr = (Cma_vr^a2)*(Cc^b2)*Vdt;   
dma_dr = (Cma_dr^a2)*(Cc^b2)*Vdt;

dm(1,1) = k4*dma_vr - (k4+k3+k5)*dma_dr;
dm(2,1) = -(2*k4+k3+k5)*dma_vr;
dm(3,1) = k3*(dma_vr+dma_dr) - (k1+k2)*das;   
dm(4,1) = 0;
dm(5,1) = k4*(dma_vr+dma_dr);
dm(6,1) = k2*das + k5*(dma_vr+dma_dr);
dm(7,1) = k1*das;

