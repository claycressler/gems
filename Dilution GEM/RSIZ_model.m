function dydt = RSIZ_model(~,y,sigma,mz,a,epsilon,mi,e,ms,r,K)
dydt = zeros(size(y));

% variables
Z = y(1);
I = y(2);
S = y(3);
R = y(4);

dydt(1) = sigma*I - mz*Z - a*Z*(S + I);
dydt(2) = epsilon*a*Z*S - mi*I;
dydt(3) = e*a*R*(S + I) - ms*S - epsilon*a*Z*S;
dydt(4) = r*R*(1 - R/K) - a*R*(S + I);
