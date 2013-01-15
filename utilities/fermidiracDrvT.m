function g = fermidiracDrvT(ev,efermi,Tbeta);
% The derivative of the Fermi-Dirac function with respect to 1/beta.
dif = ev - efermi;
if( numel(efermi) > 1 )
	error('efermi has to be a scalar');
end
g = fermidirac( ev, efermi, Tbeta ) .* Tbeta.^2 .* dif;
ind1 = find(dif > 0 );
ind2 = find(dif <= 0 );
g(ind1) = g(ind1) ./ ( 1 + exp( -Tbeta * dif(ind1) ) );
g(ind2) = g(ind2) .* exp( Tbeta * dif(ind2) ) ./ ...
	( 1 + exp( Tbeta * dif(ind2) ) );
