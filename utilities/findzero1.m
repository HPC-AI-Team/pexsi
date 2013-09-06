function [xHist, fHist, gHist] = findzero1( fFunc, gFunc, x0, ...
	xmin, xmax, fmin, fmax, gmin, gmax, maxiter, fTol )
% FINDZERO1 uses the derivative information for root finding of a
% monotonic function.
%
% Lin Lin
% 12/28/2012

if( fmin >= 0 || fmax <= 0 )
	error('Problem with the initial interval.');
end

x     = x0;
xHist = [];
fHist = [];
gHist = [];

for iter = 1 : maxiter
	fVal = fFunc(x);
	gVal = gFunc(x);

	xHist = [xHist; x];
	fHist = [fHist; fVal];
	gHist = [gHist; gVal];

	% Shrink the interval
	if( abs( fVal ) < fTol )
		fprintf('Found zero in %d steps, xsol=%25.15f.\n', iter, x );
		break;
	end	
	if( fVal < 0 )
		xmin = x;
		fmin = fVal;
		gmin = gVal;
	else
		xmax = x;
		fmax = fVal;
		gmax = gVal;
	end

	% FIXME MAGIC NUMBER
	if( abs( fVal ) > 1e-10 )

		% Cubic interpolation and root finding
		% coef(1)*(x-xmin)^3 + coef(2)*(x-xmin)^2 + coef(3) * (x-xmin) +
		% coef(4)
		coef(3) = gmin;
		coef(4) = fmin;
		A   = [  (xmax-xmin)^2    (xmax-xmin); 
		3*(xmax-xmin)      2];
		Ainv = - [2                 -(xmax - xmin);
		-3 * (xmax-xmin)   (xmax - xmin)^2] / (xmax-xmin)^3;
		rhs = [ fmax - fmin - gmin * (xmax - xmin);
		gmax - gmin ];

		coef(1:2) = Ainv * rhs;

		rootspoly = roots(coef);

		x = rootspoly( find(rootspoly <= xmax-xmin & rootspoly >= 0) );
		if( isempty(x) )
			xHist
			error('Could not find the root.');
		end
		x = real(x(1) + xmin);

	else

		% Newton's method
		x = x - fVal / gVal;
		
	end

end
