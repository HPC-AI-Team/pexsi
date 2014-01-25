function z = DosInertia( x, y )
%% DOSINERTIA computes the density of states (DOS) from inertia count.
%
%  z = DosInertia(x,y) returns the density of states (DOS) from the
%  inertia count, or cumulative density of states (CDOS) y given at a
%  sequence of points x.  The computation is done via polynomial
%  interpolation (pchip) and differentiation.

pp = pchip( x, y );
dpp = fnder( pp );
z = ppval( dpp, x );

