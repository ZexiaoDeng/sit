% x =
%    2.307810312375409
%    0.941286896105092
%    2.000000000000000
% 
% fval =
%   10.918853812662503

clc ;
clear ;
close all ;

p = [ 0 ; 0 ; 0 ; ] ;
q = [ 3 ; 4 ; 2 ; ] ;

M.p = p ;
M.q = q ;

GDC.M = M ;

[ x, fval, output ] = sit_solver03( GDC )

g1_cst( x ) - g2_cst( x )


