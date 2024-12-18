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

rep.l = p ;
rep.u = q ;

P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep   
CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep

plot( P ) ;
hold on


GDC.M = M ;

[ x, fval, output ] = sit_solver04( GDC )


