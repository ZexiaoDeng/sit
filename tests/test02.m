% ====================
%   min     f(x) = x1^2 + x2^2 - x3^2 + 8
%   s.t.    x1^2 - x2 + x3^2 >= 0
%           x1 + x2^2 - x3^2 <= 20
%           -x1 - x2^2 + 2 = 0
%           x2 + 2*x3^2 = 3
%           x1, x2, x3 >= 0
%

clc ;
clear ;
close all ;

Aineq = [] ;
bineq = [] ;
Aeq   = [] ;
beq   = [] ;
lb    = zeros( 3, 1 ) ;
ub    = [] ;


options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','sqp-legacy', ...
                        'PlotFcns' ,'optimplotfval' ) ;
x0 = zeros( 3, 1 ) ;

% x0 = rand( 3, 1 ) ;

[ x       , ...
  fval    , ...
  exitflag, ...
  output   ] = fmincon( @oracle ,        ...
                        x0      ,        ...
                        Aineq   , bineq, ...
                        Aeq     , beq  , ...
                        lb      , ub   , ...
                        @nonlcon,        ...
                        options         )

[ c, ceq ] = nonlcon( x )

% ====================
% 目标函数
% ====================
function f = oracle( x )
    % 目标函数
    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    f = x1^2 + x2^2 - x3^2 + 8 ;
    
end

% ====================
% 非线性约束函数
% ====================
function [ c, ceq ] = nonlcon( x )
    % 非线性等式和不等式约束函数
    % c(x) <= 0
    % ceq( x ) = 0
    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    
    c1 =  -( x1^2 - x2 + x3^2 ) ;
    c2 =  x1 + x2^2 - x3^2 - 20 ;
    c  = [ c1 ; c2 ; ] ;
    
    ceq1 = -x1 - x2^2 + 2 ;
    ceq2 = x2 + 2*x3^2 - 3 ;
    ceq = [ ceq1 ; ceq2 ; ] ;    % ceq(x) = 0

end








