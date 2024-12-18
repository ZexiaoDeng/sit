% ==========================================
%   min     f(x) = x1^2 + x2^2 - x3^2 + 8
%   s.t.    -( x1^2 - x2 + x3^2 ) <= 0
%           -( x1 + x2^2 - x3^2 - 20 ) <= 0
%           -x1 - x2^2 + 2 <= 0
%           -( x2 + 2*x3^2 - 3 ) <= 0
%           x1, x2, x3 >= 0
%
% x =
%    0.898026312177311
%    4.806451257198426
%    2.000000000000000
% 
% fval =
%   27.908424945187111
% ===========================================

clc ;
clear ;
close all ;

Aineq = [] ;
bineq = [] ;
Aeq   = [] ;
beq   = [] ;
lb    = zeros( 3, 1 ) ;
ub    = [ 1 ; 5 ; 2 ; ] ;

options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','sqp-legacy', ...
                        'PlotFcns' ,'optimplotfval' ) ;
% options = optimoptions( 'fmincon'  , ...
%                         'Display'  ,'iter', ...
%                         'Algorithm','active-set', ...
%                         'PlotFcns' ,'optimplotfval' ) ;
options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','interior-point', ...
                        'MaxFunctionEvaluations', 1e4, ...
                        'MaxIterations', 1e4, ...
                        'OptimalityTolerance', 1e-6 , ...
                        'ConstraintTolerance', 1e-15, ...
                        'StepTolerance', 1e-10, ...
                        'PlotFcns' ,'optimplotfval' ) ;
                    
x0 = ones( 3, 1 ) ;
% x0 = [ 0.5 ; 4.5 ; 2 ; ] ;

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

[ c, ceq ] = nonlcon( x ) ;
c

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
    c2 =  -( x1 + x2^2 - x3^2 - 20 ) ;
    c3 = -x1 - x2^2 + 2 ;
    c4 = -( x2 + 2*x3^2 - 3 ) ;
    
    c  = [ c1 ; c2 ; c3 ; c4 ; ] ;
    
    ceq = [] ;    % ceq(x) = 0

end











