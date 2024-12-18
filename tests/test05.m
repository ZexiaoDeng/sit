

clc ;
clear ;
close all ;

Aineq = [] ;
bineq = [] ;
Aeq   = [] ;
beq   = [] ;
lb    = [ 0 ; 0 ; ] ;
ub    = [ 6 ; 6 ; ] ;

% options = optimoptions( 'fmincon'  , ...
%                         'Display'  ,'iter', ...
%                         'Algorithm','active-set', ...
%                         'PlotFcns' ,'optimplotfval' ) ;
options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','interior-point', ...
                        'PlotFcns' ,'optimplotfval' ) ;


x0 = [ 2 ; 4 ; ] ;
% x0 = [ 3.0 ; 2.0 ; ] ;
% x0 = [ 0 ; 0 ; ] ;

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

% ====================
% 目标函数
% ====================
function f = oracle( x )

    % 目标函数
    x1 = x(1) ;
    x2 = x(2) ;

    f = x2^2 + x1 ;
    
end

% ====================
% 非线性约束函数
% ====================
function [ c, ceq ] = nonlcon( x )
    % 非线性等式和不等式约束函数
    % c(x) <= 0
    % ceq( x ) = 0
    
    c = [  g_cst( x ) ; ...
          -h_cst( x ) ; ] ;
    ceq = [] ;    % ceq(x) = 0

end

function g = g_cst( x )

    x1 = x(1) ;
    x2 = x(2) ;
    g = max( [ ( x1 - 3 )^3 + x2 - 3 ; ...
                 x1^2 + x2^2 - 36    ; ] ) ;

end

function h = h_cst( x )

    x1 = x(1) ;
    x2 = x(2) ;
    
    h = 0.75*x1 + x2 - 5.5 ;
    
%     h = 0.75*x1^2 + x2^2 - 5.5 ;
    
end





