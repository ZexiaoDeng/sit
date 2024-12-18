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
lb    = [ 0 ; 0 ; 0 ; -inf ; ] ;
ub    = [ 1 ; 5 ; 2 ; 4 ; ] ;

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
                        'PlotFcns' ,'optimplotfval' ) ;
                    
x0 = ones( 4, 1 ) ;

% x0 = [ 0.898026312090722 ;
%    4.806451257207263 ;
%    2.000000000000000 ; ...
%    -4 ] ;

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
function f1 = oracle( x )
    % 目标函数
    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    t  = x(4) ;
    
    f1 = x1^2 + x2^2 - t + 8 ;
    
    return ;
    
end

% ====================
% 非线性约束函数
% ====================
function g = g_cst( x )

    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    t  = x(4) ;

    g1 = x2 ;               h1 = x1^2 + x3^2 ;
    g2 = x3^2 - x1 + 20 ;   h2 = x2^2 ;
    g3 = -x1 + 2 ;          h3 = x2^2 ;
    g4 = -x2 + 3 ;          h4 = 2*x3^2 ;
    g5 = t ;                h5 = x3^2 ;

%     h  = h1 + h2 + h3 + h4 + h5 ;
%     [ g, idx ] = max( [ h1 + h2 + h3 + h4 + h5 ; ...
%                         g1 + h2 + h3 + h4 + h5 ; ...
%                         h1 + g2 + h3 + h4 + h5 ; ...
%                         h1 + h2 + g3 + h4 + h5 ; ...
%                         h1 + h2 + h3 + g4 + h5 ; ...
%                         h1 + h2 + h3 + h4 + g5 ; ] ) ;

    [ g, idx ] = max( [ g1 + h2 + h3 + h4 + h5 ; ...
                        h1 + g2 + h3 + h4 + h5 ; ...
                        h1 + h2 + g3 + h4 + h5 ; ...
                        h1 + h2 + h3 + g4 + h5 ; ...
                        h1 + h2 + h3 + h4 + g5 ; ] ) ;
    return ;
    
end

function h = h_cst( x )

    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    t  = x(4) ;

    g1 = x2 ;               h1 = x1^2 + x3^2 ;
    g2 = x3^2 - x1 + 20 ;   h2 = x2^2 ;
    g3 = -x1 + 2 ;          h3 = x2^2 ;
    g4 = -x2 + 3 ;          h4 = 2*x3^2 ;
    g5 = t ;                h5 = x3^2 ;

    h  = h1 + h2 + h3 + h4 + h5 ;
%     [ g, idx ] = max( [ h1 + h2 + h3 + h4 + h5 ; ...
%                         g1 + h2 + h3 + h4 + h5 ; ...
%                         h1 + g2 + h3 + h4 + h5 ; ...
%                         h1 + h2 + g3 + h4 + h5 ; ...
%                         h1 + h2 + h3 + g4 + h5 ; ...
%                         h1 + h2 + h3 + h4 + g5 ; ] ) ;
%     [ g, idx ] = max( [ g1 + h2 + h3 + h4 + h5 ; ...
%                         h1 + g2 + h3 + h4 + h5 ; ...
%                         h1 + h2 + g3 + h4 + h5 ; ...
%                         h1 + h2 + h3 + g4 + h5 ; ...
%                         h1 + h2 + h3 + h4 + g5 ; ] ) ;
    return ;
    
end

% ====================
% 非线性约束函数
% ====================
function [ c, ceq ] = nonlcon( x )
    
    c  = g_cst(x) - h_cst(x) ;
    
    ceq = [] ;    % ceq(x) = 0

end











