% ============================================
% Convex Analysis and Global Optimization, 2nd Edition, Hoang Tuy,
% Spring, 2016, P212
% =============================================
%   min        ( 3          + x1*x3      )*( x1*x2*x3*x4 + 2*x1*x3            + 2                            )^( 2/3 )
%
%   s.t.    -3*( 2*x1*x2    + 3*x1*x2*x4 )*( 2*x1*x3     + 4*x1*x4            - x2                           )
%           -  ( x1*x3      + 3*x1*x2*x4 )*( 4*x3*x4     + 4*x1*x3*x4         + x1*x3   - 4*x1*x2*x4         )^( 1/3 )
%           +3*( x4         + 3*x1*x3*x4 )*( 3*x1*x2*x3  + 3*x1*x4            + 2*x3*x4 - 3*x1*x2*x4         )^( 1/4 ) <= -309.219315
%
%           -2*( 3*x3       + 3*x1*x2*x3 )*( x1*x2*x3    + 4*x2*x4            - x3*x4                        )^2
%           +  ( 3*x1*x2*x3              )*( 3*x3        + 2*x1*x2*x3         + 3*x4                         )^4 - ( x2*x3*x4 + x1*x3*x4 )*( 4*x1 - 1 )^( 3/4 )
%           -3*( 3*x3*x4    + 2*x1*x3*x4 )*( x1*x2*x3*x4 + x3*x4 - 4*x1*x2*x3 - 2*x1                         )^4       <= -78243.910551
%
%           -3*( 4*x1*x3*x4              )*( 2*x4        + 2*x1*x2            - x2      - x3                 )^2
%           +2*( x1*x2*x4   + 3*x1*x3*x4 )*( x1*x2       + 2*x2*x3            + 4*x2    - x2*x3*x4   - x1*x3 )^4 <= 9618
%
%           0 <= xi <= 5, i = 1, ..., 4
%
clc ;
clear ;
close all ;

Aineq = [] ;
bineq = [] ;
Aeq   = [] ;
beq   = [] ;
lb    =  zeros( 4, 1 ) ;
ub    = 5*ones( 4, 1 ) ;

options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','sqp-legacy', ...
                        'PlotFcns' ,'optimplotfval' ) ;
x0 = 5*ones( 4, 1 ) ;

% options = optimoptions( 'fmincon'  , ...
%                         'Display'  ,'iter', ...
%                         'Algorithm','interior-point', ...
%                         'MaxFunctionEvaluations', 1e4, ...
%                         'MaxIterations', 1e4, ...
%                         'OptimalityTolerance', 1e-8 , ...
%                         'ConstraintTolerance', 1e-15, ...
%                         'StepTolerance', 1e-10, ...
%                         'PlotFcns' ,'optimplotfval' ) ;
% x0    = zeros( 4, 1 ) ;          % 初值


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
% 目标函数
function f = oracle( x )
    % 目标函数
    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    x4 = x(4) ;
    f = ( 3 + x1*x3 )*( x1*x2*x3*x4 + 2*x1*x3 + 2)^( 2/3 ) ;
    
end

% 非线性约束函数
function [ c, ceq ] = nonlcon( x )
    % 非线性等式和不等式约束函数
    % c(x) <= 0
    % ceq( x ) = 0
    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    x4 = x(4) ;
    
    c1 =  -3*( 2*x1*x2    + 3*x1*x2*x4 )*( 2*x1*x3     + 4*x1*x4            - x2                           ) ...
          -  ( x1*x3      + 3*x1*x2*x4 )*( 4*x3*x4     + 4*x1*x3*x4         + x1*x3   - 4*x1*x2*x4         )^( 1/3 ) ...
          +3*( x4         + 3*x1*x3*x4 )*( 3*x1*x2*x3  + 3*x1*x4            + 2*x3*x4 - 3*x1*x2*x4         )^( 1/4 )  + 309.219315 ;

    c2 =  -2*( 3*x3       + 3*x1*x2*x3 )*( x1*x2*x3    + 4*x2*x4            - x3*x4                        )^2 ...
          +  ( 3*x1*x2*x3              )*( 3*x3        + 2*x1*x2*x3         + 3*x4                         )^4 - ( x2*x3*x4 + x1*x3*x4 )*( 4*x1 - 1 )^( 3/4 ) ...
          -3*( 3*x3*x4    + 2*x1*x3*x4 )*( x1*x2*x3*x4 + x3*x4 - 4*x1*x2*x3 - 2*x1                         )^4        + 78243.910551 ;

    c3 =  -3*( 4*x1*x3*x4              )*( 2*x4        + 2*x1*x2            - x2      - x3                 )^2 ...
          +2*( x1*x2*x4   + 3*x1*x3*x4 )*( x1*x2       + 2*x2*x3            + 4*x2    - x2*x3*x4   - x1*x3 )^4 - 9618 ;
      
    c = [ c1 ; c2 ; c3 ; ] ;
    
    ceq = [] ;    % ceq(x) = 0

end