% ====================
% g1 约束函数
% ====================
function g1 = g1_cst( x )
    
    x1 = x(1) ;
    x2 = x(2) ;
    x3 = x(3) ;
    
    g1 = ( x1 - 3 )^2 + ( x2 - 2 )^2 ;
    
end