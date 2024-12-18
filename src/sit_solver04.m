function [ x, fval, output ] = sit_solver04( GDC )


% =====================
% 求解器参数设置
% =====================
epsilon = 1e-3 ;
eta     = 1e-3 ; 
options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'none', ...
                        'Algorithm','interior-point' ) ;
path = './bt-1.3' ;
addpath( path ) ;
% 输出参数      
M1 = GDC.M ;
a  = M1.p ;
b  = M1.q ;

n = length( a ) ;           % 决策变量 x 的维数

[ w, gamma0, output ] = set_gamma0( M1 ) ;

if output.exitflag == 1
    xBar   = w                   ;
    gamma0 = f_obj( xBar ) - eta ;
end

% ======================
% Step 0
% ======================
gamma = gamma0 ;

node.M      = M1   ;
node.beta_M = -inf ;
node.xk     = a    ;
node.yk     = b    ; 
                    
Pk_Stack(1) = { node } ;       % Mk 剖分后得到的集合 collection
Rk_Stack    = {}       ;       % 剩余集合 collection

k = 0 ;
branch = [] ;


while ( true && k <= 1e4 )
    
    k = k + 1
    
    % =====================================
    % Step 1: box reduction & bounding
    % =====================================
    % (1) 盒子收缩( box reduction )
    while ( ~isempty( Pk_Stack ) )
        
        len_P             = length( Pk_Stack ) ;
        node              = Pk_Stack{ len_P }  ;
        Pk_Stack( len_P ) = []                 ;
        
        M    = node.M ;
        p    = M.p    ;
        q    = M.q    ;
        
        % ====================
        % 基本凹规划( BCP )
        % ====================
        rep.l = p ;
        rep.u = q ;

        P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep   
        CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
        
        yopt = zeros( size( CH.V, 2 ), 1 ) ;
        for idx_k = 1: size( CH.V, 2 )
            yopt( idx_k, 1 ) = g2_cst( CH.V( :, idx_k ) ) ;
        end
        [ ~, idx_yopt ] = max( yopt ) ;
        y = CH.V( :, idx_yopt ) ;
        
        gBarx   = @(x) g1_cst(x) - g2_cst(y) ;
        
        pp = zeros( n, 1 ) ;
        qq = zeros( n, 1 ) ;
        
        for idx_j = 1: n
            % ====================
            % 凸规划( CP )
            % ====================
            [ ~       , ...
              ppi     , ...
              exitflag1 ] = fmincon( @(x)  x( idx_j ),     ...
                                    3*ones( n, 1 )  ,     ...
                                    []              , [], ...
                                    []              , [], ...
                                    p               , q , ...
                                    @piqi_cst       ,     ...
                                    options             ) ;
            if exitflag1 <= 0
                % 剪枝
                break ;
            end
            % ====================
            % 凸规划( CP )
            % ====================
            [ ~       , ...
              qqi     , ...
              exitflag2 ] = fmincon( @(x) -x( idx_j ),     ...
                                   3*ones( n, 1 )  ,     ...
                                   []              , [], ...
                                   []              , [], ...
                                   p               , q , ...
                                   @piqi_cst       ,     ...
                                   options               ) ;
           if exitflag2 <= 0
                % 剪枝
                break ;
           end
           
           pp( idx_j ) = min( ppi, -qqi ) ;     % 避免出现数值问题
           qq( idx_j ) = max( ppi, -qqi ) ;

        end
        
        if exitflag1<=0 || exitflag2 < 0
            % 剪枝
            continue ;
        end
        
        redM.p = pp ;
        redM.q = qq ;
        
        % ===============================================
        % (2) 定界( bounding )
        %     交替方法计算下边界 beta_M ( lower_bound )
        % ===============================================
        % 一般凹规划( GCP )
        % ====================
        rep.l = redM.p ;
        rep.u = redM.q ;

        P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep   
        CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
        
        opt.color = [ 0.5, 1, 0.05*k ] ;
        plot( P, opt ) ;
        drawnow
        
        yMopt = zeros( size( CH.V, 2 ), 1 ) ;
        for idx_k = 1: size( CH.V, 2 )
            yMopt( idx_k, 1 ) = g2_cst( CH.V( :, idx_k ) ) ;
        end
        [ ~, idx_yMopt ] = max( yMopt ) ;
        yM = CH.V( :, idx_yMopt ) ;
        
%         [ yM, ...
%           ~ , ...
%           exitflag ] = fmincon( @(x) -g2_cst(x),           ...
%                              3*ones( n, 1 ) ,           ...
%                              []             , []      , ...
%                              []             , []      , ...
%                              redM.p         , redM.q  , ...
%                              @f_cst         ,           ...
%                              options                    ) ;
%         if exitflag <= 0
%             % 剪枝情况
%             continue ;
%         end

        % ====================
        % 凸规划( CP )
        % ====================
        [ xM, ...
          ~ , ...
          exitflag ] = fmincon( @(x)  g1_cst(x),           ...
                             3*ones( n, 1 ) ,           ...
                             []             , []      , ...
                             []             , []      , ...
                             redM.p         , redM.q  , ...
                             @f_cst         ,           ...
                             options                    ) ;
        if exitflag <= 0
            % 剪枝情况
            continue ;
        end
        
        beta_M = g1_cst( xM ) - g2_cst( yM ) ;  % 上边界 beta_M
        
        % ========================
        % Step 2: 剪枝和收集
        % ========================
        if beta_M > -epsilon
            % 剪枝操作
            continue ;
        else
            % 注意更新
            node.M            = redM  ;
            node.xk           = xM ;
            node.yk           = yM ;
            node.beta_M       = beta_M ;
            
            len_R = length( Rk_Stack ) ;
            Rk_Stack( len_R + 1 ) = { node } ;
        end % if
    end % while
    
    % ==================================================
    % Step 3: 停止准则( terminate )
    % (1) 问题(P)是 epsilon-实质不可行情况
    % (2) xBar 为问题(P)的 ( epsilon, eta )-实质最优解
    % ==================================================
    if ( isempty( Rk_Stack ) )
        if ( gamma == gamma0 )
            fprintf( 'problem (P) is epsilon-infeasible!\n' ) ;
            x               = inf ;
            fval            = nan ;
            output.exitflag = -1  ;
            output.message  = 'problem (P) is epsilon-infeasible!' ;
            break ;
        else
            x               = xBar          ;
            fval            = f_obj( xBar ) ;
            output.exitflag = 1             ;
            output.message  = '( epsilon, eta )-optimal solution is found!' ;
            break ;
        end
    end
    
    % ==================================================
    % Step 4: 
    % 判断 xk 是否满足 f( xk ) >= gamma
    % 最优界优先搜索
    % ==================================================
    len_R     = length( Rk_Stack )
    
    branch = [ branch ; len_R ] ;
    
    beta_list = zeros( len_R, 1 ) ;
    for idx = 1: len_R
        beta_list = Rk_Stack{ idx }.beta_M ;
    end
    [ betak, idx_opt ] = min( beta_list ) ;
    nodek = Rk_Stack{ idx_opt } ;       % 提取最优节点
    
    Rk_Stack( idx_opt ) = [] ;         % 删除已经分析节点
    
    xk = nodek.xk
    yk = nodek.yk
    Mk = nodek.M ;
    
    % 停止准则
    if ( g1_cst( xk ) - g2_cst( xk ) ) < -epsilon
        % ===============
        % Step 5
        % ===============
        if ( gamma == gamma0 ) || ( f_obj( xk ) < f_obj( xBar ) )
            xBar  = xk ;
            gamma = f_obj( xBar ) - eta ;     % 更新 gamma
        end % if
    end % if
    
    % =========================
    % Step 6
    % 分支操作( bounding )
    % 自适应对分矩形算法
    % ==========================
    [ ~, jk ] = max( abs( yk - xk ) ) ;
    vk = ( xk + yk )/2 ;
    
    M_minus.p = Mk.p ;
    M_minus.q = Mk.q ;
    M_minus.q( jk ) = vk( jk ) ;
    
    M_plus.p = Mk.p ;
    M_plus.q = Mk.q ;
    M_plus.p( jk ) = vk( jk ) ;
    
    left_child_node  = nodek ;
    right_child_node = nodek ;
    left_child_node.M  = M_minus ;
    right_child_node.M = M_plus  ;
    
    Pk_Stack( 1 ) = { left_child_node  } ;
    Pk_Stack( 2 ) = { right_child_node } ;
    
    rep.l = M_plus.p ;
    rep.u = M_plus.q ;

    P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep   
    CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep

    opt.color = [ 0.05*k, 0.5, 0.2 ] ;
    plot( P, opt ) ;
    drawnow
    
end % while

plot( branch, '-rs' ) ;
drawnow

x               = xBar          ;
fval            = f_obj( xBar ) ;
output.k        = k             ;
output.exitflag = 1             ;
output.message  = '( epsilon, eta )-optimal solution is found!' ;

% ====================
% 非线性约束函数
% ====================
function [ c, ceq ] = piqi_cst( x )
    % 非线性等式和不等式约束函数
    c   = [ f_obj(x) - gamma   ; ...
            gBarx(x) - epsilon ; ] ;
    ceq = [] ;    % ceq(x) = 0

end

function [ c, ceq ] = f_cst( x )
    % 非线性等式和不等式约束函数
    c   = f_obj(x) - gamma ;
    ceq = [] ;

end



end


