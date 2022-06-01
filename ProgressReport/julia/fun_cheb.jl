# CHEB  compute D = differentiation matrix, x = Chebyshev grid

function fun_cheb(N)
  if N==0
    D=0
    x=1
    return
  end
    x = cos.(pi*(0:N)/N)'
    c = [2; ones(N-1,1); 2].*(-1).^[0:N;]
    X = repeat(x',1, N+1)
    dX = X - X'
    D  = (c*(1 ./ c)')./(dX+(I(N+1)))  # off-diagonal entries
    D  = D - Diagonal(vec(sum(D',dims=1)))                 # diagonal entries
  return D,x
end
