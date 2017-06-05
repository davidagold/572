module Fuller

using Distributions
export FullerFit, fuller, test

struct FullerFit{T<:Real}
  δ̂::Vector{T}
  σ̂::T
  Λ̂::Matrix{T}
end

# `Pr`oject
Pr(X, B, Xt = X') = X * ((Xt * X) \ (Xt * B))

function fuller(y::AbstractVector{T}, X::AbstractArray{T},
                Z::AbstractArray{T}, W::AbstractArray{T}, C = 1) where T<:Real
  y = y - Pr(W, y)
  X = X - Pr(W, X)
  Z = Z - Pr(W, Z)
  fuller(y, X, Z)
end

function fuller(y::AbstractVector{T}, X::AbstractArray{T},
                Z::AbstractArray{T}, C = 1) where T<:Real

  N, K = size(Z); G = size(X, 2)
  # if W != nothing # Partial out controls
  #   y = y - Pr(W, y)
  #   X = X - Pr(W, X)
  #   Z = Z - Pr(W, Z)
  # end
  Xt = X'; XX = Xt * X
  Zt = Z'; ZZ = Zt * Z

  Q = Matrix(ZZ) \ Matrix(Zt)
  Pdiag = zeros(T, N)
  for i in 1:N
    Zti = Zt[:,i]
    Qi = Q[:,i]
    for j in 1:K
      Pdiag[i] += Zti[j] * Qi[j]
    end
  end
  # P = Z * (Matrix(ZZ) \ Matrix(Zt)); Pdiag = diag(P)
  X̄ = [y X]; X̄t = X̄'
  X̄PX̄ = (X̄t * Z) * (Matrix(ZZ) \ (Matrix(Zt)*X̄))
  a = ((X̄t * X̄) \ X̄PX̄) |> x -> eigvals(x, scale=false) |> minimum
  â = (a - (1 - a) * C/(N - 1 - K)) / (1 - (1 - a) * C/(N - 1 - K))

  # Compute Fuller estimator
  PX = Pr(Matrix(Z), X); Py = Pr(Matrix(Z), y)
  XPX = Xt * PX
  # XZ = Xt * Z
  # QX = ZZ \ (Zt * X)
  # XPX = XZ * QX
  âXX = â * Xt*X
  # Qy = ZZ \ (Zt * y)
  XPy = Xt * Py
  # XPy = XZ * Qy
  âXy = â * Xt*y
  δ̂ = (XPX - âXX) \ (XPy - âXy)

  # StaNdard errors as given in [Hansen, Hausman, Newey 04; p4]
  # Note: Υ is an Upsilon; T is a type parameter, not the number of observations
  u = y - X * δ̂
  ut = u'; u2 = u.^2
  l2u = sum(u2)
  σ̂2 = l2u / (N - G)
  σ̂ = sqrt(σ̂2)

  # uZ = ut * Z
  # Qu = ZZ \ (Zt * u)
  # ã = uZ * Qu ./ l2u
  ã = ut * Pr(Matrix(Z), u) ./ l2u


  Υ = Pr(Matrix(Z), X);
  # Υ = Z * QX
  # Υ[i,j] = Zt[:,i] * Q[:,j]

  X̃ = X - u * (ut * X) ./ l2u
  X̃t = X̃'; PX̃ = Pr(Matrix(Z), X̃)
  V̂ = X̃ - PX̃
  κ = sum(Pdiag.^2)/N
  τ = K/N

  H = Xt * PX - ã .* XX
  Σ̂B = σ̂2*(1 - ã)^2 .* (X̃' * PX̃) + ã^2 .* X̃t * (X̃ - PX̃)
  # Compute Â (Note: Â0 is result of inner sum)
  Â0 = zeros(T, G)
  for j in 1:G
    V̂j = V̂[:,j]
    for i in 1:N
      Â0[j] += u2[i] * V̂j[i] / N
    end
  end
  Â = zeros(T, G, G)
  Pminusτ = Pdiag .- τ
  for j in 1:G
    Υj = Υ[:,j]
    for l in 1:G
      for i in 1:N
        Â[l,j] += Pminusτ[i] * Υj[i] * Â0[l]
        # Υi =
        # Â[l,j] += Pminusτ[i] * Υj[i] * Â0[l]
      end
    end
  end
  # Compute B̂
  B̂ = zeros(T, G, G)
  u2minusσ̂2 = u2 .- σ̂2
  ν = K * (κ-τ) / (N * (1 - 2*τ + κ*τ))
  for j in 1:G
    V̂j = V̂[:,j]
    for l in 1:G
      V̂l = V̂[:,l]
      for i in 1:N
        B̂[l,j] += u2minusσ̂2[i] * V̂j[i] * V̂l[i] * ν
      end
    end
  end
  # Compute Λ̂
  Ĥinv = inv(H)
  Σ̂ = Σ̂B + Â + Â' + B̂
  Λ̂ = Ĥinv * Σ̂ * Ĥinv
  return FullerFit([δ̂...], σ̂, Λ̂)
end

function test(N=100, G=10, K=50, s=5)
  _π = map(x -> .7^(x-1), 1:K)
  π = [ shuffle(_π) for j in 1:G ] |> x -> reduce(hcat,x)
  Z = rand(MvNormal(zeros(K), eye(K)), N)'
  V = rand(MvNormal(zeros(G), eye(G)), N)'
  X = Z * π .+ V
  δ = map(x -> .7^(x-1), 1:G)
  y = X * δ + rand(Normal(0,.1), N)
  W = ones(N)

  fuller(y, X, sparse(Z), W)
end

end # module Fuller
