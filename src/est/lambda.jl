using Lasso

function σ̂_iter(y, X; γ=.05, c=1.1, ψ=.1, ϵ=.1, K=10, L=100, verbose=true)
  n, p = size(X)
  σ̂_k0 = Inf
  x̃ = map(j -> cor(y, X[:,j]), 1:p) |> A -> findmax(A)[2] |> j -> X[:,j]
  σ̂_k1 = ψ * std(y .- x̃*((x̃'*x̃)\(x̃'*y)), corrected=false)
  if verbose; println("σ̂ (k = ", 0,") = ", σ̂_k1); end
  ν = .2 * std(y)
  k = 1
  while (abs(σ̂_k1 - σ̂_k0) > ν) & (k <= K)
    σ̂_k0 = σ̂_k1
    λ = c * σ̂_k1 * sqrt(2*log(2*p)/n)
    fit = lasso(y, X, λ = λ, verbose = true)

    σ̂_k1 = mean((y .- fit.β0 - X * fit.β).^2)
    k = k + 1
    if verbose; println("σ̂ (k = ", k,") = ", σ̂_k1); end
  end
  return σ̂_k1
end

function _λ(y, X, c = 1.1, γ = .05)
  n, p = size(X)
  σ̂ = σ̂_iter(y, X; c=c, γ=γ)
  λ_thr = c * σ̂ * sqrt(2*log(2*p)/n)
  return λ_thr
end

# function Λ̂
