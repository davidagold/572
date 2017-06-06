using Lasso

function σ̂_iter(y::AbstractVector{T}, X::AbstractArray{T};
                γ=.05, c=1.1, ψ=.1, ϵ=.1, K=10, L=100, verbose=true) where T<:Real
  n, p = size(X)
  σ̂_k0 = Inf
  x̃ = map(j -> cor(y, X[:,j]), 1:p) |> A -> findmax(A)[2] |> j -> X[:,j]
  σ̂_k1 = ψ * std(y .- x̃*((x̃'*x̃)\(x̃'*y)), corrected=false)
  if verbose; println("σ̂ (k = ", 0,") = ", σ̂_k1); end
  ν = .2 * std(y)
  k = 1

  nλ = 50
  ȳ, σ2_X, Xy, XX = Lasso.pre_lasso(y, X, false)
  λ_max = maximum(abs, Xy)/n
  λ_min = ϵ * λ_max
  λs = exp.(range(log(λ_max), -(log(λ_max)-log(λ_min))/nλ, nλ))


  while (abs(σ̂_k1 - σ̂_k0) > ν) & (k <= K)
    σ̂_k0 = σ̂_k1
    λ = c * σ̂_k1 * sqrt(2*log(2*p)/n)
    λid = map(x -> abs(x-λ), λs) |> A -> findmin(A)[2]
    # λs[λid] = λ
    splice!(λs, λid+1:λid, λ)

    βs = zeros(T, p, nλ)
    for k in 2:λid
      if verbose; println("k = ", k); end
      λ = λs[k]
      β = βs[:,k-1]
      @inbounds for l in 1:L
        for j in shuffle(1:p)
          xjr = Xy[j] - Lasso.inner_prod(XX, j, β)
          Gj = xjr/n + σ2_X[j]*β[j]
          β[j] = Lasso.st(Gj, λ)/σ2_X[j]
        end
      end
      βs[:,k] = β
    end
    β = βs[λid]
    σ̂_k1 = mean((y .- ȳ .- X * β).^2)
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
