module CVLasso

using Lasso
using Distributions
export cvlasso

function cvlasso(y::AbstractVector{T}, X::AbstractMatrix{T};
                  nfolds::Int = 10, λs::Vector{T} = Vector{T}(),
                  Υ::Vector{T} = ones(T, size(X,2)),
                  ϵ::Float64 = .1, nλ::Int = 20, L::Int = 100,
                  standardize = false, verbose::Bool = true) where T <: Real
  #something
  n, p = size(X)
  ȳ, σ2_X, Xy, XX = Lasso.pre_lasso(y, X, standardize)

  if isempty(λs)
    λ_max = λ_max = maximum(abs, Xy)/n
    λ_min = ϵ * λ_max
    λs = exp.(range(log(λ_max), -(log(λ_max)-log(λ_min))/nλ, nλ))
  end

  perfold::Int = floor(n / nfolds)
  ids = shuffle(1:n)
  folds = [ ids[(a * perfold + 1):min((a + 1) * perfold, n)] for a in 0:(nfolds-1) ]

  fits = Vector{LassoFit}(nfolds)
  for r in 1:nfolds
    train = setdiff(1:n, folds[r])
    fits[r] = lasso(y[train], X[train,:]; λs = λs, Υ = Υ, standardize=standardize)
    if verbose; info(@sprintf("%2.f%% done fitting", 100*r/nfolds)); end
  end
  mses = [ mean((y[folds[r]] .- fits[r].β0 .- X[folds[r],:] * fits[r].βs[:,k]).^2)
           for r in 1:nfolds, k in 1:nλ ]
  mse_min, λcv_id = map(k->mean(mses[:,k]), 1:nλ) |> findmin

  βs = zeros(T, p, nλ)
  for k in 2:nλ
    if verbose; println("k = ", k); end
    λ = λs[k]
    λΥ = λ .* Υ
    β = βs[:,k-1]
    @inbounds for l in 1:L
      for j in shuffle(1:p)
        xjr = Xy[j] - Lasso.inner_prod(XX, j, β)
        Gj = xjr/n + σ2_X[j]*β[j]
        β[j] = Lasso.st(Gj, λΥ[j])/σ2_X[j]
      end
    end
    βs[:,k] = β
  end
  return LassoFit(λs, λs[λcv_id], Υ, βs, βs[:,λcv_id], ȳ, λcv_id, standardize)
end

function test(n=50, p=100, s=5; standardize=false, ϵ=.1)
  β = [ones(s); zeros(p-s)]
  X = rand(MvNormal(zeros(p), eye(p)), n)
  X = X'
  y = X * β + rand(Normal(0,.1), n)

  cvlasso(y, X, Υ = rand(p), ϵ=ϵ)
end

end
