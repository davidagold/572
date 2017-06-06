module Lasso

using Distributions

export LassoFit, lasso

st(x::T, λ::T) where T = sign(x) * max(abs(x) - λ, zero(T))

function inner_prod(X::Matrix{T}, j::Int, Y::Vector{T}) where T<:Real
  res = zero(T)
  @inbounds @simd for i in 1:size(X, 1)
    res += X[i,j] * Y[i]
  end
  res
end
function inner_prod(X::AbstractSparseMatrix{T}, j::Int, Y::Vector{T}) where T<:Real
  res = zero(T)
  @inbounds for i in 1:size(X, 1)
    res += X[i,j] * Y[i]
  end
  res
end

struct LassoFit{T}
  λs::Vector{T}
  λ::T
  Υ::Vector{T}
  βs::Matrix{T}
  β::Vector{T}
  β0::T
  λid::Int
  standardized::Bool
end

function pre_lasso(y::AbstractVector{T}, X::AbstractArray{T}, standardize) where T<:Real
  n, p = size(X)
  ȳ = mean(y)
  y = y - ȳ
  μ_X = map(j -> mean(X[:,j]), 1:p)
  σ2_X = Vector{T}(p)
  for j in 1:p
    σ2_X[j] = mean(X[:,j].^2) - μ_X[j]^2
  end
  σ_X = sqrt.(σ2_X)
  if standardize
    X = X ./ σ_X'
    σ2_X = ones(p)
  end

  Xy = X' * y; XX = Matrix(X' * X)
  MM = Array{T}(p, p)
  μ_Xj_by_σ_Xj = μ_X ./ σ_X
  for j in 1:p
    t = μ_Xj_by_σ_Xj[j]
    for k in j:p
      MM[j,k] = n * t * μ_Xj_by_σ_Xj[k]
    end
  end
  MM = Symmetric(MM)
  # M = ones(T, n, p) .* (μ_X ./ σ_X)'
  XM = Array{T}(p, p)
  for j in 1:p
    Xjsum = sum(X[:,j])
    for k in 1:p
      XM[j,k] = Xjsum * μ_Xj_by_σ_Xj[k]
    end
  end
  My = Vector{T}(p)
  for j in 1:p
    My[j] = sum(y) * μ_Xj_by_σ_Xj[j]
  end
  # XM = X' * M
  # My = M' * y
  Ay = Xy - My
  AA = XX - XM - XM' + MM
  return ȳ, σ2_X, Ay, AA
end

# function pre_lasso(y::AbstractVector{T}, X::AbstractMatrix{T}, standardize) where T<:Real
#   n, p = size(X)
#   ȳ = mean(y)
#   y = y - ȳ
#   μ_X = map(j -> mean(X[:,j]), 1:p)
#   σ2_X = map(j -> std(X[:,j], mean = μ_X[j], corrected=false), 1:p)
#   σ_X = sqrt.(σ2_X)
#   if standardize
#     X = (X .- μ_X')./σ_X'
#     σ2_X = ones(T, p)
#   else
#     X = X .- μ_X'
#   end
#   Xy = X' * y; XX = X' * X
#   return ȳ, σ2_X, Xy, XX
# end

function lasso(y::AbstractVector{T}, X::AbstractMatrix{T};
               λ::T = zero(T), Υ::Vector{T} = ones(T, size(X,2)),
               ϵ::Float64 = -1., nλ::Int = 100, λs::Vector{T} = Vector{T}(),
               L::Int = 100, standardize::Bool = false,
               verbose::Bool = false) where T<:Real

  n, p = size(X)
  ȳ, σ2_X, Xy, XX = pre_lasso(y, X, standardize)
  ϵ = ifelse(ϵ < 0, ifelse(n > p, .0001, .01), ϵ)
  if isempty(λs)
    λ_max = maximum(abs, Xy)/n
    λ_min = ϵ * λ_max
    λs = exp.(range(log(λ_max), -(log(λ_max)-log(λ_min))/nλ, nλ))
  else
    nλ = length(λs)
  end
  if λ > 0
    λid = map(x -> abs(x-λ), λs) |> A -> findmin(A)[2]
    λs[λid] = λ
  else
    λ = maximum(λs)
    λid = 1
  end

  βs = zeros(T, p, nλ)
  for k in 2:nλ
    if verbose; println("k = ", k); end
    λ = λs[k]
    λΥ = λ .* Υ
    β = βs[:,k-1]
    @inbounds for l in 1:L
      for j in shuffle(1:p)
        xjr = Xy[j] - inner_prod(XX, j, β)
        Gj = xjr/n + σ2_X[j]*β[j]
        β[j] = st(Gj, λΥ[j])/σ2_X[j]
      end
    end
    βs[:,k] = β
  end
  return LassoFit(λs, λs[λid], Υ, βs, βs[:,λid], ȳ, λid, standardize)
end


function test(n=50, p=100, s=5; standardize=false, verbose = false)
  β = [ones(s); zeros(p-s)]
  X = rand(MvNormal(zeros(p), eye(p)), n)
  X = X'
  y = X * β + rand(Normal(0,.1), n)

  # lasso(y, sparse(X))
  lasso(y, X, verbose=verbose)
end

end # Lasso
