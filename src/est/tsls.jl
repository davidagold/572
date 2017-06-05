module TSLS

struct TSLSFit{T<:Real}
  δ̂::Vector{T}
  σ̂::T
  Λ̂::Matrix{T}
end


function tsls(y, X, Z)
  n, K = size(Z)
  ZZ = Z' * Z
  Xt = X'
  num = Xt * Z * (ZZ \ Z'*y)
  denom = Xt * Z * (ZZ \ Z'*X)
  δ̂ = denom \ num

  σ̂2 = sum((y - X * δ̂).^2) / (n-1)
  σ̂ = sqrt(σ̂)
  Λ̂ = σ̂ * inv(denom)

  return TSLSFit(δ̂, σ̂, Λ̂)
end

end
