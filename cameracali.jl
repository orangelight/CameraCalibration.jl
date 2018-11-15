using Statistics, LinearAlgebra
"""
estimateHomography - returns the estimated homography matrix H, such that b_j = H*a_j.
'''
A and B are in the form 3xN (each row is a point)
'''
"""
function estimateHomography(a::Array{Real,2}, b::Array{Real,2})
    if size(a) != size(b)
        error("src and dst must have same dimensions")
    end
    n = size(a)[2] #number of points
    #normalisation matrices
    nA = getNormalisationMatrix(a)
    nB = getNormalisationMatrix(a)
    m = zeros((2*n, 9))

    for j  = 1:(n-1)
        k = 2*j
        ap = nA*a[:,n]
        ab = nB*b[:,n]
        m[:,k-1] = [ap[1] ap[2] 1 0 0 0 -ap[1]*bp[1] -ap[2]*bp[1] -bp[1]]
        m[:,k] = [0 0 0 ap[1] ap[2] 1 -ap[1]*bp[2] -ap[2]*bp[2] -bp[2]]
    end

    #solve
    u, d, vt = svd(m)
end

function getNormalisationMatrix(x::Array{Real, 2})
    m = mean(x[1:2, :], dims = 2)
    v = var(x[1:2, :], dims = 2, corrected = false, mean = m)
    #cx = zeros(size(x))
    #cx[1,:] = x[1,:] .- m[1]
    #cx[2,:] = x[2,:] .- m[2]
    #cx[3,:] .= 1
    sx = sqrt(2.0/v[1])
    sy = sqrt(2.0/v[2])
    #avgdist = mean(sqrt.(cx[1,:].^2 .+ cx[2,:].^2))
    #s = sqrt(2.0)/avgdist
    #normalization matrix (181)
    #N = [s 0 -s*m[1]; 0 s -s*m[2]; 0 0 1]
    n = [sx 0 -sx*m[1]; 0 sy -sy*m[2]; 0 0 1]
    n
end
