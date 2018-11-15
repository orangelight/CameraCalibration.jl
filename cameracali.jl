"""
Based on:
Wilhelm Burger: Zhang’s Camera Calibration Algorithm: In-Depth Tutorial and Implementation,
Technical Report HGB16-05, University of Applied Sciences Upper Austria, School of
Informatics, Communications and Media, Dept. of Digital Media, Hagenberg, Austria, May
2016. DOI: 10.13140/RG.2.1.1166.1688, http://staff.fh-hagenberg.at/burger/
"""

using Statistics, LinearAlgebra
"""
estimateHomography - returns the estimated homography matrix H, such that b_j = H*a_j.
'''
A and B are in the form 3xN (each row is a point)
'''
"""
function estimateHomography(a::Array{<:Real,2}, b::Array{<:Real,2})
    if size(a) != size(b)
        error("src and dst must have same dimensions")
    end
    n = size(a)[2] #number of points
    #normalisation matrices
    nA = getNormalisationMatrix(a)
    nB = getNormalisationMatrix(b)
    m = zeros((2*n, 9))

    #construct matrix for solving
    for j  = 1:n
        k = 2*j
        ap = nA*a[:,j]
        bp = nB*b[:,j]
        # M_{k,*} = (x′_a, y′_a, 1, 0, 0, 0, −x′_ax′_b , −y′_ax′_b, −x′_b )
        m[k-1,:] = [ap[1] ap[2] 1 0 0 0 -ap[1]*bp[1] -ap[2]*bp[1] -bp[1]]
        # M_{k+1,*} = (0, 0, 0, x′_a , y′_a , 1, −x′_ay′_b, −y′_ay′_b , −y′_b )
        m[k,:] = [0 0 0 ap[1] ap[2] 1 -ap[1]*bp[2] -ap[2]*bp[2] -bp[2]]
    end

    #solve
    u, d, v = svd(m, full = true)
    h = v[:,9] #last column of v
    hp = reshape(h, 3, 3)'
    return inv(nB)*hp*nA
end
"""
getNormalisationMatrix - returns normalisation matrix of points x
'''

'''
"""
function getNormalisationMatrix(x::Array{<:Real, 2})
    m = mean(x[1:2, :], dims = 2)
    v = var(x[1:2, :], dims = 2, corrected = false, mean = m)
    #scale values (182)
    sx = sqrt(2.0/v[1])
    sy = sqrt(2.0/v[2])
    #construct normalisation matrix (181)
    n = [sx 0 -sx*m[1]; 0 sy -sy*m[2]; 0 0 1]
    n
end

function refineHomography(h::Array{<:Real, 2}, a::Array{<:Real, 2}, b::Array{<:Real, 2})

end


function calibrate(x::Array{<:Real,2}, u::Array{Array{<:Real,2},1})
    hListInit = getHomographies(x,u)
    aInit = getCameraIntrinsics(hListInit)
    wInit = getExtrinsics(aInit,hListInit)
end

function getHomographies(x::Array{<:Real,2}, u::Array{Array{<:Real,2},1})
    m = size(u)[2]
    hList = Array{Float64,2}[]
    for i = 1:m
        hInit = estimateHomography(x, u[i])
        # TODO refineHomography
        push!(hList, hInit)
    return hList
end

function getCameraIntrinsics(hs::Array{Array{<:Real,2},1})
    m = size(hs)[1]
    v = zeros((2*m, 6))

    for i = 1:m
        v[2*i-1,:] = getIntrinsicRowVector(1,2,hs[i])
        v[2*i,:] = getIntrinsicRowVector(1,1,hs[i]) - getIntrinsicRowVector(2,2,hs[i])
    end

    #solve
    u, s, v = svd(m, full = true)
    b = reshape(v[:,9], 3, 3)'

    w = b[1]*b[3]*b[6] - (b[2]^2)*b[6] - b[1]*(b[5]^2) + 2*b[2]*b[4]*b[5] - b[3]*(b[4]^2) # (104)
    d = b[1]*b[3] - b[2]^2 # (105)
    α = sqrt(w/(d*b[1])) # (99)
    β = sqrt(w/d^2 *b[1]) # (100)
    γ = sqrt(w/(d^2 * b[1])) * b[2] # (101)
    uc = (b[2]*b[5] - b[3]*b[4])/d
    vc = (b[2]*b[4] - b[1]*b[5])/d
    return [α γ uc
            0 β vc
            0 0 1]
end

function getIntrinsicRowVector(p::Int64, q::Int64,  h::Array{<:Real,2})
    return [h[1,p]*h[1,q]
            h[1,p]*h[2,q]+h[2,p]*h[1,q]
            h[2,p]*h[2,q]
            h[3,p]*h[1,q]+h[1,p]*h[3,q]
            h[3,p]*h[2,q]+h[2,p]*h[3,q]
            h[3,p]*h[3,q]]'
end

function getExtrinsics(a::Array{<:Real,2},hs::Array{Array{<:Real,2},1})

end
"""
tests
asample = [182 535 171 537; 350 358 553 563; 1 1 1 1]
bsample = [0 888 0 888; 0 0 500 500; 1 1 1 1]
"""
