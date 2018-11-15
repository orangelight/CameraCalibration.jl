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
	println(h)
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

"""
tests
asample = [182 535 171 537; 350 358 553 563; 1 1 1 1]
bsample = [0 888 0 888; 0 0 500 500; 1 1 1 1]
"""
