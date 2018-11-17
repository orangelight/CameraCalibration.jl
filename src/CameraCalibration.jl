"""
Based on:
Wilhelm Burger: Zhang’s Camera Calibration Algorithm: In-Depth Tutorial and Implementation,
Technical Report HGB16-05, University of Applied Sciences Upper Austria, School of
Informatics, Communications and Media, Dept. of Digital Media, Hagenberg, Austria, May
2016. DOI: 10.13140/RG.2.1.1166.1688, http://staff.fh-hagenberg.at/burger/
"""
module CameraCalibration
using Statistics, LinearAlgebra
export estimateHomography, calibrate, getHomographies, getCameraIntrinsics, getExtrinsics, getCameraIntrinsicsB
"""
estimateHomography - returns the estimated homography matrix H, such that b_j = H*a_j.
'''
A and B are in the form 3xN (each column is a point)
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
    H = [h[1] h[2] h[3];
         h[4] h[5] h[6]
         h[7] h[8] h[9]]

    H = inv(nB)*H*nA
    return H./H[3,3]
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
    return [sx 0 -sx*m[1]; 0 sy -sy*m[2]; 0 0 1]
end

function refineHomography(h::Array{<:Real, 2}, a::Array{<:Real, 2}, b::Array{<:Real, 2})

end


function calibrate(x::Array{<:Real,2}, u::Array{Array{T1,2},1}) where T1 <: Real
    hListInit = getHomographies(x,u)
    aInit = getCameraIntrinsics(hListInit)
    #wInit = getExtrinsics(aInit,hListInit)
end

function getHomographies(x::Array{<:Real,2}, u::Array{Array{T1,2},1}) where T1 <: Real
    m = size(u)[1]
    hList = Array{Float64,2}[]
    for i = 1:m
        hInit = estimateHomography(x, u[i])
        # TODO refineHomography
        push!(hList, hInit)
    end
    return hList
end

function getCameraIntrinsics(hs::Array{Array{T1,2},1}) where T1 <: Real
    m = size(hs)[1]
    V = zeros((2*m, 6))

    for i = 1:m
        V[2*i-1,:] = getIntrinsicRowVector(1,2,hs[i])
        V[2*i,:] = getIntrinsicRowVector(1,1,hs[i]) - getIntrinsicRowVector(2,2,hs[i])
    end

    #solve
    u, s, v = svd(V, full = true)
    print(size(v))
    b = v[:,6]

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

function getCameraIntrinsicsB(hs::Array{Array{T1,2},1}) where T1 <: Real
    m = size(hs)[1]
    V = zeros((2*m, 6))

    for i = 1:m
        V[2*i-1,:] = getIntrinsicRowVector(1,2,hs[i])
        V[2*i,:] = getIntrinsicRowVector(1,1,hs[i]) .- getIntrinsicRowVector(2,2,hs[i])
    end

    #solve
    u, s, v = svd(V, full = true)
    b = v[:,6]
    if b[1] < 0 || b[3] < 0 || b[6] < 0
        b = -b
    end
    bInit = [b[1] b[2] b[4]
             b[2] b[3] b[5]
             b[4] b[5] b[6]]
    C = cholesky(bInit)
    return inv(C.L)' * C.L[3,3]
end

function getIntrinsicRowVector(p::Int64, q::Int64,  h::Array{<:Real,2})
    return [h[1,p]*h[1,q]
            h[1,p]*h[2,q]+h[2,p]*h[1,q]
            h[2,p]*h[2,q]
            h[3,p]*h[1,q]+h[1,p]*h[3,q]
            h[3,p]*h[2,q]+h[2,p]*h[3,q]
            h[3,p]*h[3,q]]'
end


"""
a is intrinsic matrix hs is list of homographies

"""
function getExtrinsics(a::Array{<:Real,2},hs::Array{Array{T1,2},1}) where T1 <: Real
    wList = Array{Float64,2}[]
    m = size(hs)[1]

    for i = 1:m
        push!(wList, estimateViewTransorm(a, hs[i]))
    end
    return wList
end

function estimateViewTransorm(a::Array{<:Real,2}, h::Array{<:Real,2})
    κ = 1/norm(inv(a)*h[:,1])
    r0 = κ*inv(a)*h[:,1]
    r1 = κ*inv(a)*h[:,2]
    r2= cross(r0,r1)
    R = zeros((3, 3))
    R[:,1] = r0
    R[:,2] = r1
    R[:,3] = r2
    # make true rotation matrix
    u, s, v = svd(R, full = true)
    R = v*u'
    t = κ*inv(a)*h[:,3]
    return hcat(R,t)
end
end  # module CameraCalib

"""
tests
asample = [182 535 171 537; 350 358 553 563; 1 1 1 1]
bsample = [0 888 0 888; 0 0 500 500; 1 1 1 1]


testobj = [0. 1. 2. 3. 4. 5. 6. 0. 1. 2. 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6;
           4 4 4 4 4 4 4 3 3 3 3 3 3 3 2 2 2 2 2 2 2 1 1 1 1 1 1 1 0 0 0 0 0 0 0;
           1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]

testimgpts = [[592.2633 653.09924 714.3591 776.29395 838.91156 902.527 967.32227 573.7549 635.1121 696.7237 759.04254 822.108 886.0505 951.25116 555.0277 616.6187 678.6967 741.484 804.88855 869.41974 934.88934 535.66956 597.74225 660.4037 723.53906 787.5381 852.47253 918.50275 515.9703 578.48724 641.5394 705.3553 769.73474 835.2958 901.7051;961.7239 976.7783 992.08234 1007.34906 1022.76776 1038.4393 1054.3014 1022.00385 1037.5725 1053.1384 1068.8923 1084.8145 1101.1063 1117.5796 1083.35 1099.186 1115.2794 1131.4943 1148.0342 1164.8348 1182.111 1145.628 1162.0836 1178.5444 1195.3885 1212.4789 1229.9805 1247.8092 1209.4084 1226.151 1243.2577 1260.5515 1278.3271 1296.4642 1315.0276; 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1], [674.8609 727.98016 779.95856 831.293 881.73804 931.8534 981.5053 716.21625 769.2318 821.40515 872.7004 923.50214 973.66504 1023.5292 757.9965 811.3215 863.4914 915.12305 965.9224 1016.38257 1066.3563 800.7805 854.0147 906.5378 958.2317 1009.37115 1059.8568 1109.9656 844.23456 897.7622 950.3341 1002.2766 1053.4996 1104.2455 1154.2904; 1205.071 1163.7655 1123.5414 1083.8586 1044.8091 1006.29535 967.9789 1258.4407 1216.5182 1175.4103 1135.056 1095.3274 1056.0822 1017.2831 1312.9585 1270.1844 1228.3265 1187.1783 1146.6536 1106.7433 1067.3729 1368.5806 1324.9633 1282.2985 1240.3256 1199.1078 1158.4658 1118.3892 1425.1616 1380.7924 1337.2587 1294.4718 1252.4662 1211.0784 1170.275; 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],[585.37946 656.9197 728.30743 799.52435 871.03687 942.7587 1014.948 582.7957 654.8577 726.54987 798.268 870.07837 942.32196 1014.87634 580.22986 652.6165 724.7922 796.873 869.2403 941.8298 1015.02026 577.3089 650.3446 722.9474 795.5354 868.39667 941.5566 1015.2981 574.2917 647.7883 721.1167 794.3217 867.6887 941.4382 1015.6252;862.67017 864.45166 866.06635 867.48676 868.6601 869.74786 870.7049 932.5665 934.32007 935.7993 937.31305 938.5743 939.77203 940.93 1003.2936 1004.8191 1006.36804 1007.7958 1009.2827 1010.6019 1011.9532 1074.7859 1076.3717 1077.7374 1079.3152 1080.7668 1082.409 1083.8876 1147.5334 1148.8306 1150.3627 1151.8218 1153.5355 1155.252 1157.0806; 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],[531.8516 573.4479 614.3756 654.57697 694.2811 733.538 772.48413 590.47864 631.7903 672.4354 712.45355 751.9731 791.1963 830.1773 649.4317 690.4351 730.7844 770.6563 810.1748 849.3436 888.2964 708.72363 749.5596 789.7047 829.52374 868.90356 908.0785 947.0286 768.6923 809.3671 849.4014 889.09753 928.42706 967.5378 1006.5149;1115.8102 1057.2235 999.3525 942.07184 885.32214 828.8219 772.6122 1155.7935 1096.5599 1038.2855 980.54877 923.42804 866.59845 810.0582 1196.1228 1136.3735 1077.4707 1019.3606 961.6923 904.5434 847.6189 1236.9254 1176.5455 1117.148 1058.4702 1000.3988 942.7334 885.48773 1278.3625 1217.3563 1157.3229 1098.0851 1039.4376 981.332 923.55176; 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],[491.51382 550.1794 607.4642 663.70624 719.2153 774.1711 828.6427 572.1291 630.21906 687.0408 742.8777 798.18335 853.0171 907.53864 652.99426 710.5357 766.96356 822.5869 877.6892 932.5778 987.1633 734.51434 791.6003 847.6954 903.2611 958.34814 1013.201 1067.6447 817.06616 873.74066 929.63696 985.0527 1040.1337 1094.8501 1149.2551;1094.3114 1014.279 935.49817 857.6356 780.40076 703.4269 626.9951 1150.6583 1069.8397 990.4713 912.06635 834.24915 756.7018 679.61334 1207.6022 1125.8734 1045.6892 966.58 888.2288 810.2001 732.53076 1265.403 1182.6146 1101.5314 1021.63574 942.5429 863.88257 785.67676 1324.1888 1240.353 1158.2549 1077.3331 997.43353 918.0722 839.4289; 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]]
"""
