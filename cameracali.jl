function find_homography(src::Array{T1,2}, dst::Array{T2,2}) where {T1<:Real, T2<:Real}
    if size(src) != size(dst)
        error("src and dst must have same dimensions")
    end
end

function normalizepoints(...)

end
