#< CUSTOM PEAK FINDER
function findmaxpeaks(x; 
                    height::Union{Nothing,<:Real,NTuple{2,<:Real}}=nothing,
                    distance::Union{Nothing,Int}=nothing)#,
                    #find_maxima::Bool=true) #where {T<:Real}
    midpts = Vector{Int}(undef, 0)
    i = 2
    imax = length(x)

    while i < imax
        if x[i-1] < x[i] #|| (!find_maxima && x[i-1] > x[i])
            iahead = i + 1
            while (iahead < imax) && (x[iahead] == x[i])
                iahead += 1
            end

            if x[iahead] < x[i] #|| (!find_maxima && x[iahead] > x[i])
                push!(midpts, (i + iahead - 1) ÷ 2)
                i = iahead
            end
        end
        i += 1
    end 

    #* Filter by height if needed
    if !isnothing(height)
        hmin, hmax = height isa Number ? (height, nothing) : height
        keepheight = (hmin === nothing || x[midpts] .>= hmin) .& (hmax === nothing || x[midpts] .<= hmax)
        midpts = midpts[keepheight]
    end

    #* Filter by distance if needed
    if !isnothing(distance)
        # priority = find_maxima ? x[midpts] : -x[midpts] # Use negative values for minima
        priority = x[midpts] 
        keep = selectbypeakdistance(midpts, priority, distance)
        midpts = midpts[keep]
    end

    midpts
    # extrema_indices = midpts
    # # extrema_heights = [0.0] # initial value of 0.0 so that getDif returns 0.0 if no peaks are found
    # # append!(extrema_heights, x[extrema_indices]) 
    # extrema_heights = x[extrema_indices]

    # extrema_indices, extrema_heights
end




"""
    findextrema(x::Vector{Float64}; min_height::Float64=0.0)

Find and filter the maxima and minima of a data vector based on height prominence.

# Arguments
- `x::Vector{Float64}`: The data vector
- `min_prominence::Float64`: Minimum height prominence for maxima and corresponding minima

# Returns
- `maxima_indices::Vector{Int}`: Indices of the filtered maxima
- `minima_indices::Vector{Int}`: Indices of the filtered minima

"""
function findextrema(x; min_prominence::Float64=0.0)
    maxima_indices = Vector{Int}(undef, 0)
    minima_indices = Vector{Int}(undef, 0)
    
    # Initial loop to identify potential maxima and minima
    for i in 2:(length(x) - 1)
        if x[i] > x[i-1] && x[i] > x[i+1]
            push!(maxima_indices, i)
        elseif x[i] < x[i-1] && x[i] < x[i+1]
            push!(minima_indices, i)
        end
    end
    
    if isempty(maxima_indices) || isempty(minima_indices)
        return maxima_indices, minima_indices
    end

    # Filter maxima and corresponding minima based on height prominence
    if min_prominence > 0.0
        filter_extrema!(x, maxima_indices, minima_indices, min_prominence)
    end
    
    # # Get values of the filtered maxima and minima
    # maxima_values = @view x[maxima_indices]
    # minima_values = @view x[minima_indices]
    
    return maxima_indices, minima_indices
end


"""
    filter_extrema!(x::Vector{Float64}, maxima_indices::Vector{Int}, minima_indices::Vector{Int}, min_prominence::Float64)

Filter out extrema based on prominence.

# Arguments
- `x::Vector{Float64}`: The original data vector
- `maxima_indices::Vector{Int}`: Indices of maxima
- `minima_indices::Vector{Int}`: Indices of minima
- `min_height::Float64`: Minimum prominence for an extremum to be considered
"""
function filter_extrema!(x, maxima_indices, minima_indices, min_prominence::Float64)
    # Calculate the prominence for each extremum
    max_prominences, min_prominences = calculate_prominence(x, maxima_indices, minima_indices)

    # Filter maxima by prominence
    keep_max_by_prominence = max_prominences .>= min_prominence
    deleteat!(maxima_indices, .!keep_max_by_prominence)

    # Filter minima by prominence
    keep_min_by_prominence = min_prominences .>= min_prominence
    deleteat!(minima_indices, .!keep_min_by_prominence)
end



"""
    calculate_prominence(x::Vector{Float64}, maxima_indices::Vector{Int}, minima_indices::Vector{Int})

Calculate the prominence of each maxima relative to its preceding minima, and vice versa.

# Arguments
- `x::Vector{Float64}`: The original data vector
- `maxima_indices::Vector{Int}`: Indices of maxima
- `minima_indices::Vector{Int}`: Indices of minima

# Returns
- `max_prominences::Vector{Float64}`: Prominence of each maxima
- `min_prominences::Vector{Float64}`: Prominence of each minima
"""
function calculate_prominence(x, maxima_indices, minima_indices)
    n_max = length(maxima_indices)
    n_min = length(minima_indices)
    max_prominences = Vector{Float64}(undef, n_max)
    min_prominences = Vector{Float64}(undef, n_min)

    # Initialize the first preceding minima or maxima
    prev_min_idx = first(maxima_indices) < first(minima_indices) ? 1 : first(minima_indices)
    prev_max_idx = first(minima_indices) < first(maxima_indices) ? 1 : first(maxima_indices)
    
    # Calculate the prominence of each maxima
    for (i, max_idx) in enumerate(maxima_indices)
        found_idx = findlast(min_idx -> min_idx < max_idx, minima_indices)
        prev_min_idx = isnothing(found_idx) ? i : minima_indices[found_idx]
        max_prominences[i] = x[max_idx] - x[prev_min_idx]
    end
    
    # Calculate the prominence of each minima
    for (i, min_idx) in enumerate(minima_indices)
        found_idx = findlast(max_idx -> max_idx < min_idx, maxima_indices)
        prev_max_idx = isnothing(found_idx) ? i : maxima_indices[found_idx]
        min_prominences[i] = x[prev_max_idx] - x[min_idx]
    end
    
    return max_prominences, min_prominences
end


"""
    selectbypeakdistance(peak_indices, priority, min_distance)

Filters out peaks that are too close to higher-priority peaks.

# Arguments
- `peak_indices::Vector{Int}`: Indices of the peaks
- `priority::Vector`: Priority levels for each peak, usually based on height
- `min_distance::Int`: The minimum distance required between peaks

# Returns
- `keep::BitVector`: A boolean vector indicating which peaks to keep
"""
function selectbypeakdistance(peak_indices, priority, min_distance::Int)
    num_peaks = length(peak_indices)
    keep = trues(num_peaks)  # Initialize with all true

    # Sort peaks by priority (usually height), in descending order
    sorted_by_priority = sortperm(priority, rev=true)

    # Loop through each peak, starting with the highest priority
    for i in num_peaks:-1:1
        current_peak_pos = sorted_by_priority[i]

        # Skip if this peak is already marked for removal
        iszero(keep[current_peak_pos]) && continue

        # Check lower-indexed peaks
        lower_idx = current_peak_pos - 1
        while (lower_idx >= 1) && (abs(peak_indices[current_peak_pos] - peak_indices[lower_idx]) < min_distance)
            keep[lower_idx] = false
            lower_idx -= 1
        end

        # Check higher-indexed peaks
        higher_idx = current_peak_pos + 1
        while (higher_idx <= num_peaks) && (abs(peak_indices[higher_idx] - peak_indices[current_peak_pos]) < min_distance)
            keep[higher_idx] = false
            higher_idx += 1
        end
    end
    return keep
end
