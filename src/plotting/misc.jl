

# Function to find continuous segments of true values in a BitVector
function find_continuous_segments(mask::BitVector)
    segments = []
    current_start = nothing

    for (index, value) in enumerate(mask)
        if value
            if current_start === nothing
                current_start = index  # Start a new segment
            end
        else
            if current_start !== nothing
                push!(segments, (current_start, index - 1))  # Save the segment
                current_start = nothing  # Reset for the next segment
            end
        end
    end

    # Save the last segment if it ended with true values
    if current_start !== nothing
        push!(segments, (current_start, length(mask)))  # Save the final segment
    end
    
    return segments
end