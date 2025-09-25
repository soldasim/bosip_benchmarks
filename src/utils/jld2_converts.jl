
function Base.convert(::Type{BosipCallback}, reconstructed::JLD2.ReconstructedMutable{:SaveCallback, (:dir, :filename, :run_idx), Tuple{String, String, Any}})
    return SaveCallback(;
        dir = reconstructed.dir,
        filename = reconstructed.filename,
        run_idx = reconstructed.run_idx,
        continued = false,
    )
end
