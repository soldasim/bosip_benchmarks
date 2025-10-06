
function check_alias(alias::Function)
    @assert string(alias) == get_name(alias())
end
