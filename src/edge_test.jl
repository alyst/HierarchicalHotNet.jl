@inline isweaker(a::Number, b::Number; rev::Bool=false) =
    (!rev && (a < b)) || (rev && (a > b))
@inline isstronger(a::Number, b::Number; rev::Bool=false) =
    (!rev && (a > b)) || (rev && (a < b))

@inline isvalidedge(w::Number;
            skipval::Union{Number, Nothing} = zero(w),
            threshold::Union{Number, Nothing} = nothing,
            rev::Bool = false) =
    (isnothing(skipval) || (w != skipval)) &&
    (isnothing(threshold) || !isweaker(w, threshold, rev=rev))

@inline defaultweight(::Type{T};
              skipval::Union{Number, Nothing} = zero(T),
              rev::Bool = false) where T =
    isnothing(skipval) ? (rev ? typemax(T) : typemin(T)) : skipval

struct EdgeTest{W, T, S, R}
    threshold::T

    function EdgeTest{W}(; skipval::Union{Nothing, Number} = zero(W),
        threshold::Union{Nothing, Number} = nothing,
        rev::Bool = false) where W <: Number
        used_threshold = isnothing(threshold) ? nothing : convert(W, threshold)
        used_skipval = isnothing(skipval) ? nothing : convert(W, skipval)
        new{W, typeof(used_threshold), used_skipval, rev}(used_threshold)
    end
end

skipval(::Type{<:EdgeTest{<:Any, <:Any, S}}) where S = S
skipval(test::EdgeTest) = skipval(typeof(test))

isreverse(::Type{<:EdgeTest{<:Any, <:Any, <:Any, R}}) where R = R
isreverse(test::EdgeTest) = isreverse(typeof(test))

@inline isvalidedge(w::Number, test::EdgeTest) =
    isvalidedge(w, skipval=skipval(test), threshold=test.threshold, rev=isreverse(test))

defaultweight(test::EdgeTest{W}) where W =
    defaultweight(W, skipval=skipval(test), rev=isreverse(test))
