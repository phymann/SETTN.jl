let
    @show @__MODULE__
    function jwf(x)
        return x^2
    end

    @show methods(jwf)
    str = "jwf"
    # @show getfield(Main, Symbol(str))

    # getfield(Main, Symbol(str))(2)
end
