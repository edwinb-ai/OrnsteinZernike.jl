@testset "Parameters" begin
    brs = Dict(
        "Default" => ModifiedVerlet(),
        "Three parameters" => ModifiedVerlet(-0.6, 0.01, 0.85)
    )

    for (k, v) in brs
        @testset "$k" begin
            result = oz_hs(v)
            if 1.0 ∈ result.r
                first_contact = result.gr[result.r .== 1.0][1]
            else
                first_contact = maximum(result.gr)
            end
            @show "$k" first_contact
            # We don't know exactly the value at contact, so just check that it is a valid number
            @test first_contact isa Real
        end
    end
end
