function bits_signed_exp(x::Float32)
    # converts the exponent bits from an unsigned integer to a signed integer
    all_bits = bitstring(x)  # unsigned exp

    bias = 127
    k = parse(Int,"0b"*all_bits[2:9])-bias    # signed exponent

    # subnormal numbers, i.e. x=0 means all exponent bits are 0
    if k == -bias   # i.e. x was 0
        s_exp_bits = "00000000"
    else
        if k < 0
            s0 = "1"
        else
            s0 = "0"
        end

        s_exp_bits = s0*bitstring(abs(k))[end-6:end]    # the remaining 7 bits of the exponent
    end
    return all_bits[1]*s_exp_bits*all_bits[10:end]
end

function bits2bitarray(x::Array{Float32,1})
    B = BitArray{2}(undef,32,length(x))
    for (xin,xi) in enumerate(x)
        for (ib,b) in enumerate(bits_signed_exp(xi))
            if b == "1"
                B[ib,xin] = true
            else
                B[ib,xin] = false
            end
        end
    end
    return B
end
