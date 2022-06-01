
function fun_normalize(p,S)

# Calculate normalization factor
f = sqrt(transpose(p)*S*p)
# Normalize such that transpose(p)*S*p = 1
p_out = p/f
# Multiply by -1 if int(p) < 0
    if (real(transpose(1 .+0*p)*S*real(p))) < 0
        p_out= -p_out
        f = -f
    end

    return p_out,f
end
