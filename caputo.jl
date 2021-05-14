using SpecialFunctions
using RHEOS
using NLopt; using BenchmarkTools;using DSP

################## L1 operator ##################################

function a(j, α)
    return (j+1)^(1-α) - j^(1-α)
end

function L1(f, t, Δt, α)

    N = size(t,1)
    l1 = zeros(Float64, N)
    g = gamma(2-α)
    
    for k in 1:N-1
        l1_k = f[k+1]
        for j in 1:k-1 # Note Julia iterates beg and end of beg:end
            l1_k = l1_k - (a(k-j-1, α) - a(k-j, α))*f[j+1]
        end
        l1_k = (l1_k - f[1]*a(k-1, α))*(Δt^(-α))/g
        l1[k+1] = l1_k
    end

    return l1
end

##################### L1-2 Operator ############################

function b(j, α)
    # j >= 0
    return (1/(2-α))*((j+1)^(2-α) - j^(2-α)) - 0.5*((j+1)^(1-α) + j^(1-α))
end

function δt(f, k, Δt)
    # the paper implementation used δₜf_{k-1/2}  
    # modified a bit for code purposes
    return (f[k]-f[k-1])/Δt
end

function δ2t(f, k, Δt)
    return (δt(f,k+1,Δt)-δt(f,k,Δt))/Δt
end

function c(j, α, k)

    if k == 1 && j == 0
        return 1
    end

    c = 0
    if j == 0
        c = a(0, α)  + b(0, α)

    elseif j >= 1 && j <= (k-2)
        c = a(j, α) + b(j, α) - b(j-1, α)

    elseif j == k-1
        c = a(j, α)  - b(j-1, α)

    end
    return c

end

function L12(f, t, Δt, α)

    N = size(t,1)
    l12 = zeros(Float64, N)
    g = gamma(2-α)

    for k in 1:N-1
        l12_k = c(0,α,k)*f[k+1]
    
        for j in 1:k-1 # Note Julia iterates beg and end of beg:end
            l12_k = l12_k - (c(k-j-1, α, k) - c(k-j, α, k))*f[j+1]
        end
        l12_k = (l12_k - f[1]*c(k-1, α, k))*(Δt^(-α))/g
        l12[k+1] = l12_k
    end

    return l12
    
end




function fractionalDiff(sigderiv, sigdoublederiv, t, dt, α)

    I1 = sigderiv[1]*t.^(1-α)
    I2 = conv(t.^(1-α), sigdoublederiv)[1:length(t)]*dt
    return (I1 + I2)/gamma(2 - α)

end
