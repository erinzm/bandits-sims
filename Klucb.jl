module Klucb

export KlucbState, getQuery, processAnswer!

const Arm = Int
const Reward = Int

type KlucbState
    N :: Int
    δ :: Float64
    xsum :: Array{Int,1}
    x2sum :: Array{Int,1}
    T :: Array{Int,1}
    μ :: Array{Float64,1}
    UCB :: Array{Float64,1}

    function KlucbState(N::Int, δ::Float64)
        new(N, δ,
            zeros(N), zeros(N), zeros(N),
            [Inf for _=1:N], [Inf for _=1:N])
    end
end

function getQuery(state::KlucbState)::Arm
    findmax(state.UCB)[2]
end

function processAnswer!(state::KlucbState, arm::Arm, reward::Reward)
    state.xsum[arm] += reward
    state.x2sum[arm] += reward^2
    state.T[arm] += 1

    state.μ[arm] = state.xsum[arm] / state.T[arm]
    state.UCB[arm] = computeUCB((state.μ[arm]-1)/2,
        log(2*state.T[arm].^2/state.δ)/state.T[arm])
end

function computeUCB(μ_hat, threshold, accuracy=(10^-6))
    lower = μ_hat
    upper = 1
    ucb = (lower+upper)/2
    while (upper-lower) > accuracy
        lower, upper, ucb = leftright(μ_hat, lower, upper, threshold)
    end

    return ucb
end

function leftright(μ_hat, lower, upper, threshold)
    if μ_hat*(1-μ_hat) != 0
        shit = (upper+lower)/2
        kl = μ_hat*log(μ_hat/shit) + (1-μ_hat)*log((1-μ_hat)/(1-shit))
        if kl >= threshold
            return lower, shit, (shit+lower)/2
        else
            return shit, upper, (shit+upper)/2
        end
    end

    if μ_hat == 0
        shit = (upper+lower)/2
        kl = (1-μ_hat)*log((1-μ_hat)/(1-shit))
        if kl >= threshold
            return lower, shit, (shit+lower)/2
        else
            return shit, upper, (shit+upper)/2
        end
    end

    if μ_hat == 1
        return 1, 1, 1
    end
end
end