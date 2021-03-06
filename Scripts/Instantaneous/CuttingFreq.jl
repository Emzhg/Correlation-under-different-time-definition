## Author: Patrick Chang
# Script file to investigate the impact of the cutting frequencies in the
# Fourier instantaneous estimators under the asynchronous case.
# Compared are the Malliavin-Mancino estimator and the jump robust
# adaptation by Cuchiero and Teichmann.

# Note: The Malliavin-Mancino estimator has the ability to deal with asynchrony
# while the adaptation by Cuchiero and Teichmann requires the synchronisation
# of the data beforehand.

using LinearAlgebra, Plots, LaTeXStrings, StatsBase, Intervals, JLD, ProgressMeter, Distributions, CSV
#---------------------------------------------------------------------------

cd("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-Inst")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-JR")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/GBM")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Merton Model")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Heston")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Bates2D")

#---------------------------------------------------------------------------
## Supporting functions

function rexp(n, mean)
    t = -mean .* log.(rand(n))
end

function synchronise(P, t1, t2, T, τ)
    # t1 = [0; t1]
    # t2 = [0; t2]

    t = collect(0:τ:T)
    n = length(t)
    p1 = zeros(n,1)
    p2 = zeros(n,1)
    for j in 1:n
        γ1 = maximum(filter(x-> x .<= t[j], t1))
        γ2 = maximum(filter(x-> x .<= t[j], t2))
        p1[j] = P[Int(floor(γ1)+1), 1]
        p2[j] = P[Int(floor(γ2)+1), 2]
    end
    p = [p1 p2]
    return p
end

#---------------------------------------------------------------------------
# Compute Impact of M
#---------------------------------------------------------------------------
## Asynchronous, Stochastic Volatility, No Jumps - Heston
# Setup
nsim = 21600
outlength = 1000
# M = 100
lam = 30
dt = 1
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

# Simulate and create the price and time matrices
P_Heston = Heston_CT(nsim, seed = 1, dt = nsim)
t = collect(1:1:nsim)

Random.seed!(1)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(2)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

# Synchronised using PTI
P_Heston_syn = synchronise(P_Heston[1], t1, t2, nsim, 1)

# Padding the matrix with NaNs
p1 = P_Heston[1][Int.(floor.(t1)) .+ 1, 1]
p2 = P_Heston[1][Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_Heston_asyn = [p1 p2]
t_Heston_asyn = [t1 t2]

#---------------------------------------------------------------------------
# M = 50

tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

MM_M50 = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = 50, N = N)
JR_M50 = MM_JR(P_Heston_syn, 50, outlength)


p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_M50[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_M50[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/CFM50_11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_M50[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_M50[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/CFM50_22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_M50[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_M50[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/CFM50_12.svg")

#---------------------------------------------------------------------------
# M = 15

MM_M15 = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = 15, N = N)
JR_M15 = MM_JR(P_Heston_syn, 15, outlength)


p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_M15[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_M15[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/CFM15_11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_M15[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_M15[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/CFM15_22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_M15[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_M15[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/CFM15_12.svg")

#---------------------------------------------------------------------------
# Compute Impact of N
#---------------------------------------------------------------------------
## Asynchronous, Stochastic Volatility, No Jumps - Heston
# Setup
nsim = 21600
outlength = 1000
# M = 100
lam = 30

# Simulate and create the price and time matrices
P_Heston = Heston_CT(nsim, seed = 1, dt = nsim)
t = collect(1:1:nsim)

Random.seed!(1)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(2)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

# Padding the matrix with NaNs
p1 = P_Heston[1][Int.(floor.(t1)) .+ 1, 1]
p2 = P_Heston[1][Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_Heston_asyn = [p1 p2]
t_Heston_asyn = [t1 t2]


#---------------------------------------------------------------------------
# Δt = 1s; M = 15
tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim


dt = 1
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

# Synchronised using PTI
P_syn_1s = synchronise(P_Heston[1], t1, t2, nsim, 1)

MM_1s = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = 15, N = N)
JR_1s = MM_JR(P_syn_1s, 15, outlength)


p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_1s[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_1s[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/CFN1s_11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_1s[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_1s[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/CFN1s_22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_1s[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_1s[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/CFN1s_12.svg")

#---------------------------------------------------------------------------
# Δt = 25s; M = 15

dt = 25
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

# Synchronised using PTI
P_syn_50s = synchronise(P_Heston[1], t1, t2, nsim, 25)

MM_50s = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = 15, N = N)
JR_50s = MM_JR(P_syn_50s, 15, outlength)


p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_50s[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_50s[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/CFN25s_11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_50s[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_50s[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/CFN25s_22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_50s[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_50s[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/CFN25s_12.svg")

#---------------------------------------------------------------------------
# Δt = 100; M = 15

dt = 100
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

# Synchronised using PTI
P_syn_100s = synchronise(P_Heston[1], t1, t2, nsim, 100)

MM_100s = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = 15, N = N)
JR_100s = MM_JR(P_syn_100s, 15, outlength)


p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_100s[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_100s[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/CFN50s_11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_100s[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_100s[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/CFN50s_22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_100s[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_100s[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/CFN50s_12.svg")
