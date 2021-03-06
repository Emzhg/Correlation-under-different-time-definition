## Author: Patrick Chang
# Script file to investigate the Fourier instantaneous estimators
# under the asynchronous case. Compared are the Malliavin-Mancino
# estimator and the jump robust adaptation by Cuchiero and Teichmann.

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

function synchronise(P, t1, t2, T, τ = 1)
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
## Asynchronous, Constant Volatility, No Jumps - GBM
# Setup
nsim = 28800
outlength = 1000
M = 100
lam = 30
dt = 1
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

mu = [0.01/28800, 0.01/28800]
sigma = [0.1/28800 sqrt(0.1/28800)*0.35*sqrt(0.2/28800);
        sqrt(0.1/28800)*0.35*sqrt(0.2/28800) 0.2/28800]

# Simulate and create the price and time matrices
P_GBM = GBM(nsim, mu, sigma)
t = collect(1:1:nsim)

Random.seed!(1)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(2)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

P_GBM_syn = synchronise(P_GBM, t1, t2, nsim)

p1 = P_GBM[Int.(floor.(t1)) .+ 1, 1]
p2 = P_GBM[Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_GBM_asyn = [p1 p2]
t_GBM_asyn = [t1 t2]


# Compute
MM_GBM = MM_inst(P_GBM_asyn, t_GBM_asyn, outlength, M = M, N = N)
JR_GBM = MM_JR(P_GBM_syn, M, outlength)

# Plot
tt = collect(0:1/outlength:1-(1/outlength))

p1 = plot(tt, MM_GBM[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n_1,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0, 0.45), dpi = 300)
plot!(p1, tt, JR_GBM[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p1, [0.1], label = L"\textrm{True } \Sigma^{11}(t)", color = :black, line=(1, [:solid]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/AsynGBM11.svg")


p2 = plot(tt, MM_GBM[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n_2,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0, 1.6), dpi = 300)
plot!(p2, tt, JR_GBM[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p2, [0.2], label = L"\textrm{True } \Sigma^{22}(t)", color = :black, line=(1, [:solid]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/AsynGBM22.svg")


p3 = plot(tt, MM_GBM[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n_1,n_2,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (-0.04, 0.053), dpi = 300)
plot!(p3, tt, JR_GBM[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p3, [sqrt(0.1)*0.35*sqrt(0.2)], label = L"\textrm{True } \Sigma^{12}(t)", color = :black, line=(1, [:solid]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/AsynGBM12.svg")


# p4 = plot(tt, MM_GBM[3] ./ sqrt.(MM_GBM[1] .* MM_GBM[2]), label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
# plot!(p4, tt, JR_GBM[3] ./ sqrt.(JR_GBM[1] .* JR_GBM[2]), label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
# hline!(p4, [0.35], label = L"\textrm{True } {\rho}^{12}(t)", color = :black, line=(1, [:solid]))
# xlabel!(p4, L"\textrm{Time}")
# ylabel!(p4, L"\textrm{Instantaneous correlation } \rho^{12}(t)")

#---------------------------------------------------------------------------
## Asynchronous, Constant Volatility, With Jumps - Merton
# Setup
nsim = 28800
outlength = 1000
M = 100
lam = 30
dt = 1
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

mu = [0.01/28800, 0.01/28800]
sigma = [0.1/28800 sqrt(0.1/28800)*0.35*sqrt(0.2/28800);
        sqrt(0.1/28800)*0.35*sqrt(0.2/28800) 0.2/28800]

a = [-0.005; -0.003]
b= [0.015; 0.02]
lambda = [100/28800, 100/28800]

# Simulate and create the price and time matrices
P_Mert = Merton(nsim, mu, sigma, lambda, a, b)
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
P_Mert_syn = synchronise(P_Mert, t1, t2, nsim)

# Padding the matrix with NaNs
p1 = P_Mert[Int.(floor.(t1)) .+ 1, 1]
p2 = P_Mert[Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_Mert_asyn = [p1 p2]
t_Mert_asyn = [t1 t2]


# Compute
MM_Mert = MM_inst(P_Mert_asyn, t_Mert_asyn, outlength, M = M, N = N)
JR_Mert = MM_JR(P_Mert_syn, M, outlength)

# Plot
tt = collect(0:1/outlength:1-(1/outlength))

p1 = plot(tt, MM_Mert[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n_1,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0, 0.45), dpi = 300)
plot!(p1, tt, JR_Mert[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p1, [0.1], label = L"\textrm{True } \Sigma^{11}(t)", color = :black, line=(1, [:solid]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/AsynMert11.svg")


p2 = plot(tt, MM_Mert[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n_2,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0, 1.6), dpi = 300)
plot!(p2, tt, JR_Mert[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p2, [0.2], label = L"\textrm{True } \Sigma^{22}(t)", color = :black, line=(1, [:solid]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/AsynMert22.svg")


p3 = plot(tt, MM_Mert[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n_1,n_2,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (-0.04, 0.053), dpi = 300)
plot!(p3, tt, JR_Mert[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p3, [sqrt(0.1)*0.35*sqrt(0.2)], label = L"\textrm{True } \Sigma^{12}(t)", color = :black, line=(1, [:solid]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/AsynMert12.svg")


# p4 = plot(tt, MM_Mert[3] ./ sqrt.(MM_Mert[1] .* MM_Mert[2]), label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
# plot!(p4, tt, JR_Mert[3] ./ sqrt.(JR_Mert[1] .* JR_Mert[2]), label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
# hline!(p4, [0.35], label = L"\textrm{True } {\rho}^{12}(t)", color = :black, line=(1, [:solid]))
# xlabel!(p4, L"\textrm{Time}")
# ylabel!(p4, L"\textrm{Instantaneous correlation } \rho^{12}(t)")

#---------------------------------------------------------------------------
## Asynchronous, Stochastic Volatility, No Jumps - Heston
nsim = 28800
outlength = 1000
M = 100
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
P_Heston_syn = synchronise(P_Heston[1], t1, t2, nsim)

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


# Compute
MM_Heston = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = M, N = N)
JR_Heston = MM_JR(P_Heston_syn, M, outlength)

tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_Heston[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n_1,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_Heston[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/AsynHeston11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_Heston[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n_2,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_Heston[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/AsynHeston22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_Heston[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n_1,n_2,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_Heston[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/AsynHeston12.svg")


# p4 = plot(t, P_Heston[4] ./ sqrt.(P_Heston[2] .* P_Heston[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]))
# plot!(p4, tt, MM_Heston[3] ./ sqrt.(MM_Heston[1] .* MM_Heston[2]), label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
# plot!(p4, tt, JR_Heston[3] ./ sqrt.(JR_Heston[1] .* JR_Heston[2]), label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
# xlabel!(p4, L"\textrm{Time}")
# ylabel!(p4, L"\textrm{Instantaneous correlation } \rho^{12}(t)")

#---------------------------------------------------------------------------
## Asynchronous, Stochastic Volatility, With Jumps - Bates
nsim = 28800
outlength = 1000
M = 100
lam = 30
dt = 1
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

# Simulate and create the price and time matrices
P_Bates = Bates_CT(nsim, seed = 1, dt = nsim)
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
P_Bates_syn = synchronise(P_Bates[1], t1, t2, nsim)

# Padding the matrix with NaNs
p1 = P_Bates[1][Int.(floor.(t1)) .+ 1, 1]
p2 = P_Bates[1][Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_Bates_asyn = [p1 p2]
t_Bates_asyn = [t1 t2]


# Compute
MM_Bates = MM_inst(P_Bates_asyn, t_Bates_asyn, outlength, M = M, N = N)
JR_Bates = MM_JR(P_Bates_syn, M, outlength)

tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

p1 = plot(t, P_Bates[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_Bates[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n_1,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_Bates[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/AsynBates11.svg")


p2 = plot(t, P_Bates[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_Bates[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n_2,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_Bates[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/AsynBates22.svg")


p3 = plot(t, P_Bates[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_Bates[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n_1,n_2,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_Bates[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/AsynBates12.svg")

# p4 = plot(t, P_Bates[4] ./ sqrt.(P_Bates[2] .* P_Bates[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]))
# plot!(p4, tt, MM_Bates[3] ./ sqrt.(MM_Bates[1] .* MM_Bates[2]), label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
# plot!(p4, tt, JR_Bates[3] ./ sqrt.(JR_Bates[1] .* JR_Bates[2]), label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
# xlabel!(p4, L"\textrm{Time}")
# ylabel!(p4, L"\textrm{Instantaneous correlation } \rho^{12}(t)")
