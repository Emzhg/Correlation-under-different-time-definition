## Author: Patrick Chang
# Script file to investigate the Fourier instantaneous estimators
# under the synchronous case. Compared are the Malliavin-Mancino
# estimator and the jump robust adaptation by Cuchiero and Teichmann.

using LinearAlgebra, Plots, LaTeXStrings, StatsBase, Intervals, JLD, ProgressMeter, Distributions, CSV
#---------------------------------------------------------------------------

cd("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-Inst")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-JR")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/LRV")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/GBM")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Merton Model")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Heston")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Bates2D")

#---------------------------------------------------------------------------
## Synchronous, Constant Volatility, No Jumps - GBM
# Setup
nsim = 28800
outlength = 1000
M = 100

mu = [0.01/28800, 0.01/28800]
sigma = [0.1/28800 sqrt(0.1/28800)*0.35*sqrt(0.2/28800);
        sqrt(0.1/28800)*0.35*sqrt(0.2/28800) 0.2/28800]

# Sim
P_GBM = GBM(nsim, mu, sigma)
t = collect(1:1:nsim)

# Compute
MM_GBM = MM_inst(P_GBM, [t t], outlength, M = M)
JR_GBM = MM_JR(P_GBM, M, outlength)

# Plot
tt = collect(0:1/outlength:1-(1/outlength))

p1 = plot(tt, MM_GBM[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0.05, 0.35), dpi = 300)
plot!(p1, tt, JR_GBM[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p1, [0.1], label = L"\textrm{True } \Sigma^{11}(t)", color = :black, line=(1, [:solid]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/SynGBM11.svg")


p2 = plot(tt, MM_GBM[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0.17, 0.45), dpi = 300)
plot!(p2, tt, JR_GBM[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p2, [0.2], label = L"\textrm{True } \Sigma^{22}(t)", color = :black, line=(1, [:solid]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/SynGBM22.svg")


p3 = plot(tt, MM_GBM[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0.02, 0.08), dpi = 300)
plot!(p3, tt, JR_GBM[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p3, [sqrt(0.1)*0.35*sqrt(0.2)], label = L"\textrm{True } \Sigma^{12}(t)", color = :black, line=(1, [:solid]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/SynGBM12.svg")


#---------------------------------------------------------------------------
## Synchronous, Constant Volatility, With Jumps - Merton
# Setup
nsim = 28800
outlength = 1000
M = 100

mu = [0.01/28800, 0.01/28800]
sigma = [0.1/28800 sqrt(0.1/28800)*0.35*sqrt(0.2/28800);
        sqrt(0.1/28800)*0.35*sqrt(0.2/28800) 0.2/28800]

a = [-0.005; -0.003]
b= [0.015; 0.02]
lambda = [100/28800, 100/28800]

# Sim
P_Mert = Merton(nsim, mu, sigma, lambda, a, b)
t = collect(1:1:nsim)
# t = (t .- minimum(t)) .* (1 / (maximum(t) - minimum(t)))

# Compute
MM_Mert = MM_inst(P_Mert, [t t], outlength, M = M)
JR_Mert = MM_JR(P_Mert, M, outlength)

# Plot
tt = collect(0:1/outlength:1-(1/outlength))

p1 = plot(tt, MM_Mert[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0.05, 0.35), dpi = 300)
plot!(p1, tt, JR_Mert[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p1, [0.1], label = L"\textrm{True } \Sigma^{11}(t)", color = :black, line=(1, [:solid]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/SynMert11.svg")


p2 = plot(tt, MM_Mert[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0.17, 0.45), dpi = 300)
plot!(p2, tt, JR_Mert[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p2, [0.2], label = L"\textrm{True } \Sigma^{22}(t)", color = :black, line=(1, [:solid]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/SynMert22.svg")


p3 = plot(tt, MM_Mert[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]), ylims = (0.02, 0.08), dpi = 300)
plot!(p3, tt, JR_Mert[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
hline!(p3, [sqrt(0.1)*0.35*sqrt(0.2)], label = L"\textrm{True } \Sigma^{12}(t)", color = :black, line=(1, [:solid]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/SynMert12.svg")

#---------------------------------------------------------------------------
## Synchronous, Stochastic Volatility, No Jumps - Heston
nsim = 28800
outlength = 1000
M = 100

P_Heston = Heston_CT(nsim, seed = 1, dt = nsim)
t = collect(1:1:nsim)

MM_Heston = MM_inst(P_Heston[1], [t t], outlength, M = M)
JR_Heston = MM_JR(P_Heston[1], M, outlength)
# LRV_Heston = LRV(P_Heston[1], 500, 500)


tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_Heston[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_Heston[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/SynHeston11.svg")


p2 = plot(t, P_Heston[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_Heston[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_Heston[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/SynHeston22.svg")


p3 = plot(t, P_Heston[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_Heston[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_Heston[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/SynHeston12.svg")

# tt = collect(0:1/500:1-(1/500))
# ttt = collect(0:1/outlength:1-(1/outlength))
# p4 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(1, [:solid]), dpi = 300)
# plot!(p4, tt, LRV_Heston[1], label = L"\textrm{LRV } \hat{\Sigma}^{11}_{M}(t)", color = :orange, line=(0.3, [:solid]), alpha = 0.5)
# plot!(p4, ttt, MM_Heston[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
# plot!(p4, ttt, JR_Heston[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
# xlabel!(p4, L"\textrm{Time}")
# ylabel!(p4, L"\textrm{Instantaneous covariance of } \Sigma^{11}(t)")

#---------------------------------------------------------------------------
## Synchronous, Stochastic Volatility, With Jumps - Bates
nsim = 28800
outlength = 1000
M = 100

P_Bates = Bates_CT(nsim, seed = 1, dt = nsim)
t = collect(1:1:nsim)

MM_Bates = MM_inst(P_Bates[1], [t t], outlength, M = M)
JR_Bates = MM_JR(P_Bates[1], M, outlength)


tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

p1 = plot(t, P_Bates[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_Bates[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p1, tt, JR_Bates[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

# savefig(p1, "Plots/Instantaneous/SynBates11.svg")


p2 = plot(t, P_Bates[3], label = L"\textrm{True } \Sigma^{22}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p2, tt, MM_Bates[2], label = L"\textrm{Estimated MM } \hat{\Sigma}^{22}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p2, tt, JR_Bates[2], label = L"\textrm{Estimated CT } \hat{\Sigma}^{22}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p2, L"\textrm{Time}")
ylabel!(p2, L"\textrm{Instantaneous variance of } \Sigma^{22}(t)")

# savefig(p2, "Plots/Instantaneous/SynBates22.svg")


p3 = plot(t, P_Bates[4], label = L"\textrm{True } \Sigma^{12}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p3, tt, MM_Bates[3], label = L"\textrm{Estimated MM } \hat{\Sigma}^{12}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
plot!(p3, tt, JR_Bates[3], label = L"\textrm{Estimated CT } \hat{\Sigma}^{12}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p3, L"\textrm{Time}")
ylabel!(p3, L"\textrm{Instantaneous covariance of } \Sigma^{12}(t)")

# savefig(p3, "Plots/Instantaneous/SynBates12.svg")
