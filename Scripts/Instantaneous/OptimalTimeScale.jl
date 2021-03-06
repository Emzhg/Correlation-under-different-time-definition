## Author: Patrick Chang
# Script file to see if there is an optimal time-scale to reconstruct
# the instantaneous estimates. The method I propose is to perform some
# quick EDA using the integrated correlation to determine a time-scale
# where the effect of the Epps effect has been removed. Then we look
# for frequencies M ≦ N/2, where N is the time-scale where the Epps effect
# is removed. The choice of M is still unclear, but because of the quick method
# I have, we can perform some quick EDA to find an M that seems to make sense.

# Note: The experiment will use a diffusion process with a
# deterministic correlation component.

# Note: Only the MM estimator will be explored, as the CT
# is not good under asynchrony.

using LinearAlgebra, Plots, LaTeXStrings, StatsBase, Intervals, JLD, ProgressMeter, Distributions, CSV
#---------------------------------------------------------------------------

cd("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-Inst")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-JR")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Correlation Estimators/Fejer/NUFFTcorrFK-FGG")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/SDEs/Heston")

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

function deterministicCor(t)
    return sin(t*pi)
end

# vol = 2x1 vector of constant volatilities
function deterministicDiffusion(n, vol; kwargs...)
    # n - simlulation length
    # mu - vector input of the drift component
    # sigma - covariance matrix of the stocks
    # startprice - starting price for the assets

    # all inputs must have appropriate dimensions

    k = 2

    kwargs = Dict(kwargs)

    if haskey(kwargs, :startprice)
        startprice = kwargs[:startprice]
    else
        startprice = fill(100.0, (k,1))
    end

    if haskey(kwargs, :seed)
        seed = kwargs[:seed]
    else
        seed = 1
    end

    if haskey(kwargs, :dt)
        dt = kwargs[:dt]
    else
        dt = 1
    end

    P = zeros(n, k)
    P[1,:] = startprice

    Random.seed!(seed)
    Z = randn(k, n-1)

    for i in 2:n
        ρ = deterministicCor(i/n)
        sigma = [vol[1] sqrt(vol[1]*vol[2])*ρ; sqrt(vol[1]*vol[2])*ρ vol[2]]
        sigma = reshape(sigma, k, k)
        sigma2 = reshape(diag(sigma), k, 1)
        # A = cholesky(sigma).L
        A = [1 0; ρ sqrt(1-ρ^2)]
        b = -sigma2./2

        z = Z[:,i-1]
        X = b./dt + (A * (z.*vol))./sqrt(dt)
        P[i,:] = P[i-1,:] .* exp.(X)
    end
    return P
end

#---------------------------------------------------------------------------
# Experiment setup - Synchronous case, base-line for deterministic correlation
nsim = 28800
outlength = 1000

t = collect(1:1:nsim)/nsim
plot(t, deterministicCor.(t))

P_syn = deterministicDiffusion(nsim, [0.01; 0.02], dt = nsim)
t_syn = collect(1:1:nsim)

M_syn = collect(1:1:100)

MM_corr_syn = zeros(length(M_syn), outlength)

for i in 1:length(M_syn)
    MM_res = MM_inst(P_syn, [t_syn t_syn], outlength, M = M_syn[i])

    MM_corr_syn[i,:] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
end

t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))

MM_syn = @animate for i in 1:length(M_syn)
    now = M_syn[i]
    p1 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
    plot!(p1,tt, MM_corr_syn[i,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$M = $now\$"))
end

gif(MM_syn, "Plots/GIFs/OT_MM_syn.gif", fps = 10)


p2 = surface(M_syn, tt, MM_corr_syn', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p3 = contour(tt, M_syn, MM_corr_syn, fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))
p4 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :black, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
plot!(p4,tt, MM_corr_syn[20,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))


# savefig(p2, "Plots/Instantaneous/OT_surf_syn.svg")
# savefig(p3, "Plots/Instantaneous/OT_HM_syn.png")
# savefig(p4, "Plots/Instantaneous/OT_opt_syn.svg")

#---------------------------------------------------------------------------
# Asynchronous case, 1/λ = 10

lam = 10

Random.seed!(1)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(2)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

# Padding the matrix with NaNs
p1 = P_syn[Int.(floor.(t1)) .+ 1, 1]
p2 = P_syn[Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_asyn_lam10 = [p1 p2]
t_asyn_lam10 = [t1 t2]

# Quick EDA to determine optimal time-scale

Δt = collect(1:1:300)
res = zeros(length(Δt), 1)
for i in 1:length(Δt)
    NN = Int.(floor.(((nsim ./Δt[i]).-1.0) ./ 2))
    res[i] = NUFFTcorrDKFGG(P_asyn_lam10, t_asyn_lam10, N = NN)[1][1,2]
    # res[i] = NUFFTcorrFKFGG(P_asyn_lam10, t_asyn_lam10, N = NN)[1][1,2]
end

p1 = plot(Δt, res, legend = :bottomright,label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}", color = :blue, line=(1, [:solid]), dpi = 300)
vline!(p1, [60], label = L"\textrm{Choice of } \Delta t")
xlabel!(p1, L"\Delta t\textrm{[sec]}")
ylabel!(p1, L"\textrm{Integrated correlation } \rho_{\Delta t}^{12}")

# savefig(p1, "Plots/Instantaneous/OT_EDA_lam10.svg")

# Optimal time-scale for 1/λ = 10 is Δt = 60 [sec]
dt = 60
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

M_lam10 = collect(1:1:Int(floor(N/2)))

MM_corr_lam10 = zeros(length(M_lam10), outlength)

for i in 1:length(M_lam10)
    MM_res = MM_inst(P_asyn_lam10, t_asyn_lam10, outlength, M = M_lam10[i], N = N)

    MM_corr_lam10[i,:] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
end



t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))

MM_asyn_lam10 = @animate for i in 1:length(M_lam10)
    now = M_lam10[i]
    p1 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
    plot!(p1,tt, MM_corr_lam10[i,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$M = $now\$"))
end

gif(MM_asyn_lam10, "Plots/GIFs/OT_MM_asyn_lam10.gif", fps = 10)


p2 = surface(M_lam10, tt, MM_corr_lam10', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p3 = contour(tt, M_lam10, MM_corr_lam10, fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))
p4 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :black, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
plot!(p4,tt, MM_corr_lam10[15,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))


# savefig(p2, "Plots/Instantaneous/OT_surf_lam10.svg")
# savefig(p3, "Plots/Instantaneous/OT_HM_lam10.png")
# savefig(p4, "Plots/Instantaneous/OT_opt_lam10.svg")


#---------------------------------------------------------------------------
# Asynchronous case, 1/λ = 20

lam = 20

Random.seed!(3)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(4)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

# Padding the matrix with NaNs
p1 = P_syn[Int.(floor.(t1)) .+ 1, 1]
p2 = P_syn[Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_asyn_lam20 = [p1 p2]
t_asyn_lam20 = [t1 t2]

# Quick EDA to determine optimal time-scale

Δt = collect(1:1:300)
res = zeros(length(Δt), 1)
for i in 1:length(Δt)
    NN = Int.(floor.(((nsim ./Δt[i]).-1.0) ./ 2))
    res[i] = NUFFTcorrDKFGG(P_asyn_lam20, t_asyn_lam20, N = NN)[1][1,2]
    # res[i] = NUFFTcorrFKFGG(P_asyn_lam20, t_asyn_lam20, N = NN)[1][1,2]
end

p1 = plot(Δt, res, legend = :bottomright, label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}", color = :blue, line=(1, [:solid]), dpi = 300)
vline!(p1, [100], label = L"\textrm{Choice of } \Delta t")
xlabel!(p1, L"\Delta t\textrm{[sec]}")
ylabel!(p1, L"\textrm{Integrated correlation } \rho_{\Delta t}^{12}")

# savefig(p1, "Plots/Instantaneous/OT_EDA_lam20.svg")

# Optimal time-scale for 1/λ = 10 is Δt = 100 [sec]
dt = 100
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

M_lam20 = collect(1:1:Int(floor(N/2)))

MM_corr_lam20 = zeros(length(M_lam20), outlength)

for i in 1:length(M_lam20)
    MM_res = MM_inst(P_asyn_lam20, t_asyn_lam20, outlength, M = M_lam20[i], N = N)

    MM_corr_lam20[i,:] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
end



t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))

MM_asyn_lam20 = @animate for i in 1:length(M_lam20)
    now = M_lam20[i]
    p1 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
    plot!(p1,tt, MM_corr_lam20[i,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$M = $now\$"))
end

gif(MM_asyn_lam20, "Plots/GIFs/OT_MM_asyn_lam20.gif", fps = 10)


p2 = surface(M_lam20, tt, MM_corr_lam20', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p3 = contour(tt, M_lam20, MM_corr_lam20, fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))
p4 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :black, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
plot!(p4,tt, MM_corr_lam20[11,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))


# savefig(p2, "Plots/Instantaneous/OT_surf_lam20.svg")
# savefig(p3, "Plots/Instantaneous/OT_HM_lam20.png")
# savefig(p4, "Plots/Instantaneous/OT_opt_lam20.svg")


#---------------------------------------------------------------------------
# Asynchronous case, 1/λ = 50

lam = 50

Random.seed!(2)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(3)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

# Padding the matrix with NaNs
p1 = P_syn[Int.(floor.(t1)) .+ 1, 1]
p2 = P_syn[Int.(floor.(t2)) .+ 1, 2]

D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1 = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2 = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_asyn_lam50 = [p1 p2]
t_asyn_lam50 = [t1 t2]

# Quick EDA to determine optimal time-scale

Δt = collect(1:1:300)
res = zeros(length(Δt), 1)
for i in 1:length(Δt)
    NN = Int.(floor.(((nsim ./Δt[i]).-1.0) ./ 2))
    res[i] = NUFFTcorrDKFGG(P_asyn_lam50, t_asyn_lam50, N = NN)[1][1,2]
    # res[i] = NUFFTcorrFKFGG(P_asyn_lam50, t_asyn_lam50, N = NN)[1][1,2]
end

p1 = plot(Δt, res, legend = :bottomright, label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}", color = :blue, line=(1, [:solid]), dpi = 300)
vline!(p1, [220], label = L"\textrm{Choice of } \Delta t")
xlabel!(p1, L"\Delta t\textrm{[sec]}")
ylabel!(p1, L"\textrm{Integrated correlation } \rho_{\Delta t}^{12}")

# savefig(p1, "Plots/Instantaneous/OT_EDA_lam50.svg")


# Optimal time-scale for 1/λ = 10 is Δt = 220 [sec]
dt = 220
N = Int.(floor.(((nsim ./dt).-1.0) ./ 2))

M_lam50 = collect(1:1:Int(floor(N/2)))

MM_corr_lam50 = zeros(length(M_lam50), outlength)

for i in 1:length(M_lam50)
    MM_res = MM_inst(P_asyn_lam50, t_asyn_lam50, outlength, M = M_lam50[i], N = N)

    MM_corr_lam50[i,:] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
end



t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))

MM_asyn_lam50 = @animate for i in 1:length(M_lam50)
    now = M_lam50[i]
    p1 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
    plot!(p1,tt, MM_corr_lam50[i,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$M = $now\$"))
end

gif(MM_asyn_lam50, "Plots/GIFs/OT_MM_asyn_lam50.gif", fps = 5)


p2 = surface(M_lam50, tt, MM_corr_lam50', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p3 = contour(tt, M_lam50, MM_corr_lam50, fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

t = collect(1:1:nsim)/nsim
tt = collect(0:1/outlength:1-(1/outlength))
p4 = plot(t, deterministicCor.(t), label = L"\textrm{True } {\rho}^{12}(t)", color = :black, line=(1, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
plot!(p4,tt, MM_corr_lam50[10,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))


# savefig(p2, "Plots/Instantaneous/OT_surf_lam50.svg")
# savefig(p3, "Plots/Instantaneous/OT_HM_lam50.png")
# savefig(p4, "Plots/Instantaneous/OT_opt_lam50.svg")
