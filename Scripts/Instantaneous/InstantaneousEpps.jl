## Author: Patrick Chang
# Script file to investigate the impact of the Epps effect on the
# instantaneous dynamics of the Malliavin-Mancino estimator and the jump robust
# adaptation by Cuchiero and Teichmann.

# Note: Investigation will be done with a diffusion process
# with deterministic correlation; sampled to induce the Epps effect.
# Then Δt will range from 1-100 seconds with three choices of M.

# Note: The Malliavin-Mancino estimator has the ability to deal with asynchrony
# while the adaptation by Cuchiero and Teichmann requires the synchronisation
# of the data beforehand.

using LinearAlgebra, Plots, LaTeXStrings, StatsBase, Intervals, JLD, ProgressMeter, Distributions, CSV
#---------------------------------------------------------------------------

cd("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-Inst")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-JR")

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

#---------------------------------------------------------------------------
# Synchronous, ranging through M, N is Nyquist
#---------------------------------------------------------------------------
nsim = 28800
outlength = 1000

Hest_syn = Heston_CT(nsim, seed = 1, dt = nsim)
t = collect(1:1:nsim)

MM_corr_syn = zeros(100, 1000)
CT_corr_syn = zeros(100, 1000)

M = collect(1:1:100)

for i in 1:length(M)
    MM_res = MM_inst(Hest_syn[1], [t t], outlength, M = M[i])
    CT_res = MM_JR(Hest_syn[1], M[i], outlength)

    MM_corr_syn[i,:] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
    CT_corr_syn[i,:] = CT_res[3] ./ sqrt.(CT_res[1] .* CT_res[2])
end


tt = collect(0:1/outlength:1-(1/outlength))

SynM = @animate for i in 1:length(M)
    now = M[i]
    plot(collect(1:1:nsim)./nsim, Hest_syn[4] ./ sqrt.(Hest_syn[2] .* Hest_syn[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]))
    plot!(tt, MM_corr_syn[i,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    plot!(tt, CT_corr_syn[i,:], label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$M = $now\$"))
end

gif(SynM, "Plots/GIFs/SynchronousVaryingM.gif", fps = 10)

# MM
p1 = surface(M, tt, MM_corr_syn', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p2 = contour(tt, M, MM_corr_syn, fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

#CT
p3 = surface(M, tt, CT_corr_syn', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p4 = contour(tt, M, CT_corr_syn, fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

p5 = plot(collect(1:1:nsim)./nsim, Hest_syn[4] ./ sqrt.(Hest_syn[2] .* Hest_syn[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Instantaneous correlation } \rho^{12}(t)", dpi = 300)
plot!(p5,tt, MM_corr_syn[100,:], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
plot!(p5,tt, CT_corr_syn[100,:], label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]), size = (1200,400))

# savefig(p1, "Plots/Instantaneous/SynM_surf_MM.svg")
# savefig(p2, "Plots/Instantaneous/SynM_HM_MM.png")
# savefig(p3, "Plots/Instantaneous/SynM_surf_CT.svg")
# savefig(p4, "Plots/Instantaneous/SynM_HM_CT.png")
# savefig(p5, "Plots/Instantaneous/SynchronousCompM.svg")

#---------------------------------------------------------------------------
# Asynchronous, ranging through M, three choices of M
#---------------------------------------------------------------------------
nsim = 28800
outlength = 1000
lam = 30

Hest_syn = Heston_CT(nsim, seed = 1, dt = nsim)
t = collect(1:1:nsim)

# Sample times
Random.seed!(1)
t1 = [0; rexp(nsim, lam)]
t1 = cumsum(t1)
t1 = filter((x) -> x < nsim, t1)

Random.seed!(2)
t2 = [0; rexp(nsim, lam)]
t2 = cumsum(t2)
t2 = filter((x) -> x < nsim, t2)

# Padding the matrix with NaNs
p1 = Hest_syn[1][Int.(floor.(t1)) .+ 1, 1]
p2 = Hest_syn[1][Int.(floor.(t2)) .+ 1, 2]

t1_asyn = t1; t2_asyn = t2
D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
if length(t1) < length(t2)
    t1_asyn = [t1; repeat([NaN], D)]
    p1 = [p1; repeat([NaN], D)]
else
    t2_asyn = [t2; repeat([NaN], D)]
    p2 = [p2; repeat([NaN], D)]
end

P_Heston_asyn = [p1 p2]
t_Heston_asyn = [t1_asyn t2_asyn]

# Set up the experiment
M = [10; 20; 50]
dt = collect(1:1:100)

MM_corr_asyn = zeros(length(dt), 1000, length(M))
CT_corr_asyn = zeros(length(dt), 1000, length(M))

@showprogress "Computing..." for i in 1:length(dt)
    N = Int(floor(((nsim / dt[i])-1.0) / 2))
    Synchronised_Hest = synchronise(Hest_syn[1], t1, t2, nsim, dt[i])
    for j in 1:length(M)
        MM_res = MM_inst(P_Heston_asyn, t_Heston_asyn, outlength, M = M[j], N = N)
        CT_res = MM_JR(Synchronised_Hest, M[j], outlength)

        MM_corr_asyn[i,:,j] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
        CT_corr_asyn[i,:,j] = CT_res[3] ./ sqrt.(CT_res[1] .* CT_res[2])
    end
end

tt = collect(0:1/outlength:1-(1/outlength))

AsynM10 = @animate for i in 1:length(dt)
    nowΔt = dt[i]
    plot(collect(1:1:nsim)./nsim, Hest_syn[4] ./ sqrt.(Hest_syn[2] .* Hest_syn[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]))
    plot!(tt, MM_corr_asyn[i,:,1], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    plot!(tt, CT_corr_asyn[i,:,1], label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$\\Delta t = $nowΔt, M=10\$"))
end

AsynM20 = @animate for i in 1:length(dt)
    nowΔt = dt[i]
    plot(collect(1:1:nsim)./nsim, Hest_syn[4] ./ sqrt.(Hest_syn[2] .* Hest_syn[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]))
    plot!(tt, MM_corr_asyn[i,:,2], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    plot!(tt, CT_corr_asyn[i,:,2], label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$\\Delta t = $nowΔt, M=20\$"))
end

AsynM50 = @animate for i in 1:length(dt)
    nowΔt = dt[i]
    plot(collect(1:1:nsim)./nsim, Hest_syn[4] ./ sqrt.(Hest_syn[2] .* Hest_syn[3]), label = L"\textrm{True } {\rho}^{12}(t)", color = :lightblue, line=(0.5, [:solid]), ylims = (-0.65, 0.9))
    plot!(tt, MM_corr_asyn[i,:,3], label = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", color = :blue, line=(1, [:solid]))
    plot!(tt, CT_corr_asyn[i,:,3], label = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", color = :red, line=(1, [:dash]))
    xlabel!(L"\textrm{Time}")
    ylabel!(L"\textrm{Instantaneous correlation } \rho^{12}(t)")
    title!(latexstring("\$\\Delta t = $nowΔt, M=50\$"))
end

gif(AsynM10, "Plots/GIFs/AsynchronousVaryingN_M10.gif", fps = 10)
gif(AsynM20, "Plots/GIFs/AsynchronousVaryingN_M20.gif", fps = 10)
gif(AsynM50, "Plots/GIFs/AsynchronousVaryingN_M50.gif", fps = 10)

# MM plots
p1 = surface(dt, tt, MM_corr_asyn[:,:,1]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p2 = surface(dt, tt, MM_corr_asyn[:,:,2]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p3 = surface(dt, tt, MM_corr_asyn[:,:,3]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))

p4 = contour(tt, dt, MM_corr_asyn[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))
p5 = contour(tt, dt, MM_corr_asyn[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))
p6 = contour(tt, dt, MM_corr_asyn[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

# savefig(p1, "Plots/Instantaneous/Asyn_surf_MM_M10.svg")
# savefig(p2, "Plots/Instantaneous/Asyn_surf_MM_M20.svg")
# savefig(p3, "Plots/Instantaneous/Asyn_surf_MM_M50.svg")
# savefig(p4, "Plots/Instantaneous/Asyn_HM_MM_M10.png")
# savefig(p5, "Plots/Instantaneous/Asyn_HM_MM_M20.png")
# savefig(p6, "Plots/Instantaneous/Asyn_HM_MM_M50.png")


# CT plots
p1 = surface(dt, tt, CT_corr_asyn[:,:,1]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p2 = surface(dt, tt, CT_corr_asyn[:,:,2]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
p3 = surface(dt, tt, CT_corr_asyn[:,:,3]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{12}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))

p4 = contour(tt, dt, CT_corr_asyn[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))
p5 = contour(tt, dt, CT_corr_asyn[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))
p6 = contour(tt, dt, CT_corr_asyn[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, clim=(-1,1))

# savefig(p1, "Plots/Instantaneous/Asyn_surf_CT_M10.svg")
# savefig(p2, "Plots/Instantaneous/Asyn_surf_CT_M20.svg")
# savefig(p3, "Plots/Instantaneous/Asyn_surf_CT_M50.svg")
# savefig(p4, "Plots/Instantaneous/Asyn_HM_CT_M10.png")
# savefig(p5, "Plots/Instantaneous/Asyn_HM_CT_M20.png")
# savefig(p6, "Plots/Instantaneous/Asyn_HM_CT_M50.png")
