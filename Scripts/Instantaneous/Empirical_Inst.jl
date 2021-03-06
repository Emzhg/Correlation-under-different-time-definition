## Author: Patrick Chang
# Script file to investigate the intraday correlation using real TAQ
# data from the JSE for two stocks in the banking sector

# The data is cleaned by aggregating trades with the same trade time using
# a VWAP average.

using JLD; using LaTeXStrings; using Plots; using Statistics; using CSV;
using Optim; using Distributions; using ProgressMeter

#---------------------------------------------------------------------------

cd("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main")

include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-Inst.jl")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Instantaneous Estimators/MM-JR.jl")
include("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")

#---------------------------------------------------------------------------
# Read in the data
#---------------------------------------------------------------------------
days = ["2019-06-24"; "2019-06-25"; "2019-06-26"; "2019-06-27"; "2019-06-28"]

JSE_Data = Matrix{Float64}[]
for i in 1:length(days)
    p = CSV.read("/home/amo/dev/sunzulab/API/multifreq/PCRBTG-VT-main/Real Data/JSE_"*days[i]*".csv", DataFrame)
    p = convert(Matrix, p[:,2:12])
    push!(JSE_Data, p)
end

# Extract SBK and FSR prices
JSE_Price_SBKFSR = Matrix{Float64}[]
for i in 1:length(days)
    temp = JSE_Data[i][:,[7;11]]
    push!(JSE_Price_SBKFSR, temp)
end

# Create the appropriate time matrix for the two assets
# for the MM estimator
JSE_Time_SBKFSR = Matrix{Float64}[]
for i in 1:length(days)
    ind1 = findall(!isnan, JSE_Data[i][:,7])
    ind2 = findall(!isnan, JSE_Data[i][:,11])

    SBKtime = fill(NaN, size(JSE_Price_SBKFSR[i])[1], 1)
    FSRtime = fill(NaN, size(JSE_Price_SBKFSR[i])[1], 1)

    SBKtime[ind1] = JSE_Data[i][:,1][ind1]
    FSRtime[ind2] = JSE_Data[i][:,1][ind2]

    Time = [SBKtime FSRtime]

    push!(JSE_Time_SBKFSR, Time)
end

#---------------------------------------------------------------------------
# Compute the Epps curves for the 5 days to pick appropriate time-scales
#---------------------------------------------------------------------------

T = 28200   # 8 hour trading day minus 10 min closing auction
dt = collect(1:1:400)

MM_Int_cor = zeros(length(dt), length(days))

for i in 1:length(days)
    p = JSE_Price_SBKFSR[i]
    t = JSE_Time_SBKFSR[i]

    for j in 1:length(dt)
        N = Int(floor(((T /dt[j])-1) / 2))
        MM_Int_cor[j,i] = NUFFTcorrDKFGG(p, t, N = N)[1][1,2]
    end
end

days = ["2019-06-24" "2019-06-25" "2019-06-26" "2019-06-27" "2019-06-28"]
weekdays = ["Monday " "Tuesday " "Wednesday " "Thursday " "Friday "]
day_labels = weekdays .* days

p1 = plot(dt, MM_Int_cor, label = day_labels, legend = :bottomright, dpi = 300)
xlabel!(p1, L"\Delta t\textrm{[sec]}")
ylabel!(p1, L"\textrm{Integrated correlation } \hat{\rho}_{\Delta t}^{SBK/FSR}")

# savefig(p1, "Plots/Instantaneous/Emp_Int_EppsCurves.svg")

#---------------------------------------------------------------------------
# Additional supporting functions for previous tick interpolation
#---------------------------------------------------------------------------
## Use Previous tick interpolaton to create the price paths for the CT estimator

# Function to clean the data such that both assets start when the slower
# asset makes its first trade. Using Previous tick interpolation to interpolate
# for the asset which traded first, if the faster asset does not trade at the
# same time as the first trade of the slower asset.
function fixdata(data)
    m = size(data)[1]
    Fixed_Data = Matrix{Float64}[]
    # Fixed_Data = Vector{Vector{Float64}}()
    for i in 1:m
        n = size(data[i])[1]
        indA1 = findall(!isnan, data[i][:,2])
        indA2 = findall(!isnan, data[i][:,3])
        indA1_min = minimum(indA1)
        indA2_min = minimum(indA2)

        if indA1_min==indA2_min
            temp_data = data[i][indA1_min:end,:]
        end

        if indA1 > indA2
            if isnan(data[i][indA1_min,3])
                data[i][indA1_min,3] = data[i][maximum(filter(x-> x < indA1_min, indA2)),3]
            end
            temp_data = data[i][indA1_min:end,:]
        else
            if isnan(data[i][indA2_min,2])
                data[i][indA2_min,2] = data[i][maximum(filter(x-> x < indA2_min, indA1)),2]
            end
            temp_data = data[i][indA2_min:end,:]
        end
        time = temp_data[:,1]
        mintime = minimum(time)
        time = time .- mintime
        temp_data[:,1] = time

        push!(Fixed_Data, temp_data)
    end
    return Fixed_Data
end

# Function to synchronise the data at a particular time-scale ??t using the
# previous tick interpolation for the CT estimator
function synchronise(data, dt, T = 28200)
    m = size(data)[1]
    Synchronised_Data = Matrix{Float64}[]

    t = collect(0:dt:T)
    n = length(t)

    for k in 1:m
        temp = data[k]
        t1 = temp[findall(!isnan, temp[:,2]),1]
        t2 = temp[findall(!isnan, temp[:,3]),1]
        p1 = zeros(n,1)
        p2 = zeros(n,1)
        for j in 1:n
            inds = findall(x -> x .<= t[j], temp[:,1])
            p1[j] = filter(!isnan, temp[inds,2])[end]
            p2[j] = filter(!isnan, temp[inds,3])[end]
        end
        p = [p1 p2]
        push!(Synchronised_Data, p)
    end
    return Synchronised_Data
end

#---------------------------------------------------------------------------
# Obtain the synchronised SBK and FSR data
Data_SBKFSR_PTI = Matrix{Float64}[]
for i in 1:length(days)
    p = JSE_Data[i][:,[1;7;11]] # Index for SBK and FSR is 7 and 11 resp.
    push!(Data_SBKFSR_PTI, p)
end
# Fix the opening trades for each day, so they line up
Data_fixed_SBKFSR_PTI = fixdata(Data_SBKFSR_PTI)

#---------------------------------------------------------------------------
# Investigate the instantaneous Epps effect
#---------------------------------------------------------------------------
outlength = 499
T = 28200   # 8 hour trading day minus 10 min closing auction
dt = collect(1:1:400)

# MM Results
MM_Epps_inst_M10 = zeros(length(dt), outlength, length(days))
MM_Epps_inst_M20 = zeros(length(dt), outlength, length(days))

for j in 1:length(days)
    p = JSE_Price_SBKFSR[j]
    t = JSE_Time_SBKFSR[j]
    for i in 1:length(dt)
        N = Int(floor(((T / dt[i])-1) / 2))
        MM_res_M10 = MM_inst(p, t, outlength, M = 10, N = N)
        MM_res_M20 = MM_inst(p, t, outlength, M = 20, N = N)

        MM_Epps_inst_M10[i,:,j] = MM_res_M10[3] ./ sqrt.(MM_res_M10[1] .* MM_res_M10[2])
        MM_Epps_inst_M20[i,:,j] = MM_res_M20[3] ./ sqrt.(MM_res_M20[1] .* MM_res_M20[2])
    end
end

# CT Results: Note that end points for Fourier method is unstable
# so we remove the point 0
CT_Epps_inst_M10 = zeros(length(dt), outlength-1, length(days))
CT_Epps_inst_M20 = zeros(length(dt), outlength-1, length(days))

# takes around 5 min to compute, first few ??t are slow but it gets faster
# for larger ??t
@showprogress "Computing..." for i in 1:length(dt)
    p = synchronise(Data_fixed_SBKFSR_PTI, dt[i])
    for j in 1:length(days)
        pp = p[j]
        CT_res_M10 = MM_JR(pp, 10, outlength)
        CT_res_M20 = MM_JR(pp, 20, outlength)

        CT_Epps_inst_M10[i,:,j] = CT_res_M10[3][2:end] ./ sqrt.(CT_res_M10[1][2:end] .* CT_res_M10[2][2:end])
        CT_Epps_inst_M20[i,:,j] = CT_res_M20[3][2:end] ./ sqrt.(CT_res_M20[1][2:end] .* CT_res_M20[2][2:end])
    end
end

## Plot the results
# MM results
tt_MM = collect(0:1/outlength:1-(1/outlength))

p1 = surface(dt, tt_MM, MM_Epps_inst_M10[:,:,1]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p2 = surface(dt, tt_MM, MM_Epps_inst_M10[:,:,2]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p3 = surface(dt, tt_MM, MM_Epps_inst_M10[:,:,3]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p4 = surface(dt, tt_MM, MM_Epps_inst_M10[:,:,4]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p5 = surface(dt, tt_MM, MM_Epps_inst_M10[:,:,5]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))

p6 = contour(tt_MM, dt, MM_Epps_inst_M10[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p7 = contour(tt_MM, dt, MM_Epps_inst_M10[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p8 = contour(tt_MM, dt, MM_Epps_inst_M10[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p9 = contour(tt_MM, dt, MM_Epps_inst_M10[:,:,4], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p10 = contour(tt_MM, dt, MM_Epps_inst_M10[:,:,5], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))

# savefig(p1, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M10_M.svg")
# savefig(p2, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M10_Tues.svg")
# savefig(p3, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M10_W.svg")
# savefig(p4, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M10_Thurs.svg")
# savefig(p5, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M10_F.svg")
#
# savefig(p6, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M10_M.png")
# savefig(p7, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M10_Tues.png")
# savefig(p8, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M10_W.png")
# savefig(p9, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M10_Thurs.png")
# savefig(p10, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M10_F.png")

p1 = surface(dt, tt_MM, MM_Epps_inst_M20[:,:,1]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p2 = surface(dt, tt_MM, MM_Epps_inst_M20[:,:,2]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p3 = surface(dt, tt_MM, MM_Epps_inst_M20[:,:,3]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p4 = surface(dt, tt_MM, MM_Epps_inst_M20[:,:,4]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p5 = surface(dt, tt_MM, MM_Epps_inst_M20[:,:,5]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))

p6 = contour(tt_MM, dt, MM_Epps_inst_M20[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p7 = contour(tt_MM, dt, MM_Epps_inst_M20[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p8 = contour(tt_MM, dt, MM_Epps_inst_M20[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p9 = contour(tt_MM, dt, MM_Epps_inst_M20[:,:,4], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p10 = contour(tt_MM, dt, MM_Epps_inst_M20[:,:,5], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))

# savefig(p1, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M20_M.svg")
# savefig(p2, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M20_Tues.svg")
# savefig(p3, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M20_W.svg")
# savefig(p4, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M20_Thurs.svg")
# savefig(p5, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_M20_F.svg")
#
# savefig(p6, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M20_M.png")
# savefig(p7, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M20_Tues.png")
# savefig(p8, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M20_W.png")
# savefig(p9, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M20_Thurs.png")
# savefig(p10, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_M20_F.png")

# CT results
tt_CT = collect(1/outlength:1/outlength:1-(1/outlength))

p1 = surface(dt, tt_CT, CT_Epps_inst_M10[:,:,1]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p2 = surface(dt, tt_CT, CT_Epps_inst_M10[:,:,2]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p3 = surface(dt, tt_CT, CT_Epps_inst_M10[:,:,3]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p4 = surface(dt, tt_CT, CT_Epps_inst_M10[:,:,4]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p5 = surface(dt, tt_CT, CT_Epps_inst_M10[:,:,5]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))

p6 = contour(tt_CT, dt, CT_Epps_inst_M10[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p7 = contour(tt_CT, dt, CT_Epps_inst_M10[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p8 = contour(tt_CT, dt, CT_Epps_inst_M10[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p9 = contour(tt_CT, dt, CT_Epps_inst_M10[:,:,4], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p10 = contour(tt_CT, dt, CT_Epps_inst_M10[:,:,5], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))

# savefig(p1, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M10_M.svg")
# savefig(p2, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M10_Tues.svg")
# savefig(p3, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M10_W.svg")
# savefig(p4, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M10_Thurs.svg")
# savefig(p5, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M10_F.svg")
#
# savefig(p6, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M10_M.png")
# savefig(p7, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M10_Tues.png")
# savefig(p8, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M10_W.png")
# savefig(p9, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M10_Thurs.png")
# savefig(p10, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M10_F.png")

p1 = surface(dt, tt_CT, CT_Epps_inst_M20[:,:,1]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p2 = surface(dt, tt_CT, CT_Epps_inst_M20[:,:,2]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p3 = surface(dt, tt_CT, CT_Epps_inst_M20[:,:,3]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p4 = surface(dt, tt_CT, CT_Epps_inst_M20[:,:,4]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p5 = surface(dt, tt_CT, CT_Epps_inst_M20[:,:,5]', xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, xticks = ([100:100:300;], ["100", "200", "300"]), yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))

p6 = contour(tt_CT, dt, CT_Epps_inst_M20[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p7 = contour(tt_CT, dt, CT_Epps_inst_M20[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p8 = contour(tt_CT, dt, CT_Epps_inst_M20[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p9 = contour(tt_CT, dt, CT_Epps_inst_M20[:,:,4], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p10 = contour(tt_CT, dt, CT_Epps_inst_M20[:,:,5], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\Delta t \textrm{[sec]}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))

# savefig(p1, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M20_M.svg")
# savefig(p2, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M20_Tues.svg")
# savefig(p3, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M20_W.svg")
# savefig(p4, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M20_Thurs.svg")
# savefig(p5, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_M20_F.svg")
#
# savefig(p6, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M20_M.png")
# savefig(p7, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M20_Tues.png")
# savefig(p8, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M20_W.png")
# savefig(p9, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M20_Thurs.png")
# savefig(p10, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_M20_F.png")

#---------------------------------------------------------------------------
# Compute the MM and CT estimator for the appropriate time-scales
#---------------------------------------------------------------------------
??t = 300    # obtained from Epps curves
Synchronised_SBKFSR_PTI = synchronise(Data_fixed_SBKFSR_PTI, ??t)
#---------------------------------------------------------------------------
## Compute the results for MM and CT estimator for various values of M
outlength = 499
N = Int(floor(((T /??t)-1) / 2))
M = collect(1:1:Int(floor(N/2)))

## MM results
MM_corr_inst = zeros(length(M), outlength, length(days))

for j in 1:length(days)
    p = JSE_Price_SBKFSR[j]
    t = JSE_Time_SBKFSR[j]
    for i in 1:length(M)
        MM_res = MM_inst(p, t, outlength, M = M[i], N = N)

        MM_corr_inst[i,:,j] = MM_res[3] ./ sqrt.(MM_res[1] .* MM_res[2])
    end
end

## CT results
# Note that end points for Fourier method is unstable
# so we remove the point 0
outlength = 499
CT_corr_inst = zeros(length(M), outlength-1, length(days))

for j in 1:length(days)
    p = Synchronised_SBKFSR_PTI[j]
    for i in 1:length(M)
        CT_res = MM_JR(p, M[i], outlength)

        CT_corr_inst[i,:,j] = CT_res[3][2:end] ./ sqrt.(CT_res[1][2:end] .* CT_res[2][2:end])
    end
end

## Plot the results
# MM results
tt_MM = collect(0:1/outlength:1-(1/outlength))
p1 = surface(M, tt_MM, MM_corr_inst[:,:,1]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p2 = surface(M, tt_MM, MM_corr_inst[:,:,2]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p3 = surface(M, tt_MM, MM_corr_inst[:,:,3]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p4 = surface(M, tt_MM, MM_corr_inst[:,:,4]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p5 = surface(M, tt_MM, MM_corr_inst[:,:,5]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))

p6 = contour(tt_MM, M, MM_corr_inst[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p7 = contour(tt_MM, M, MM_corr_inst[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p8 = contour(tt_MM, M, MM_corr_inst[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p9 = contour(tt_MM, M, MM_corr_inst[:,:,4], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p10 = contour(tt_MM, M, MM_corr_inst[:,:,5], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))

p11 = plot(tt_MM, MM_corr_inst[10,:,:], label=day_labels, legend = :bottomleft, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", dpi = 300)

# savefig(p1, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_dt300_M.svg")
# savefig(p2, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_dt300_Tues.svg")
# savefig(p3, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_dt300_W.svg")
# savefig(p4, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_dt300_Thurs.svg")
# savefig(p5, "Plots/Instantaneous/Emp_MM_Inst_Surf_Epps_dt300_F.svg")
#
# savefig(p6, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_dt300_M.png")
# savefig(p7, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_dt300_Tues.png")
# savefig(p8, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_dt300_W.png")
# savefig(p9, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_dt300_Thurs.png")
# savefig(p10, "Plots/Instantaneous/Emp_MM_Inst_HM_Epps_dt300_F.png")


# CT results
tt_CT = collect(1/outlength:1/outlength:1-(1/outlength))

p1 = surface(M, tt_CT, CT_corr_inst[:,:,1]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p2 = surface(M, tt_CT, CT_corr_inst[:,:,2]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p3 = surface(M, tt_CT, CT_corr_inst[:,:,3]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p4 = surface(M, tt_CT, CT_corr_inst[:,:,4]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))
p5 = surface(M, tt_CT, CT_corr_inst[:,:,5]', xlabel = L"\textrm{M}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", fc=ColorGradient([:red,:yellow,:blue]), camera = (65, 55), dpi = 300, yticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1), zlims = (-1, 1))

p6 = contour(tt_CT, M, CT_corr_inst[:,:,1], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p7 = contour(tt_CT, M, CT_corr_inst[:,:,2], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p8 = contour(tt_CT, M, CT_corr_inst[:,:,3], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p9 = contour(tt_CT, M, CT_corr_inst[:,:,4], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))
p10 = contour(tt_CT, M, CT_corr_inst[:,:,5], fill = true, fc = ColorGradient([:red,:yellow,:blue]), ylabel = L"\textrm{M}", xlabel = L"\textrm{Time}", dpi = 150, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), clim=(-1,1))

p11 = plot(tt_CT, CT_corr_inst[10,:,:], label=day_labels, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Estimated CT } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)")

# savefig(p1, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_dt300_M.svg")
# savefig(p2, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_dt300_Tues.svg")
# savefig(p3, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_dt300_W.svg")
# savefig(p4, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_dt300_Thurs.svg")
# savefig(p5, "Plots/Instantaneous/Emp_CT_Inst_Surf_Epps_dt300_F.svg")
#
# savefig(p6, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_dt300_M.png")
# savefig(p7, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_dt300_Tues.png")
# savefig(p8, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_dt300_W.png")
# savefig(p9, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_dt300_Thurs.png")
# savefig(p10, "Plots/Instantaneous/Emp_CT_Inst_HM_Epps_dt300_F.png")

p1 = plot(tt_MM, MM_corr_inst[10,:,2], label=L"\textrm{Estimated MM Tuesday 2019-06-25}", legend = :bottom, xticks = ([0:0.5:1;], ["09:00", "12:55", "16:50"]), xlabel = L"\textrm{Time}", ylabel = L"\textrm{Estimated MM } \hat{\rho}^{SBK/FSR}_{\Delta t}(t)", dpi = 300, color = :blue, line=(1, [:solid]))
plot!(p1, tt_MM, MM_corr_inst[10,:,3], label=L"\textrm{Estimated MM Wednesday 2019-06-26}", color = :red, line=(1, [:solid]))
plot!(p1, tt_CT, CT_corr_inst[10,:,2], label=L"\textrm{Estimated CT Tuesday 2019-06-25}", color = :blue, line=(1, [:dash]))
plot!(p1, tt_CT, CT_corr_inst[10,:,3], label=L"\textrm{Estimated CT Wednesday 2019-06-26}", color = :red, line=(1, [:dash]))

# savefig(p1, "Plots/Instantaneous/Emp_Inst_Epps_Comp_dt300_M10.svg")
