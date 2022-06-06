## Modified to fit purposes of this project
using LinearAlgebra, Plots, LaTeXStrings, StatsBase, Intervals, JLD, ProgressMeter, Distributions, CSV, DataFrames, Dates, Pipe, ColorSchemes
#---------------------------------------------------------------------------

cd("C:/Users/User/Desktop/Dauphine/Semestre 4/Projet - GQ/PCRBTG-VT-main/PCRBTG-VT-main")

include("../Functions/Instantaneous Estimator/MM-Inst.jl") ## JSP: added the Instantaneous function

# Read in the data
prices = CSV.read("Real Data/JSE_prices_2019-06-24_2019-06-28.csv",DataFrame)
times = CSV.read("Real Data/JSE_times_2019-06-24_2019-06-28.csv",DataFrame)
volume = CSV.read("Real Data/JSE_volume_2019-06-24_2019-06-28.csv",DataFrame)

# Remove extra column
prices = prices[:,2:end]; times = times[:,2:end]; volume = volume[:,2:end]
# Pull out banking stocks
prices = prices[:,[6:8; 10]]; times = times[:,[6:8; 10]]; volume = volume[:,[6:8; 10]]
# Get the names of tickers
tickers = names(prices)

#---------------------------------------------------------------------------
## Supporting functions
function MakeMMHYData(data)
    # Pull out the data
    A1 = data[1]; A2 = data[2]
    # Extract the prices and time
    p1 = A1[:,1]; p2 = A2[:,1]
    t1 = A1[:,2]; t2 = A2[:,2]
    # Re-scale the times
    T0 = min(t1[1], t2[1])
    t1 = t1 .- T0; t2 = t2 .- T0
    # Combine and sort the times
    tt = unique([t1; t2])
    tt = sort(tt)
    # Initialize price and time matrix
    n = length(tt)
    P = fill(NaN, n, 2); t = fill(NaN, n, 2)
    # Loop through the unique times and order the events
    for i in 1:length(tt)
        # Logic checks
        if tt[i] ∈ t1
            # Get index
            ind = findall(x -> x == tt[i], t1)
            # Current time was in t1
            P[i,1] = p1[ind][1]
            t[i,1] = i
        end
        if tt[i] ∈ t2
            # Get index
            ind = findall(x -> x == tt[i], t2)
            # Current time was in t1
            P[i,2] = p2[ind][1]
            t[i,2] = i
        end
    end
    return P, t
end

function DataSplit(A1::String, A2::String, P::DataFrame, t::DataFrame, V::DataFrame)
    # Filter out the pair of interest
    p1 = filter(!isnan, P[:,A1]); p2 = filter(!isnan, P[:,A2])
    t1 = filter(!isnan, t[:,A1]); t2 = filter(!isnan, t[:,A2])
    V1 = filter(!isnan, V[:,A1]); V2 = filter(!isnan, V[:,A2])
    # Convert the times to dates for index extraction later on
    A1dates = Date.(unix2datetime.(t1))
    A2dates = Date.(unix2datetime.(t2))
    dates_unique = unique(A1dates)
    # Initialize storage
    data = Dict()
    # Loop through each day
    for i in 1:length(dates_unique)
        # Extract the data
        date_indsA1 = findall(x -> x == dates_unique[i], A1dates)
        date_indsA2 = findall(x -> x == dates_unique[i], A2dates)
        # Data for the day
        day_data = []
        push!(day_data, [p1[date_indsA1] t1[date_indsA1] V1[date_indsA1]])
        push!(day_data, [p2[date_indsA2] t2[date_indsA2] V2[date_indsA2]])
        # Add to dictionary
        push!(data, i => day_data)
    end
    return data
end

function getCTcorrs(tickers, outlength; P=prices, t=times, V=volume)
     # Compute the number of pairwise comparisons
     npairs = Int(factorial(length(tickers)) / (factorial(2) * factorial(length(tickers)-2)))
     # Time scale of investigation
     dt = collect(1:1:300)
     # Initialize the estimates   
     MMinst = zeros(length(dt),npairs,outlength) ## JSP: initialize the inst. estimate
    # Set up ind for storage
    ind = 1
     # Loop through the pairs
     @showprogress "Computing..." for i in 1:(length(tickers)-1)
        for j in (i+1):length(tickers)
            # Split the data into separate days
            data = DataSplit(tickers[i], tickers[j], P, t, V)
            # Initialize temporary storage matricies for the estimates
            MMinsttemp = zeros(length(dt), 5,outlength) ## JSP: initialize temporary storage matrix for inst. estimate
            # Loop through the different days
            for k in 1:length(data)
                # Extract data for the day
                day_data = data[k]
                # Create the MM and HY dataset
                MMHYData = MakeMMHYData(day_data)
                # Loop through the different time scales
                for l in 1:length(dt)
                    # Get N for MM
                    N = Int(floor((28200/dt[l] - 1)/2))
                    # Compute correlations
                    MMinsttemp_aux=MM_inst(MMHYData[1],MMHYData[2],outlength,M= 10,N=N)##JSP: estimate co-vol.
                    MMinsttemp[l,k,:] = MMinsttemp_aux[3]./ MMinsttemp_aux[1].* MMinsttemp_aux[2]##JSP: Estimate corr.
                end
            end
            MMinst[:,ind,:]= mean(MMinsttemp,dims=2)
            ind +=1
        end
     end
     return MMinst
end

###--------------------------------------------------------------------------
##                  EXECUTION DU PROG     ##
###--------------------------------------------------------------------------
outlength = 500
dt = collect(1:1:300)
tt = collect(0:1/outlength:1-(1/outlength))
res = getCTcorrs(tickers, outlength)
#---------------------------------------------------------------------------
pairnames = Matrix{Union{Nothing, String}}(nothing, 1, 6)
let inds = 1
    for i in 1:(length(tickers)-1)
        for j in (i+1):length(tickers)
            # push!(pairnames, "$(tickers[i])"*"/"*"$(tickers[j])")
            pairnames[inds] = "$(tickers[i])"*"/"*"$(tickers[j])"
            inds += 1
        end
    end
end

p1 = surface(dt, tt, res[:,1,:]'.*1000, xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
title!(p1,pairnames[1])
savefig(p1, "Plots/ET_INSTANT_1.svg")

p2 = surface(dt, tt, res[:,2,:]'.*1000, xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
title!(p2,pairnames[2])
savefig(p2, "Plots/ET_INSTANT_2.svg")

p3 = surface(dt, tt, res[:,3,:]'.*1000, xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
title!(p3,pairnames[3])
savefig(p3, "Plots/ET_INSTANT_3.svg")

p4 = surface(dt, tt, res[:,4,:]'.*1000, xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
title!(p4,pairnames[4])
savefig(p4, "Plots/ET_INSTANT_4.svg")

p5 = surface(dt, tt, res[:,5,:]'.*1000, xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
title!(p5,pairnames[5])
savefig(p5, "Plots/ET_INSTANT_5.svg")

p6 = surface(dt, tt, res[:,6,:]'.*1000, xlabel = L"\Delta t \textrm{[sec]}", ylabel = L"\textrm{Time}", zlabel = L"\textrm{Estimated MM } \hat{\rho}^{12}_{\Delta t}(t)", camera = (65, 55), dpi = 300, clim=(-1,1), zlims = (-1, 1))
title!(p6,pairnames[6])
savefig(p6, "Plots/ET_INSTANT_6.svg")