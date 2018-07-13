using MendelBase
using SnpArrays
using IHT
using DataFrames
using Distances
using StatsBase
#using Plotly
using PyPlot
#using StatPlots
using Plots
using StatPlots
using BenchmarkTools
using StandardizedMatrices
using Distributions
using JLD
using NPZ
#export MendelIHT

include("MendelIHT_utilities.jl")

"""
This is the wrapper function for the Iterative Hard Thresholding analysis option in Open Mendel.
"""
function MendelIHT(control_file = ""; args...)
    const MENDEL_IHT_VERSION :: VersionNumber = v"0.1.0"
    #
    # Print the logo. Store the initial directory.
    #
    print(" \n \n")
    println("     Welcome to OpenMendel's")
    println("      IHT analysis option")
    println("        version ", MENDEL_IHT_VERSION)
    print(" \n \n")
    println("Reading the data.\n")
    initial_directory = pwd()
    #
    # The user specifies the analysis to perform via a set of keywords.
    # Start the keywords at their default values.
    #
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    #
    # Define some keywords unique to this analysis option.
    #
    keyword["data_type"] = ""
    keyword["predictors"] = ""
    keyword["manhattan_plot_file"] = ""
    #
    # Process the run-time user-specified keywords that will control the analysis.
    # This will also initialize the random number generator.
    #
    process_keywords!(keyword, control_file, args)
    #
    # Check that the correct analysis option was specified.
    #
    lc_analysis_option = lowercase(keyword["analysis_option"])
    if (lc_analysis_option != "" &&
      lc_analysis_option != "iht")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
    end
    keyword["analysis_option"] = "Iterative Hard Thresholding"
    #
    # Read the genetic data from the external files named in the keywords.
    #
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
    #
    # Execute the specified analysis.
    #
    println(" \nAnalyzing the data.\n")
##
    phenotype = convert(Array{Float64,1}, pedigree_frame[:Trait])
    k = keyword["predictors"]
    result = L0_reg(snpdata, phenotype, k)
    return result
##
    # execution_error = iht_gwas(person, snpdata, pedigree_frame, keyword)
    # if execution_error
    #   println(" \n \nERROR: Mendel terminated prematurely!\n")
    # else
    #   println(" \n \nMendel's analysis is finished.\n")
    # end
    # #
    # # Finish up by closing, and thus flushing, any output files.
    # # Return to the initial directory.
    # #
    # close(keyword["output_unit"])
    # cd(initial_directory)
    # return nothing
end #function IterHardThreshold

"""
This function performs IHT on GWAS data.
"""
function L0_reg(
    x        :: SnpData,
    y        :: Vector{Float64},
    k        :: Int;
    v        :: IHTVariable = IHTVariables(x, y, k),
    mask_n   :: Vector{Int} = ones(Int, size(y)),
    tol      :: Float64 = 1e-4,
    max_iter :: Int     = 100,
    max_step :: Int     = 50,
)

    # start timer
    tic()

    # first handle errors
    k        >= 0            || throw(ArgumentError("Value of k must be nonnegative!\n"))
    max_iter >= 0            || throw(ArgumentError("Value of max_iter must be nonnegative!\n"))
    max_step >= 0            || throw(ArgumentError("Value of max_step must be nonnegative!\n"))
    tol      >  eps(Float64) || throw(ArgumentError("Value of global tol must exceed machine precision!\n"))

    # initialize return values
    mm_iter   = 0                 # number of iterations of L0_reg
    mm_time   = 0.0               # compute time *within* L0_reg
    next_loss = oftype(tol,Inf)   # loss function value

    # initialize floats
    current_obj = oftype(tol,Inf) # tracks previous objective function value
    the_norm    = 0.0             # norm(b - b0)
    scaled_norm = 0.0             # the_norm / (norm(b0) + 1)
    μ           = 0.0             # Landweber step size, 0 < tau < 2/rho_max^2

    # initialize integers
    i       = 0                   # used for iterations in loops
    mu_step = 0                   # counts number of backtracking steps for mu

    # initialize booleans
    converged = false             # scaled_norm < tol?

    #convert bitarrays to Float64 genotype matrix, standardize each SNP, and add intercept
    snpmatrix = convert(Array{Float64,2}, x.snpmatrix)
    const RED     = "\x1b[31m"
    const GREEN   = "\x1b[32m"
    const YELLOW  = "\x1b[33m"
    const BLUE    = "\x1b[34m"
    const PURPLE  = "\x1b[35m"
    const BOLD    = "\x1b[1m"
    const DEFAULT = "\x1b[0m"
    println(BLUE*"Gordon 124: new code for Prior Weights work"*DEFAULT)
    println("Note: numbers after my name should be close to line numbers in original code.")

    xxx = SnpArray("gwas 1 data")
    xxx = convert(Array{Float64,2}, xxx)
    println("Gordon 124: contents of snpmatrix from SnpArray()")
    println(xxx[1,1:10])

    println("Gordon 127: contents of snpmatrix before zscore()")
    println(snpmatrix[1,1:10])
    println(snpmatrix[2,1:10])
    println(snpmatrix[3,1:10])
    println(snpmatrix[4,1:10])
    println(snpmatrix[5,1:10])
    SMz = StandardizedMatrix(snpmatrix)         # WORKS !!!
    println("\nGordon 130: contents of snpmatrix_ztemp using NEW StandardizedMatrix()")
    println(SMz[1,1:10])
    println(SMz[2,1:10])
    println(SMz[3,1:10])
    println(SMz[4,1:10])
    println(SMz[5,1:10])
    snpmatrix_ztemp = StatsBase.zscore(snpmatrix, 1)  # standardize
    println("\nGordon 130: contents of snpmatrix_ztemp after zscore()")
    println(snpmatrix_ztemp[1,1:10])
    println(snpmatrix_ztemp[2,1:10])
    println(snpmatrix_ztemp[3,1:10])
    println(snpmatrix_ztemp[4,1:10])
    println(snpmatrix_ztemp[5,1:10])
    snpmatrix_ztemp2 = StatsBase.zscore(snpmatrix[:,1:10], 1)  # standardize
    println(size(snpmatrix[1:3,:]))
    println(size(snpmatrix_ztemp2))
    println("\nGordon 130: Just verifying that zscore() is being applied by column (per SNP).")
    println(snpmatrix_ztemp2[1,1:10])
    println(snpmatrix_ztemp2[2,1:10])
    println(snpmatrix_ztemp2[3,1:10])
    println(snpmatrix_ztemp2[4,1:10])
    println(snpmatrix_ztemp2[5,1:10])
    println("\n")

    # Create my own snpmatrix for Prior Weights work
    my_snpmatrix = zeros(snpmatrix)
    copy!(my_snpmatrix, snpmatrix)
    my_snpmatrix = deepcopy(snpmatrix)  # ask Ben about copying matricies
    println(my_snpmatrix[1,1:20])       # all 0,1,2 as Ben said

    my_snpsum = sum(my_snpmatrix, 1)
    print("Gordon - size(my_snpsum) = ")
    println(size(my_snpsum))
    print("Gordon - my_snpsum = ")
    println(my_snpsum[1,1:10] / 2200)       #
    print("Gordon - max(my_snpsum) = ")
    println(maximum(my_snpsum))
    println("Gordon - describe(my_snpsum[1,:] / 2200) - need for making my_snpweights_p")
    describe(my_snpsum[1,:] / 2200)

    # GORDON - CALCUATE WEIGHTS
    my_snpweights_p = my_snpsum / 2200      # p_hat
    my_snpweights_p /= 2    # a bug fix for p is MAF s/b max .5 per Ben
                            # also is was getting some NaN in my_diagm_snpweights
                            # because of some inf in my_snpweights
                            # that I imagine came from some √(0.0*(1.0-1.0)) errors
    my_snpweights_q = 1 - my_snpweights_p   # q_hat
    my_snpweights = my_snpweights_p + my_snpweights_q   # just verifying p + q == 1 OK!
    my_snpweights = my_snpweights_p .* my_snpweights_q   # just verifying P .* q max is 0.25 OK!
    # next line is 1/pq
    my_snpweights = my_snpweights .\ 1      # not what we want, but this works! to get invervse of each element
    my_snpweights = sqrt(my_snpweights_p .* my_snpweights_q)   # just verifying sqrtm(p .* q) == 0.5 OK!
    # next line is 2 * sqrt(p(1-p)) from Hua Zhou 2011 paper pg. 110
    my_snpweights = 2 * sqrt(my_snpweights_p .* (1 - my_snpweights_p))   # just verifying 2 * sqrtm(p .* q) == 1.0 OK!
    println("\ndescribe(my_snpweights_p)")
    describe(my_snpweights_p[1,:])
    # THESE COMMENTED DESCRIPTIONS ARE FROM BEFORE I CUT P IN HALF !!!
    # Mean:           0.823944
    # Minimum:        0.000000
    # 1st Quartile:   0.731522
    # Median:         0.889327
    # 3rd Quartile:   0.973884
    # Maximum:        1.000000
    println("\ndescribe(my_snpweights_q)")
    describe(my_snpweights_q[1,:])
    #println("describe(1 - my_snpweights_p)")   # SAME AS Q
    #describe(1 .- my_snpweights_p[1,:])
    println("\ndescribe(my_snpweights)")
    describe(my_snpweights[1,:])
    my_snpweights = my_snpweights .\ 1      # this works! to get reciprocal of each element
    println("\ndescribe(my_snpweights .\\ 1)")
    describe(my_snpweights[1,:])
    # THESE COMMENTED DESCRIPTIONS ARE FROM BEFORE I CUT P IN HALF !!!
    # Mean:           Inf
    # Minimum:        1.000000
    # 1st Quartile:   1.026816
    # Median:         1.124446
    # 3rd Quartile:   1.367013
    # Maximum:        Inf
    println()
    # =================================================================
    data = [
      [
        "y" => my_snpsum[1,:],
        "boxpoints" => "all",
        "jitter" => 0.3,
        "pointpos" => -1.8,
        "type" => "box"
      ]
    ]
#    response = Plotly.plot(data, ["filename" => "box-plot-jitter", "fileopt" => "overwrite"])
#    plot_url = response["url"]

#    boxplot([1,1,2,2,3,4,5,10], # Each column/cell is one box
#    	notch=true, # Notched center
#    	whis=0.75, # Whisker length as a percent of inner quartile range
#    	widths=0.25, # Width of boxes
#    	vert=false, # Horizontal boxes
#    	sym="rs") # Symbol color and shape (rs = red square)
    StatPlots.boxplot(rand(10000,1))

    # h = PyPlot.plt.hist(x,nbins)
    #%matplotlib inline
    #default(show=true)
    #plot(rand(100,4), show=true)
    StatPlots.boxplot(rand(100,3)) # Plots.boxplot -> The applicable method may
                                #be too new: running in world age 22852, while current world is 22853
    #Plots.boxplot(rand(100,3)) # Plots.boxplot -> The applicable method may
    Plots.savefig("plot222.png")
    Plots.savefig("plot333.png")
    #display(plt)
    show()

    # =================================================================
    println(BLUE*"HERE ARE THE SUMMARY STATS FOR THE WEIGHTS"*DEFAULT)
    print("Gordon - size(my_snpweights) = ")    # (1, 10000)
    println(size(my_snpweights))
    print("Gordon - my_snpweights = ")
    println(my_snpweights[1,1:10])
    print("Gordon - max(my_snpweights) = ")
    println(maximum(my_snpweights))
    describe(my_snpweights[1,:])

    println(BLUE*"THE SUMMARY STATS FOR THE WEIGHTS ARE PLOTTED IN plot_my_snpweights.png"*DEFAULT)
    StatPlots.boxplot(my_snpweights[1,:])
    Plots.savefig("plot_my_snpweights.png")

    #PyPlot.boxplot(my_snpweights)
    #Plots.savefig("plot253-my_snpweights.png")
    #PyPlot.plot(my_snpweights)
    #Plots.savefig("plot255-my_snpweights.png")

    npzwrite("data.npz",my_snpweights)

    save("data.jld", "data", my_snpweights)
    #load("data.jld")["data"]

    #my_snpweights = [ones(size(my_snpweights, 1)) my_snpweights]   # done later
    my_diagm_snpweights = diagm(my_snpweights[1,:])
    println()
    print("Gordon - size(my_snpweights) = ")
    println(size(my_snpweights))
    print("Gordon - my_diagm_snpweights[1,1] ,[2,2], [3,3] = ")
    println(my_diagm_snpweights[1,1], ",  ", my_diagm_snpweights[2,2], ",  ", my_diagm_snpweights[3,3])


    StatPlots.boxplot(my_snpweights[1,:])
    Plots.savefig("plot279_my_snpweights.png")
    StatPlots.boxplot(rand(100,10))
    StatPlots.boxplot(my_snpweights[1,:])
    StatPlots.plot(Normal(3,5),lw=3)
    show()
    println()
    println("SNP\t0\t1\t2\tmean\tsd\tcount\tMAF\tp\tq\t2√(pq)\tweight\tcountmap")
    for i = 1:10
        (μ_snp, σ_snp) = mean_and_std(snpmatrix[:,i])
        @printf("snp[%1d]:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t", i%10,
            (0 - μ_snp)/σ_snp, (1 - μ_snp)/σ_snp, (2 - μ_snp)/σ_snp, μ_snp, σ_snp)

        @printf("%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", sum(snpmatrix[:,i]), sum(snpmatrix[:,i]) / 2200,
                                                    my_snpweights_p[i], my_snpweights_q[i],
                                                    (2 * sqrt(my_snpweights_p .* (1 - my_snpweights_p)))[i],
                                                    my_diagm_snpweights[i,i])
        #print(countmap(snpmatrix[:,i]))
        print()
    end

    println()
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")

    println(RED*"HERE I ASSIGN THE WEIGHTS before the zscore !!!"*DEFAULT)
    my_snpmatrix_temp = deepcopy(my_snpmatrix)
    # NEXT LINE APPLIES THE SNPWEIGHTS JUST FINE, BUT THE DON'T STICK THROUGH STANDARDIZE
    # SO I HAD TO REAPPLY THEM
    A_mul_B!(my_snpmatrix, my_snpmatrix_temp, my_diagm_snpweights)# put weights on first to keep it sparce
    println("Here is the first my_snpmatrix with the weights in!!!")
    println("Gordon testing 130: size(my_snpmatrix) = $(size(my_snpmatrix))")
    println("Gordon testing 130: size(my_snpmatrix_temp) = $(size(my_snpmatrix_temp))")
    println("Gordon testing 130: typeof(my_snpmatrix) = $(typeof(my_snpmatrix))")
    println("Gordon testing 130: typeof(my_snpmatrix_temp) = $(typeof(my_snpmatrix_temp))")
    println("Gordon testing 130: peek at first 3 rows of my_snpmatrix BEFORE THE WEIGHTS")
    println(my_snpmatrix_temp[1,1:10])
    println(my_snpmatrix_temp[2,1:10])
    println(my_snpmatrix_temp[3,1:10])

    println("===============================================================")

    println("Gordon testing 130: peek at first 3 rows of my_snpmatrix AFTER THE WEIGHTS")
    println(my_snpmatrix[1,1:10])
    println(my_snpmatrix[2,1:10])
    println(my_snpmatrix[3,1:10])

    my_snpmatrix = StatsBase.zscore(my_snpmatrix, 1)    # standardize
    println(RED*"USE NEXT LINE TO DISABLE WEIGHTS from above !!!"*DEFAULT)
    my_snpmatrix = StatsBase.zscore(my_snpmatrix_temp, 1)    # USE THIS LINE TO DISABLE WEIGHTS from above !!!
    #copy!(my_snpmatrix_temp, my_snpmatrix)

    println("===============================================================")

    println("Gordon testing 130: peek at first 3 rows of my_snpmatrix AFTER ZSCORE BUT BEFORE THE WEIGHTS")
    println(my_snpmatrix[1,1:10])
    println(my_snpmatrix[2,1:10])
    println(my_snpmatrix[3,1:10])
    println("Gordon testing 130: describe(my_snpmatrix[:,2]) after zscore() but before the weights")
    describe(my_snpmatrix[:,2])

    println(RED*"HERE I ASSIGN THE WEIGHTS after the zscore !!!"*DEFAULT)
    my_snpmatrix_temp_after_zscore = deepcopy(my_snpmatrix)
    # I wish I could put the weights on before zscore to keep it sparce
    #   but zscore nuked the weights, which is its job
    A_mul_B!(my_snpmatrix, my_snpmatrix_temp_after_zscore, my_diagm_snpweights)


    println("===============================================================")

    println("Gordon testing 130: peek at first 3 rows of my_snpmatrix AFTER ZSCORE AND AFTER THE WEIGHTS")
    println(my_snpmatrix[1,1:10])
    println(my_snpmatrix[2,1:10])
    println(my_snpmatrix[3,1:10])
    println("Gordon testing 130: describe(my_snpmatrix[:,2]) after zscore() and after the weights")
    describe(my_snpmatrix[:,2])
    # NEXT LINE ENEDED UP CAUSING INDEX AND BROADCAST ERRORS LIKE my_snpmatrix[1,:]=5 AND v.xk[:, :] .= snpmatrix[:, v.idx]
    # PLUS I'M NOT SURE IT HELPED ANYWAY, SO GO BACK TO OLD STANDARIZE FUNCTION
    # *** turns out the problem was NaN's in my_snpmatrix so this StandardizedMatrix call could be good after all
    # *** but won't adding the column of 1's after we invoke it, break it ???
    #my_snpmatrix = StandardizedMatrix(my_snpmatrix)     # NEW WAY !!! THIS LINE IS GOOD, BUT I MOVED IT DOWN A PAGE
    my_snpmatrix = [ones(size(my_snpmatrix, 1)) my_snpmatrix]
    println()
    println("Gordon testing 130: size(my_snpmatrix) = $(size(my_snpmatrix)) after adding the column of ones.")


    println("===============================================================")

    println(RED*"SETUP after zscore TEST DATA HERE !!! "*DEFAULT)
    TEST_DATA = [3,5]
    if 1 in TEST_DATA           # Note: this test data DID crash the program at 26 iterations with a dimensions error
                                #       I didn't follow up on it becuase it worked great before and now I'm working
                                #       on adding the weights AFTER the zscore
        println("my_snpmatrix[:,2] = a_number")
        my_snpmatrix[:,2] = 5     #39 interations, just playing w/o weights
        my_snpmatrix[:,2] = 4     #27 interations
        my_snpmatrix[:,2] = 3     #22 interations
        my_snpmatrix[:,2] = 2     #16 interations
        my_snpmatrix[:,2] = 1     #10 interations and 2 betas I could see
        my_snpmatrix[:,2] = 0     #9  interations, identical results as bare toy data when I started
        my_snpmatrix[:,2] =-1     #10 interations, same result as positive 1
        println(my_snpmatrix[2,1:11])
        println(my_snpmatrix[3,1:11])
        println("\n")
    end

    fill!(v.xb, 0.0) #initialize β = 0 vector, so Xβ = 0
    println("Gordon testing 134: size(v.xb) = $(size(v.xb))")
    copy!(v.r, y)    #redisual = y-Xβ = y  CONSIDER BLASCOPY!
    v.r[mask_n .== 0] .= 0 #bit masking? idk why we need this yet
    println("Gordon testing 137: size(v.r) = $(size(v.r))")

    # calculate the gradient v.df = -X'(y - Xβ) = X'(-1*(Y-Xb)). Future gradient
    # calculations are done in iht!. Note the negative sign will be cancelled afterwards
    # when we do b+ = P_k( b - μ∇f(b)) = P_k( b + μ(-∇f(b))) = P_k( b + μ*v.df)

    # Can we use v.xk instead of snpmatrix?
    At_mul_B!(v.df, my_snpmatrix, v.r)
    print("Gordon testing 143: size(v.df), size(my_snpmatrix), size(v.r) = ")
    println(size(v.df), size(my_snpmatrix), size(v.r))

    println("===============================================================")
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")

    # snpmatrix = use_A2_as_minor_allele(x.snpmatrix) #to compare results with using PLINK
    snpmatrix = StatsBase.zscore(snpmatrix, 1)  # standardize
    #snpmatrix = StandardizedMatrix(snpmatrix)  # NEW WAY !!! doesn't work all the way
    snpmatrix = [ones(size(snpmatrix, 1)) snpmatrix]

    copy!(snpmatrix, my_snpmatrix)  # NEW WAY !!!
    println()
    println("Gordon 130: size(snpmatrix) = $(size(snpmatrix))")
    println("Gordon 130: contents of snpmatrix after zscore()")
    println(snpmatrix[1,1:10])
    println(snpmatrix[2,1:10])
    println(snpmatrix[3,1:10])
    println("\n")
    #
    # Begin IHT calculations
    #
    fill!(v.xb, 0.0) #initialize β = 0 vector, so Xβ = 0
    println("Gordon 134: size(v.xb) = $(size(v.xb))")
    copy!(v.r, y)    #redisual = y-Xβ = y  CONSIDER BLASCOPY!
    v.r[mask_n .== 0] .= 0 #bit masking? idk why we need this yet
    println("Gordon 137: size(v.r) = $(size(v.r))")

    # calculate the gradient v.df = -X'(y - Xβ) = X'(-1*(Y-Xb)). Future gradient
    # calculations are done in iht!. Note the negative sign will be cancelled afterwards
    # when we do b+ = P_k( b - μ∇f(b)) = P_k( b + μ(-∇f(b))) = P_k( b + μ*v.df)

    # Can we use v.xk instead of snpmatrix?
    At_mul_B!(v.df, snpmatrix, v.r)
    print("Gordon 143: size(v.df), size(snpmatrix), size(v.r) = ")
    println(size(v.df), size(snpmatrix), size(v.r))

    tic()

    for mm_iter = 1:max_iter
        println()
        print_with_color(:red,"=========== mm_iter LOOP ==========================================================\n")
        println("Gordon 145: mm_iter = $mm_iter")
        # save values from previous iterate
        copy!(v.b0, v.b)   # b0 = b   CONSIDER BLASCOPY!
        copy!(v.xb0, v.xb) # Xb0 = Xb  CONSIDER BLASCOPY!
        loss = next_loss

        #calculate the step size μ. Can we use v.xk instead of snpmatrix?
        println("Gordon 153: size(y) = $(size(y))")
        print_with_color(:red,"Gordon 153:, k = $k\n")
        (μ, μ_step) = iht!(v, snpmatrix, y, k, nstep=max_step, iter=mm_iter)

        # iht! gives us an updated x*b. Use it to recompute residuals and gradient
        v.r .= y .- v.xb
        v.r[mask_n .== 0] .= 0 #bit masking, idk why we need this yet

        At_mul_B!(v.df, snpmatrix, v.r) # v.df = X'(y - Xβ) Can we use v.xk instead of snpmatrix?
        print("Gordon 159: size(v.df), size(snpmatrix), size(v.r) = ")
        println(size(v.df), size(snpmatrix), size(v.r))

        # update loss, objective, gradient, and check objective is not NaN or Inf
        next_loss = sum(abs2, v.r) / 2
        !isnan(next_loss) || throw(error("Objective function is NaN, aborting..."))
        !isinf(next_loss) || throw(error("Objective function is Inf, aborting..."))

        # track convergence
        the_norm    = chebyshev(v.b, v.b0) #max(abs(x - y))
        print("Gordon 167: size(v.b), size(v.b0) = ")
        println(size(v.b), size(v.b0))
        scaled_norm = the_norm / (norm(v.b0, Inf) + 1)
        converged   = scaled_norm < tol

        if converged
            mm_time = toq()   # stop time
            #PyPlot.boxplot(rand(100,6))

            println("IHT converged in " * string(mm_iter) * " iterations")
            println("It took " * string(mm_time) * " seconds to converge")
            println("The estimated model is stored in 'beta'")
            println("There are " * string(countnz(v.b)) * " non-zero entries of β")
            return IHTResults(mm_time, next_loss, mm_iter, copy(v.b))
        end

        if mm_iter == max_iter
            mm_time = toq() # stop time
            throw(error("Did not converge!!!!! The run time for IHT was " *
                string(mm_time) * "seconds"))
        end
        print_with_color(:red,"=========== mm_iter LOOP END ====\n")
    end
end #function L0_reg

"""
Calculates the IHT step β+ = P_k(β - μ ∇f(β)).
Returns step size (μ), and number of times line search was done (μ_step).

This function updates: b, xb, xk, gk, xgk, idx
"""
function iht!(
    v         :: IHTVariable,
    snpmatrix :: Matrix{Float64},
    y         :: Vector{Float64},
    k         :: Int;
    iter      :: Int = 1,
    nstep     :: Int = 50,
)
    println()
    print_with_color(:green,"=========== BEGIN iht!() ==============================\n")
    println("Gordon 202: BEGIN iht!()")
    # compute indices of nonzeroes in beta and store them in v.idx (also sets size of v.gk)
    _iht_indices(v, k)

    # fill v.xk, which stores columns of snpmatrix corresponding to non-0's of b
    println("Gordon 206: ready to broadcast???")
    v.xk[:, :] .= snpmatrix[:, v.idx]
    print("Gordon 207: size(v.xk), size(snpmatrix), size(v.idx) = ")
    println(size(v.xk), size(snpmatrix), size(v.idx))

    # fill v.gk, which store only k largest components of gradient (v.df)
    # fill_perm!(v.gk, v.df, v.idx)  # gk = g[v.idx]
    v.gk .= v.df[v.idx]
    println("Gordon 211: size(v.gk) = $(size(v.gk))")

    # now compute X_k β_k and store result in v.xgk
    A_mul_B!(v.xgk, v.xk, v.gk)
    print("Gordon 214: size(v.xgk), size(v.xk), size(v.gk) = ")
    println(size(v.xgk), size(v.xk), size(v.gk))

    # warn if xgk only contains zeros
    all(v.xgk .≈ 0.0) && warn("Entire active set has values equal to 0")

    #compute step size and notify if step size too small
    μ = norm(v.gk, 2)^2 / norm(v.xgk, 2)^2
    isfinite(μ) || throw(error("Step size is not finite, is active set all zero?"))
    μ <= eps(typeof(μ)) && warn("Step size $(μ) is below machine precision, algorithm may not converge correctly")
    println("Gordon 222: before _iht_gradstep, μ = $μ")

    #Take the gradient step and compute ω. Note in order to compute ω, need β^{m+1} and xβ^{m+1} (eq5)
    _iht_gradstep(v, μ, k)
    #ω = compute_ω!(v, snpmatrix) #is snpmatrix required? Or can I just use v.x
    #these next 5 lines temporarily replace the previous line's call to compute_ω()
        A_mul_B!(v.xb, snpmatrix, v.b)
        print("Gordon 224: size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0) = ")
        println(size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0))
        ω = sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)
        println("Gordon 224: ω = $ω")
    println("Gordon 226: after  _iht_gradstep and compute_ω!(), ω = $ω")

    #compute ω and check if μ < ω. If not, do line search by halving μ and checking again.
    μ_step = 0
    for i = 1:nstep
        #exit loop if μ < ω where c = 0.01 for now
        if _iht_backtrack(v, ω, μ); break; end
        println()
        println("=========== BEGIN BACKTRACK LOOP ===============")
        println("Gordon 232: after _iht_backtrack(), i = $i")

        #if μ >= ω, step half and warn if μ falls below machine epsilon
        μ /= 2
        μ <= eps(typeof(μ)) && warn("Step size equals zero, algorithm may not converge correctly")
        println("Gordon 236: μ = $μ")

        # recompute gradient step
        copy!(v.b, v.b0)
        _iht_gradstep(v, μ, k)

        # re-compute ω based on xβ^{m+1}
        A_mul_B!(v.xb, snpmatrix, v.b)
        print("Gordon 244: size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0) = ")
        println(size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0))
        ω = sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)
        println("Gordon 244: ω = $ω")


        μ_step += 1
        println("=========== END BACKTRACK LOOP =====")
    end
    print_with_color(:green,"=========== END iht!() =====\n")
    return (μ, μ_step)
end
