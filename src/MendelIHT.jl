using MendelBase
using SnpArrays
using IHT
using DataFrames
using Distances
using StatsBase

#for Gordon
#using PyPlot
#using Plots
using StatPlots
#using BenchmarkTools
#using Distributions
using JLD
using NPZ
#export MendelIHT
#using RDatasets
#using Gadfly

const RED     = "\e[1;31m" # BOLD
const GREEN   = "\x1b[32m"
const YELLOW  = "\x1b[33m"
const BLUE    = "\e[1;34m" # BOLD
const PURPLE  = "\x1b[35m"
const BOLD    = "\x1b[1m"
const DEFAULT = "\x1b[0m"
const YELLOW_BG = "\e[43m"
const CYAN_BG = "\e[46m"
#= Background
Value	Color
\e[40m	Black
\e[41m	Red
\e[42m	Green
\e[43m	Yellow
\e[44m	Blue
\e[45m	Purple
\e[46m	Cyan
\e[47m	White
=#
USE_INTERCEPT = true
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
    #k = 9 # DEBUG NO Intercept
    result = L0_reg(snpdata, phenotype, k)
    println("===============================================================")
    println("============= MORE INFORMATION ABOUT SELECTED BETAS ===========")
    println("============= IN MENDELIHT ====================================")
    println("===============================================================")
    println("time returned is $(result.time)")
    found = find(result.beta .!= 0.0)
    print("found: ")
    println(found)
    println("result.beta, the betas(coefficients), are the result of the IHT")
    print("result.beta[found]: ")
    println(result.beta[found])
    println("result.beta[1:10] $(size(result.beta)) = $(result.beta[1:10])")
    describe(result.beta)
    snp = 1
    println("SNP[$(snp)] SNP        = $(snp_definition_frame[snp,:SNP])")
    println("SNP[$(snp)] Chromosome = $(snp_definition_frame[snp,:Chromosome])")
    println("SNP[$(snp)] BasePairs  = $(snp_definition_frame[snp,:Basepairs])")
    println("SNP[$(snp)] Allele1    = $(snp_definition_frame[snp,:Allele1])")
    println("SNP[$(snp)] Allele2    = $(snp_definition_frame[snp,:Allele2])")
    println()
if true        # NO NEED FOR SNP NAMES SINCE THEY ARE FAKE? or are they BMI related
    run(`powershell pwd`)   # works !!!
    rs = Nullable{Int64}[]
    for snp_index in found
        if USE_INTERCEPT && false # intercept is at the end now
            if snp_index == 1
                continue
            end
            snp = snp_index - 1 # adjust for intercept
        else
            if snp_index == 10001
                continue
            end
            snp = snp_index
        end

        #aaaa=readstring(`bash -c 'curl -X GET --header "Accept: application/json" "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/268"'`);
        #println("back from readstring()")
        #run(`bash -c curl -X GET --header "Accept: application/json" "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/$(snp)"` |> "$(snp).txt")
        r = tryparse(Int64, snp_definition_frame[snp,:SNP][3:end])
        push!(rs, r)
        println("SNP[$(snp)] = $(snp_definition_frame[snp,:SNP])")
        println("SNP[$(snp)] = $(snp_definition_frame[snp,:Chromosome])")
        println("SNP[$(snp)] = $(snp_definition_frame[snp,:Basepairs])")
        println("SNP[$(snp)] = $(snp_definition_frame[snp,:Allele1])")
        println("SNP[$(snp)] = $(snp_definition_frame[snp,:Allele2])")
        println()
    end
    println("SNP names to pass to NCBI:")
    println(rs)
end
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


    # snpmatrix = use_A2_as_minor_allele(x.snpmatrix) #to compare results with using PLINK
    snpmatrix = StatsBase.zscore(snpmatrix, 1)
    println("Gordon testing 130,510: Set Intercept here or NOT, and in MendelIHT_utilities.jl")
    println("Note: Set USE_INTERCEPT = false to drop intercept column.")
    #USE_INTERCEPT = true  # its global now
    if USE_INTERCEPT # && false
        snpmatrix = [snpmatrix ones(size(snpmatrix, 1)) ]
    end
    println("Note: Set USE_WEIGHTS = true to use weights.")
    USE_WEIGHTS = true
    #USE_WEIGHTS = false
    if USE_WEIGHTS
        hold_snpmatrix = deepcopy(snpmatrix)
        my_snpMAF, my_snpweights, snpmatrix = calculatePriorWeightsforIHT(x,y,k,v)
        snpmatrix = deepcopy(hold_snpmatrix) # cancel snpweighting here, but keep my_snpweights for later
    end
    if USE_INTERCEPT # && false
        # this puts the intercept on the end at 10,001 - I like it there better anyway
        hold_my_snpweights = deepcopy(my_snpweights)
        my_snpweights = [ones(size(my_snpweights, 1)) my_snpweights]   # done later
        my_diagm_snpweights = diagm(my_snpweights[1,:]) # Diagonal(my_snpweights)
        my_snpweights = deepcopy(hold_my_snpweights) # no one else wants to see an intercept in the snpweights
    else
        my_diagm_snpweights = diagm(my_snpweights[1,:]) # Diagonal(my_snpweights)
    end
    # Use the next 2 lines to put the weights back in here if needed
    # but for now I want to try applying the weight after the _iht_gradstep()
    hold_snpmatrix = deepcopy(snpmatrix)
    A_mul_B!(snpmatrix, hold_snpmatrix, my_diagm_snpweights)

    #
    # Begin IHT calculations
    #
    fill!(v.xb, 0.0) #initialize β = 0 vector, so Xβ = 0
    copy!(v.r, y)    #redisual = y-Xβ = y  CONSIDER BLASCOPY!
    v.r[mask_n .== 0] .= 0 #bit masking? idk why we need this yet

    # calculate the gradient v.df = -X'(y - Xβ) = X'(-1*(Y-Xb)). Future gradient
    # calculations are done in iht!. Note the negative sign will be cancelled afterwards
    # when we do b+ = P_k( b - μ∇f(b)) = P_k( b + μ(-∇f(b))) = P_k( b + μ*v.df)

    # Can we use v.xk instead of snpmatrix?
    At_mul_B!(v.df, snpmatrix, v.r)

    tic()   # duplicate clock here to skip timing debug statements in weight function

    for mm_iter = 1:max_iter
        print_with_color(:red,"\n=========== mm_iter LOOP (mm_iter = $(mm_iter)) ===================================================\n")
        # save values from previous iterate
        copy!(v.b0, v.b)   # b0 = b   CONSIDER BLASCOPY!
        copy!(v.xb0, v.xb) # Xb0 = Xb  CONSIDER BLASCOPY!
        loss = next_loss

        #calculate the step size μ. Can we use v.xk instead of snpmatrix?
        (μ, μ_step) = iht!(v, snpmatrix, y, k, nstep=max_step, iter=mm_iter)

        # iht! gives us an updated x*b. Use it to recompute residuals and gradient
        v.r .= y .- v.xb
        v.r[mask_n .== 0] .= 0 #bit masking, idk why we need this yet

        At_mul_B!(v.df, snpmatrix, v.r) # v.df = X'(y - Xβ) Can we use v.xk instead of snpmatrix?

        # update loss, objective, gradient, and check objective is not NaN or Inf
        next_loss = sum(abs2, v.r) / 2
        !isnan(next_loss) || throw(error("Objective function is NaN, aborting..."))
        !isinf(next_loss) || throw(error("Objective function is Inf, aborting..."))

        # track convergence
        the_norm    = chebyshev(v.b, v.b0) #max(abs(x - y))
        scaled_norm = the_norm / (norm(v.b0, Inf) + 1)
        converged   = scaled_norm < tol

        if converged
            mm_time = toq()   # stop time

            println("===============================================================")
            println("============= INFORMATION ABOUT SELECTED BETAS ================")
            println("============= IN L0_reg =======================================")
            println("===============================================================")

            println("IHT converged in " * string(mm_iter) * " iterations")
            println("It took " * string(mm_time) * " seconds to converge")
            println("The estimated model is stored in 'beta'")
            println("There are " * string(countnz(v.b)) * " non-zero entries of β")
            found = find(v.b .!= 0.0)
            print("found: ")
            println(found)
            println("v.b, the betas(coefficients), are the result of the IHT")
            print("v.b[found]: ")
            println(v.b[found])
            println("v.b[1:10] $(size(v.b)) = $(v.b[1:10])")
            describe(v.b)

            println(BLUE*"THE SUMMARY STATS FOR THE Betas ARE PLOTTED IN bar_b.png"*DEFAULT)
            StatPlots.bar(v.b)
            Plots.savefig("bar_b1.png")

            println()
            if USE_INTERCEPT # && false
                my_snpweights_intercept = [ones(size(my_snpweights, 1)) my_snpweights]
                println("my_snpweights_intercept for found:")
                println(size(my_snpweights_intercept))
                println(my_snpweights_intercept[1,found])

                println()
                my_snpMAF_intercept = [ones(size(my_snpMAF, 1)) my_snpMAF]
                println("my_snpMAF_intercept for found:")
                println(size(my_snpMAF_intercept))
                println(my_snpMAF_intercept[1,found])
            else
                println("my_snpweights for found:")
                println(size(my_snpweights))
                println(my_snpweights[1,found])

                println()
                println("my_snpMAF for found:")
                println(size(my_snpMAF))
                println(my_snpMAF[1,found])
            end
            println()
            println("y, the phenotype for 2200 people, is given in the file 'gwas 1 data.fam'")
            print("y: ")
            print(size(y))
            println(y[1:10])
            #println(y[found])
            describe(y)
            println(BLUE*"THE Phenotype for 2200 people are plotted IN bar_phenotype.png"*DEFAULT)
            StatPlots.bar(y)
            Plots.savefig("bar_phenotype.png")


            println()
            println("snpmatrix holds allele counts, zscored, and weighted for 2200 people")
            println("here is a sample of the first 10 SNPs for the first 3 people")
            println("note: first SNP is likely 1 for the intercept")
            println("snpmatrix[1,1:10] = $(snpmatrix[1,1:10])")
            println("snpmatrix[2,1:10] = $(snpmatrix[2,1:10])")
            println("snpmatrix[3,1:10] = $(snpmatrix[3,1:10])")

            println()
            println("v.xb is x*b for 2200 people, where x is weighted_snpmatrix and b are the coefficients v.b")
            print("v.xb: ")
            print(size(v.xb))
            println("v.xb[1:10] = $(v.xb[1:10])")
            println("type(v.xb) = $(typeof(v.xb))")
            describe(v.xb)

            #=
            println()
            print("(y .- v.xb): ")      # same as v.r
            println(size(y .- v.xb))
            describe((y .- v.xb))
            #println((y .- v.xb)[found])
            =#
            println()
            println("v.r are the residuals for 2200 data points (people)")
            println("we are trying to minimize a function of this number")
            println("v.r = y .- v.xb")
            println("the function is 'next_loss = sum(abs2, v.r) / 2'")
            print("v.r: ")
            println(size(v.r))
            describe(v.r)
            #println(v.r[found])
            println(BLUE*"THE SUMMARY STATS FOR THE Residuals ARE PLOTTED IN bar_r1.png"*DEFAULT)
            StatPlots.bar(v.r)
            Plots.savefig("bar_r1.png")

            println("THE REGRESSION RESULTS ARE SAVED IN bar_b.png")
            StatPlots.bar(found,v.b[found])
            Plots.savefig("bar_b2.png")
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
    print_with_color(:green,"\n=========== BEGIN iht!() ==============================\n")
    # compute indices of nonzeroes in beta and store them in v.idx (also sets size of v.gk)
    _iht_indices(v, k)

    # fill v.xk, which stores columns of snpmatrix corresponding to non-0's of b
    println("Gordon 206: ready to broadcast???")
    v.xk[:, :] .= snpmatrix[:, v.idx]

    # fill v.gk, which store only k largest components of gradient (v.df)
    # fill_perm!(v.gk, v.df, v.idx)  # gk = g[v.idx]
    v.gk .= v.df[v.idx]

    # now compute X_k β_k and store result in v.xgk
    A_mul_B!(v.xgk, v.xk, v.gk)

    # warn if xgk only contains zeros
    all(v.xgk .≈ 0.0) && warn("Entire active set has values equal to 0")

    #compute step size and notify if step size too small
    μ = norm(v.gk, 2)^2 / norm(v.xgk, 2)^2
    isfinite(μ) || throw(error("Step size is not finite, is active set all zero?"))
    μ <= eps(typeof(μ)) && warn("Step size $(μ) is below machine precision, algorithm may not converge correctly")
    println("Gordon 222: μ = $(μ), before _iht_gradstep")

    #Take the gradient step and compute ω. Note in order to compute ω, need β^{m+1} and xβ^{m+1} (eq5)
    _iht_gradstep(v, μ, k)
    #ω = compute_ω!(v, snpmatrix) #is snpmatrix required? Or can I just use v.x
    #these next 10 lines temporarily replace the previous line's call to compute_ω()
#A_mul_B!(Y, A, B) → Y
#Calculates the matrix-matrix or matrix-vector product A⋅B and stores the result in Y,
    A_mul_B!(v.xb, snpmatrix, v.b)
    found = find(v.b .!= 0.0)
    print("found: ")
    println(found)
    println("v.b, the betas(coefficients), are the result of the IHT")
    print("v.b[found]: ")
    println(v.b[found])
    print("v.xb: ")
    describe(v.xb)
    ω = sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)
    println("Gordon 226: μ = $(μ), after _iht_gradstep")
    println("Gordon 226: ω = $(ω), after  _iht_gradstep and compute_ω!()")
    print("Gordon 226: size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0) = ")
    println(size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0))

    #compute ω and check if μ < ω. If not, do line search by halving μ and checking again.
    μ_step = 0
    for i = 1:nstep
        #exit loop if μ < ω where c = 0.01 for now
        if _iht_backtrack(v, ω, μ); break; end
        println("\n=========== BEGIN BACKTRACK LOOP (i = $i) ===============")
        #println("Gordon 232: after _iht_backtrack(), i = $i")

        #if μ >= ω, step half and warn if μ falls below machine epsilon
        μ /= 2
        μ <= eps(typeof(μ)) && warn("Step size equals zero, algorithm may not converge correctly")
        println("Gordon 236: μ = $μ")

        # recompute gradient step
        copy!(v.b, v.b0)
        _iht_gradstep(v, μ, k)

#A_mul_B!(Y, A, B) → Y
#Calculates the matrix-matrix or matrix-vector product A⋅B and stores the result in Y,
        # re-compute ω based on xβ^{m+1}
        A_mul_B!(v.xb, snpmatrix, v.b)
        #print("Gordon 244: size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0) = ")
        #println(size(v.xb), size(snpmatrix), size(v.b), size(v.b0), size(v.xb0))
        ω = sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)
        println("Gordon 244: ω = $ω")


        μ_step += 1
        println("=========== END BACKTRACK LOOP =====")
    end

    print_with_color(:green,"=========== END iht!() =====\n")
    return (μ, μ_step)
end

"""
Calculates the Prior Weighting for IHT.
Returns a weight matrix (snpweights)
    and a weighted SNP matrix (snpmatrix).

This function updates: hopefully nothing???
"""
function calculatePriorWeightsforIHT(
    x        :: SnpData,
    y        :: Vector{Float64},
    k        :: Int,
    v        :: IHTVariable
)
    #convert bitarrays to Float64 genotype matrix, standardize each SNP, and add intercept
    snpmatrix = convert(Array{Float64,2}, x.snpmatrix)
    println("===============================================================")
    println("============= BEGIN CODE FOR PRIOR WEIGHTS FUNCTION ===========")
    println("===============================================================")
    #myplot = Gadfly.plot(dataset("HistData", "ChestSizes"), x="Chest", y="Count", Geom.bar)
    #myplot = Gadfly.plot(x = 1:10, y = 2:11, Geom.line)
    #draw(PNG("first_gadfly_bar_chart.png", 3inch, 3inch), myplot)

    println(BLUE*"Gordon 124: new code for Prior Weights work"*DEFAULT)
    println("Note: numbers after my name should be close to line numbers in original code.")
    USE_WEIGHTS = true
    println("Note: Set USE_WEIGHTS = true to use weights.")
    #USE_INTERCEPT = true  # its global now
    println("Note: Set USE_INTERCEPT = false to drop intercept column.")
    println("Note:      must adjust intercept in MendelIHT_utilities.jl ALSO!")

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

    print(YELLOW_BG)
    println("===============================================================")
    println("========== PRIOR WEIGHTS - BEGIN CALCULATIONS =================")
    println("===============================================================")
    print(DEFAULT)
    # Create my own snpmatrix for Prior Weights work
    my_snpmatrix = zeros(snpmatrix)
    copy!(my_snpmatrix, snpmatrix)
    my_snpmatrix = deepcopy(snpmatrix)  # ask Ben about copying matricies
    println(my_snpmatrix[1,1:20])       # all 0,1,2 as Ben said

    # ALLELE_MAX is times 2 * SNPcount because each PERSON's SNP has 2 alleles
    ALLELE_MAX = 2 * size(my_snpmatrix,1)
    println("ALLELE_MAX = $(ALLELE_MAX)")
    # MAF = Minor Allele Frequency
    # my_snpMAF is the sum the allele counts 0,1, or 2 for each SNP (column)
    #   in the snpmatrix / ALLELE_MAX
    # Note: Since it is a MINOR allele the max should be 0.500000, OK!

    #= development code, no longer needed OK to delete
    println("======== SIDE STEP - Look for outliers ========================")
    #  look for outliers
    my_snpsum = sum(my_snpmatrix, 1)
    println(find(my_snpsum .< 110))
    print("Gordon - size(my_snpsum) = ")
    println(size(my_snpsum))
    my_snpsum = sort(my_snpsum[1,:])
    print("Gordon - my_snpsum = ")
    println(my_snpsum[1:100])# / ALLELE_MAX)       #
    print("Gordon - max(my_snpsum) = ")
    println(maximum(my_snpsum))
    println("Gordon - describe(my_snpsum[1,:] / ALLELE_MAX) - need for making my_snpweights_p")
    describe(my_snpsum / ALLELE_MAX)
    my_snpMAF = sum(my_snpmatrix, 1) / ALLELE_MAX
    println(find(my_snpMAF .< 0.025))
    =#

    # Can't kill the outliers here, weighting will just try to bring them back, see outlier code below

    # Minor Allele Freq (MAF) = sum of alleles / (2 * SNPcount)
    my_snpMAF = sum(my_snpmatrix, 1) / ALLELE_MAX
    print("Gordon - size(my_snpMAF) = ")
    println(size(my_snpMAF))
    print("Gordon - my_snpMAF = ")
    println(my_snpMAF[1,1:10])
    print("Gordon - max(my_snpMAF) = ")
    println(maximum(my_snpMAF))
    describe(my_snpMAF[1,:])


    # GORDON - CALCUATE CONSTANT WEIGHTS - another weighting option
    my_snpweights_const = copy(my_snpMAF) # only to allocate my_snpweights_const
    for i = 1:size(my_snpweights_const,2)
        my_snpweights_const[1,i] = 1.0   # was 0.5
    end
    println(my_snpweights_const[1,1:10])
    println("\ndescribe(my_snpweights_const)")
    describe(my_snpweights_const[1,:])

    # GORDON - CALCULATE WEIGHTS BASED ON p=MAF, 1/(2√pq) SUGGESTED BY BEN AND HUA ZHOU
    my_snpweights_p = my_snpMAF      # p_hat
    #my_snpweights_p /= 2    # NOTE: OK! THIS IS BECAUSE 1 PERSON = 2 chromosomes
                            # 2019-0816 fixed ALLELE_COUNT to ALLELE_MAX
                            # a bug fix for p is MAF s/b max .5 per Ben
                            # also is was getting some NaN in my_diagm_snpweights
                            # because of some inf in my_snpweights
                            # that I imagine came from some √(0.0*(1.0-1.0)) errors
    my_snpweights_q = 1 - my_snpweights_p   # q_hat
    my_snpweights = my_snpweights_p + my_snpweights_q   # just verifying p + q == 1 OK!
    my_snpweights_pq = my_snpweights_p .* my_snpweights_q   # just verifying P .* q max is 0.25 OK!
    # next line is 1/pq
    my_snpweights = my_snpweights .\ 1      # not what we want, but this works! to get invervse of each element
    my_snpweights = sqrt(my_snpweights_p .* my_snpweights_q)   # just verifying sqrtm(p .* q) == 0.5 OK!
    # next line is 2 * sqrt(p(1-p)) from Hua Zhou 2011 paper pg. 110
    my_snpweights = 2 * sqrt(my_snpweights_p .* (1 - my_snpweights_p))   # just verifying 2 * sqrtm(p .* q) == 1.0 OK!
    my_snpweights_huazhou = my_snpweights
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
    println("\ndescribe(my_snpweights 2√pq)")
    describe(my_snpweights[1,:])
    my_snpweights = my_snpweights .\ 1      # this works! to get reciprocal of each element
    my_snpweights_huazhou_reciprocal = my_snpweights
    println("\ndescribe(my_snpweights 2√pq .\\ 1)")
    describe(my_snpweights[1,:])
    # THESE COMMENTED DESCRIPTIONS ARE FROM BEFORE I CUT P IN HALF !!!
    # Mean:           Inf
    # Minimum:        1.000000
    # 1st Quartile:   1.026816
    # Median:         1.124446
    # 3rd Quartile:   1.367013
    # Maximum:        Inf

    print(CYAN_BG)
    println("===============================================================")
    println("=========== SELECT WEIGHT FUNCTION HERE =======================")
    println("===============================================================")
    print(DEFAULT)
    # DECIDE NOW WHICH WEIGHTS TO APPLY !!!

    if true
        my_snpweights = copy(my_snpweights_const)    # Ben/Kevin this is currently at 1.0 for testing null effect
        println("Selected my_snpweights_const")
    end
    if false
        my_snpweights = copy(my_snpweights_huazhou_reciprocal)
        println("Selected my_snpweights_huazhou_reciprocal")
    end

    #my_snpweights = copy(my_snpMAF)
    println(RED*"SELECT WEIGHT FUNCTION HERE !!!"*DEFAULT)
    println("\ndescribe(my_snpweights)")
    describe(my_snpweights[1,:])

    println("===============================================================")
    println("========== REMOVE OUTLIERS HERE ===============================")
    println("========== AND TRIM DATA POINTS ===============================")
    println("===============================================================")
    #   DO I REALLY WANT TO KILL OUTLIERS HERE - I LIKE THEM, MUST BE A BETTER WAY TO KEEP THEM
    NULLweight = 1 # was 0 or 1 (Ben says don't weight to zero, one is good enough)
    cutoff = 0.025
    cutoff = 0.001
    #cutoff = 0.025
    found = find(my_snpMAF .< cutoff)    # below 2.5% causes me problems with convergence
    println(found)
    println(my_snpMAF[found])
    println(RED*"Setting weight = $(NULLweight) for $(size(found)) outliers with MAF below $(cutoff) cutoff."*DEFAULT)
    my_snpweights[found] = NULLweight
    println(my_snpweights[found])
    println("Here are the rest of the SNPs with MAF < .05")
    found = find(my_snpMAF .< 0.05)    # below 2.5% causes me problems with convergence
    println(found)
    println(my_snpMAF[found])
    println(RED*"NOT Setting weight = $(NULLweight) for $(size(found)) RARE KEEPERS with MAF below $(0.05) cutoff."*DEFAULT)
    println(my_snpweights[found])

    # trim the top too ???
    TRIM_TOP = 0.1
    TRIM_TOP = 0.0
    if TRIM_TOP > 0.0
        cutoff_top = TRIM_TOP
        found_top = find(my_snpMAF .> cutoff_top)    # below 2.5% causes me problems with convergence
        #println(found_top)
        #println(my_snpMAF[found_top])
        println(RED*"Setting weight = 0 for $(size(found_top)) data points with MAF above $(cutoff_top) cutoff."*DEFAULT)
        my_snpweights[found_top] = 0.0
        #println(my_snpweights[found_top])
    end

    println("===============================================================")
    println("===============================================================")
    println("====== APPLYING THE APRIORI WEIGHTS AS IF BY APRIORI ==========")
    println("===============================================================")
    println(RED*"HERE I ASSIGN THE APRIORI WEIGHTS AS IF BY APRIORI !!!"*DEFAULT)
    # APRIORI
    # DEVELOPERS NOTE: TO REALLY SEE THE EFFECT OF Weighting
    #   TURN APRIORI ON AND SET WEIGHT ALGO TO CONSTANT = 1.0
    # [899, 1881, 3775, 3982, 4210, 4399, 5794, 6612, 6628, 7024] # 10 here, 12 below
    if 99 in [98]
        #for m2 in [899, 1881, 3775, 3982, 4210, 4399, 5794, 6612, 6628, 7024, 8960, 9468]
        for m2 in [898, 1880, 2111, 3774, 3981, 4398, 5793, 6611, 6627, 7023] # [898, 1880, 2111, 3774, 3981, 4398, 6611, 6627, 7023]
            m1 = 0.8
            println("Apply adjustment for APRIORI weight at m2 = $(m2), m1 = $(m1)")
            print("Before = $(my_snpweights[m2-1])")
            my_snpweights[m2-1] = m1
            println(", after = $(my_snpweights[m2-1])")
        end
        for m2 in [100:500...] # [898, 1880, 2111, 3774, 3981, 4398, 6611, 6627, 7023]
            m1 = 1.59
            println("Apply adjustment for APRIORI weight at m2 = $(m2), m1 = $(m1)")
            print("Before = $(my_snpweights[m2-1])")
            my_snpweights[m2-1] = m1
            println(", after = $(my_snpweights[m2-1])")
        end
        for m2 in [1100:1500...] # [898, 1880, 2111, 3774, 3981, 4398, 6611, 6627, 7023]
            m1 = 1.61
            println("Apply adjustment for APRIORI weight at m2 = $(m2), m1 = $(m1)")
            print("Before = $(my_snpweights[m2-1])")
            my_snpweights[m2-1] = m1
            println(", after = $(my_snpweights[m2-1])")
        end
    else
        println("None")
    end
    #    my_diagm_snpweights[7024-1,7024-1] = 2.0 # BEN'S IDEA FOR LOWEST LOSS

    println("===============================================================")
    println("============ MAKE GRAPHS OF SELECTED WEIGHT FUNCTION ==========")
    println("===============================================================")

    println(BLUE*"THE SUMMARY STATS FOR THE WEIGHTS ARE PLOTTED IN histogram_my_snpweights.png"*DEFAULT)
    StatPlots.histogram(my_snpweights[1,:])
    Plots.savefig("histogram_my_snpweights.png")


    println(BLUE*"THE SUMMARY STATS FOR THE WEIGHTS ARE PLOTTED IN histogram_my_snpweights_p.png"*DEFAULT)
    StatPlots.histogram(my_snpweights_p[1,:])
    Plots.savefig("histogram_my_snpweights_p.png")


    println(BLUE*"THE SUMMARY STATS FOR THE WEIGHTS ARE PLOTTED IN bar_my_snpweights.png"*DEFAULT)
    StatPlots.bar(my_snpweights[1,:])
    Plots.savefig("bar_my_snpweights.png")

    println(BLUE*"THE SUMMARY STATS FOR THE WEIGHTS ARE PLOTTED IN boxplot_my_snpweights.png"*DEFAULT)
    StatPlots.boxplot(my_snpweights[1,:])
    Plots.savefig("boxplot_my_snpweights.png")

    println()

    # =================================================================
    # =================================================================
    println(BLUE*"HERE ARE THE SUMMARY STATS FOR THE WEIGHTS"*DEFAULT)
    print("Gordon - size(my_snpweights) = ")    # (1, 10000)
    println(size(my_snpweights))
    print("Gordon - my_snpweights = ")
    println(my_snpweights[1,1:10])
    print("Gordon - max(my_snpweights) = ")
    println(maximum(my_snpweights))
    describe(my_snpweights[1,:])

    # save the variables to disk for external analysis
    # specifically intended for making inline plots in Jupyter
    npzwrite("my_snpweights.npz",my_snpweights)
    npzwrite("my_snpweights_p.npz",my_snpweights_p)
    npzwrite("my_snpweights_q.npz",my_snpweights_q)
    npzwrite("my_snpweights_pq.npz",my_snpweights_pq)
    npzwrite("my_snpweights_huazhou.npz",my_snpweights_huazhou)
    npzwrite("my_snpweights_huazhou_reciprocal.npz",my_snpweights_huazhou_reciprocal)


    save("my_snpweights.jld", "data", my_snpweights)
    #load("data.jld")["data"]



    println("===============================================================")
    println("============ MAKE THE WEIGHTS MATRIX DIAGONAL =================")
    println("===============================================================")
    #my_snpweights = [ones(size(my_snpweights, 1)) my_snpweights]   # done later
    my_diagm_snpweights = diagm(my_snpweights[1,:]) # Diagonal(my_snpweights)
    println()
    print("Gordon - size(my_snpweights) = ")
    println(size(my_snpweights))
    print("Gordon - size(my_diagm_snpweights) = ")
    println(size(my_diagm_snpweights))
    print("Gordon - my_diagm_snpweights[1,1] ,[2,2], [3,3] = ")
    println(my_diagm_snpweights[1,1], ",  ", my_diagm_snpweights[2,2], ",  ", my_diagm_snpweights[3,3])


    println("===============================================================")
    println("============ MAKE A CHART OF THE WEIGHT FUNCTION VARIABLES ====")
    println("===============================================================")
    println()
    println("SNP\t0\t1\t2\tmean\tsd\tcount\tMAF\tp\tq\t2√(pq)\tweight\tcountmap")
    for i = 1:10
        (μ_snp, σ_snp) = mean_and_std(snpmatrix[:,i])
        @printf("snp[%1d]:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t", i%10,
            (0 - μ_snp)/σ_snp, (1 - μ_snp)/σ_snp, (2 - μ_snp)/σ_snp, μ_snp, σ_snp)

        @printf("%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", sum(snpmatrix[:,i]), my_snpMAF[1,i],
                                                    my_snpweights_p[i], my_snpweights_q[i],
                                                    (2 * sqrt(my_snpweights_p .* (1 - my_snpweights_p)))[i],
                                                    my_diagm_snpweights[i,i])
        #print(countmap(snpmatrix[:,i]))
        print()
    end
    DEBUG_LEVEL = [98]
    if 99 in DEBUG_LEVEL    # special MOVE APRIORI winners up front for easier visualization
        my_snpmatrix[:,1] = my_snpmatrix[:,898]
        my_snpmatrix[:,2] = my_snpmatrix[:,1880]
        my_snpmatrix[:,3] = my_snpmatrix[:,3981]
        my_snpmatrix[:,4] = my_snpmatrix[:,4398]
        my_snpmatrix[:,5] = my_snpmatrix[:,5793]
        my_snpmatrix[:,6] = my_snpmatrix[:,6611]
        my_snpmatrix[:,7] = my_snpmatrix[:,6627]
        my_snpmatrix[:,8] = my_snpmatrix[:,7023]
        my_snpmatrix[:,9] = my_snpmatrix[:,9467]
    end
    my_snpmatrix_temp = deepcopy(my_snpmatrix) # protect my_snpmatrix for use below
    #copy!(my_snpmatrix_temp, my_snpmatrix)

    println("===============================================================")
    println("===============================================================")
    println("========== NOW APPLY THE ZSCORE TO SNPMATRIX ==================")
    println("===============================================================")
    my_snpmatrix = StatsBase.zscore(my_snpmatrix, 1)

    println("===============================================================")
    println("Gordon testing 130: peek at first 3 rows of my_snpmatrix AFTER ZSCORE BUT BEFORE THE WEIGHTS")
    println(my_snpmatrix[1,1:10])
    println(my_snpmatrix[2,1:10])
    println(my_snpmatrix[3,1:10])
    println("Gordon testing 130: describe(my_snpmatrix[:,2]) after zscore() but before the weights")
    describe(my_snpmatrix[:,2])

    println("===============================================================")
    println("===============================================================")
    println("========== APPLYING THE WEIGHTS AFTER ZSCORE ==================")
    println("===============================================================")
    println(RED*"HERE I ASSIGN THE WEIGHTS after the zscore !!!"*DEFAULT)
    my_snpmatrix_temp_after_zscore = deepcopy(my_snpmatrix)
    # I wish I could put the weights on before zscore to keep it sparce
    #   but zscore nuked the weights, which is its job
#    my_diagm_snpweights[7024-1,7024-1] = 2.0 # BEN'S IDEA FOR LOWEST LOSS
    A_mul_B!(my_snpmatrix, my_snpmatrix_temp_after_zscore, my_diagm_snpweights)

    println("===============================================================")
    println("===============================================================")
    println("========== CHANGING THE ACTUAL SNPs based on APRIORI  =========")
    println("===============================================================")
    println(RED*"HERE I CHANGE THE ACTUAL SNPs based on APRIORI !!!"*DEFAULT)
    # APRIORI
    # [899, 1881, 3775, 3982, 4210, 4399, 5794, 6612, 6628, 7024] # 10 here, 12 below
    describe(my_snpmatrix[:,899-1])
    describe(my_snpmatrix[:,3775-1])
    describe(my_snpmatrix[:,7024-1])
    #m1 = median(my_snpmatrix[:,7024-1])
    #m1 = maximum(my_snpmatrix[:,7024-1])
    #m1 = mean(my_snpmatrix[:,7024-1])
    #println("m1 = $(m1)")
    #my_snpmatrix[:,7024-1] = m1
    if 99 in []
        for m2 in [899, 1881, 3775, 3982, 4210, 4399, 5794, 6612, 6628, 7024, 8960, 9468]
            m1 = maximum(my_snpmatrix[:,m2-1])
            m1 = m1 * 0.8
            println("Apply adjustment to ACTUAL SNPs based on APRIORI m2 = $(m2), m1 = $(m1)")
            println("duplicates: $(find(my_snpweights .== m1))")  # matches something in m2, didn't find any
            println("snpweight = $(my_snpweights_p[m2-1])")
            my_snpmatrix[:,m2-1] = m1
            print("Before = 'weighted value for each person'")
            println(", after = 'constant value $(m1) for every person'")
        end
    else
        println("None")
    end

    # look for duplicates (put all the maximums in a list and sort it)
    println(find(my_snpweights .== maximum(my_snpmatrix[:,7024-1])))  # matches 7024
    println(maximum(my_snpmatrix))
    list1 = Float64[] # => 0-element Float64 Array
    nonzcount = zcount = 0
    for m3 in 1:size(my_snpmatrix,2)
        m4 = maximum(my_snpmatrix[:,m3])
        if m4 in list1
            my_snpmatrix[:,m3] = 0  # 0 out all dupicate records
            zcount += 1
        else
            push!(list1, m4)
            nonzcount += 1
        end
    end
    println("zcount = $(zcount)")       #  599
    println("nonzcount = $(nonzcount)") # 9401
    sort!(list1)
    #println(list1) # long list 9401 items
    #=
    duplicated(x) = foldl(
      (d,y)->(x[y,:] in d[1] ? (d[1],push!(d[2],y)) : (push!(d[1],x[y,:]),d[2])),
      (Set(), Vector{Int}()),
      1:size(x,1))[2]
    =#
    println("===============================================================")
    println("Gordon testing 130: peek at first 3 rows of my_snpmatrix AFTER ZSCORE AND AFTER THE WEIGHTS")
    println(my_snpmatrix[1,1:10])
    println(my_snpmatrix[2,1:10])
    println(my_snpmatrix[3,1:10])
    println("Gordon testing 130: describe(my_snpmatrix[:,2]) after zscore() and after the weights")
    describe(my_snpmatrix[:,2])

    println("===============================================================")
    println("========== CONSIDER DROPPING THE INTERCEPT ====================")
    println("========== TO KEEP THE MATRIX SPARSER =========================")
    println("===============================================================")
    println("Gordon testing 130,421: Set Intercept here or NOT, and in MendelIHT_utilities.jl")
    if USE_INTERCEPT # && false
        my_snpmatrix = [ones(size(my_snpmatrix, 1)) my_snpmatrix]
    end
    println()
    println("Gordon testing 130: size(my_snpmatrix) = $(size(my_snpmatrix)) after adding the column of ones for Intercept.")

    println("===============================================================")
    println("========== SETUP TEST DATA HERE ===============================")
    println("===============================================================")
    println(RED*"SETUP optional 'after zscore' TEST DATA HERE !!! "*DEFAULT)
    TEST_DATA = [3,5]
    if 1 in TEST_DATA           # Note: Designed for use when weights are OFF
                                #       this test data DID crash the program at 26 iterations with a dimensions error
                                #       I didn't follow up on it because it worked great before and now I'm working
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

    println()
    println("===============================================================")
    println("===============================================================")
    println("========== MAKE MANUAL ADJUSTMENTS IN SNPMATRIX HERE ==========")
    println("===============================================================")
    println(RED*"SETUP optional 'manual adjustments' TEST DATA HERE !!! "*DEFAULT)

    if 2 in TEST_DATA           # special - outlier DROP for 2√pq .\ 1
        println(RED*"LOOK FOR INTERESTING DATA HERE !!! "*DEFAULT)
        println(find(my_snpweights .> 3.0))

        println(find(my_snpweights .== maximum(my_snpweights)))  # 10.012523
        my_snpweights[7265] = 0
        my_snpmatrix[:,7265+1] = 0
        println(find(my_snpweights .== maximum(my_snpweights)))  #
        my_snpweights[6603] = 0
        my_snpmatrix[:,6603+1] = 0
        println(find(my_snpweights .== maximum(my_snpweights)))  #
        my_snpweights[6157] = 0
        my_snpmatrix[:,6157+1] = 0
        println(find(my_snpweights .== maximum(my_snpweights)))  #
        my_snpweights[4161] = 0
        my_snpmatrix[:,4161+1] = 0
        println(find(my_snpweights .== maximum(my_snpweights)))  #
        my_snpweights[4241] = 0
        my_snpmatrix[:,4241+1] = 0
        println(find(my_snpweights .== maximum(my_snpweights)))  #



        println(find(my_snpweights .== 1.0))    # [607, 4887, 9536, 9971]

        describe(my_snpmatrix_temp_after_zscore[:,7265])
        describe(my_snpmatrix[:,7265+1])
        describe(my_snpmatrix[:, 607+1])
        describe(my_snpmatrix[:,4887+1])
        describe(my_snpmatrix[:,9536+1])
        describe(my_snpmatrix[:,9971+1])
        describe(my_snpmatrix_temp_after_zscore[:,9971]) # max is 7 !!!!
        describe(my_snpmatrix_temp_after_zscore[:,9971+1]) # max is 7 !!!!
    end

    if 99 in DEBUG_LEVEL        # special MOVE APRIORI winners up front for easier visualization
        println()
        println("===============================================================")
        println("===============================================================")
        println("========== MOVE THE WINNERS UP FRONT ==========================")
        println("===============================================================")
        println(RED*"SETUP optional 'manual adjustments' TEST DATA HERE !!! "*DEFAULT)
    end
    if 99 in DEBUG_LEVEL        # special MOVE APRIORI winners up front for easier visualization
        my_snpmatrix[:,10] = 1  # fake intercept, works good
        my_snpmatrix[:,10] = 0  # cancel fake intercept
        my_snpmatrix[:,11] = 0
        my_snpmatrix[:,898] = 0
        my_snpmatrix[:,1880] = 0
        my_snpmatrix[:,3981] = 0
        my_snpmatrix[:,4398] = 0
        my_snpmatrix[:,5793] = 0
        my_snpmatrix[:,6611] = 0
        my_snpmatrix[:,6627] = 0
        my_snpmatrix[:,7023] = 0
        my_snpmatrix[:,9467] = 0
        describe(my_snpmatrix[:,1]) # with zcore and weight
        describe(my_snpmatrix[:,2])
        describe(my_snpmatrix[:,3])
        describe(my_snpmatrix[:,4])
        describe(my_snpmatrix[:,10])

        println(BLUE*"THE BOXPLOT FOR THE WINNERS ARE PLOTTED IN boxplot_my_winners.png"*DEFAULT)
        #StatPlots.boxplot(y=rand(1:4,1000),group=randn(1000))
        StatPlots.boxplot(rand(1:4,1000),randn(1000))
        println("Size of rand(1:4,1000) = $(size(rand(1:4,1000)))")
        println("Size of rand(1000) = $(size(rand(1000)))")
        # ALL PEOPLE, FIRST 50 SNPs
        plotdata = my_snpmatrix[:,1:50]
        println("Size of my_snpweights = $(size(my_snpweights))")
        println("Size of my_snpweights[1,:] = $(size(my_snpweights[1,:]))")
        println("Size of my_snpweights[:,1] = $(size(my_snpweights[:,1]))")
        StatPlots.boxplot(plotdata)
        Plots.savefig("boxplot_my_winners.png")
        # ALL SNPs, FIRST 50 PEOPLE
        #=
        plotdata = my_snpmatrix[1:50,:]
        println("Size of my_snpweights = $(size(my_snpweights))")
        println("Size of my_snpweights[1,:] = $(size(my_snpweights[1,:]))")
        println("Size of my_snpweights[:,1] = $(size(my_snpweights[:,1]))")
        StatPlots.boxplot(plotdata)
        Plots.savefig("boxplot_first50people.png")
        =#
        # ALL SNPs, FIRST 50 PEOPLE
        plotdata = transpose(my_snpmatrix)[:,1:2200]
        println("Size of my_snpweights = $(size(my_snpweights))")
        println("Size of my_snpweights[1,:] = $(size(my_snpweights[1,:]))")
        println("Size of my_snpweights[:,1] = $(size(my_snpweights[:,1]))")
        StatPlots.boxplot(plotdata)
        Plots.savefig("boxplot_first50people_transpose2200.png")


    end
    #= DUPLICATE CODE was used for offline testing of weight functions - OK TO DELETE
    fill!(v.xb, 0.0) #initialize β = 0 vector, so Xβ = 0
    println("Gordon testing 134: size(v.xb) = $(size(v.xb))")
    copy!(v.r, y)    #redisual = y-Xβ = y  CONSIDER BLASCOPY!
    println("Gordon: size(y) = $(size(y))")
    v.r[mask_n .== 0] .= 0 #bit masking? idk why we need this yet
    println("Gordon testing 137: size(v.r) = $(size(v.r))")

    # calculate the gradient v.df = -X'(y - Xβ) = X'(-1*(Y-Xb)). Future gradient
    # calculations are done in iht!. Note the negative sign will be cancelled afterwards
    # when we do b+ = P_k( b - μ∇f(b)) = P_k( b + μ(-∇f(b))) = P_k( b + μ*v.df)

    # Can we use v.xk instead of snpmatrix?
    print("Gordon testing 143: size(v.df), size(my_snpmatrix), size(v.r) = ")
    println(size(v.df), size(my_snpmatrix), size(v.r))
    At_mul_B!(v.df, my_snpmatrix, v.r)
    =#
    println()
    println("===============================================================")
    println("===============================================================")
    println("========== PRIOR WEIGHTS - END CALCULATIONS ===================")
    println("===============================================================")
    println("===============================================================")
    println("===============================================================")
    return my_snpMAF, my_snpweights, my_snpmatrix
end
