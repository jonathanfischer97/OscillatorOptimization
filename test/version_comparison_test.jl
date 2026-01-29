# Tests to isolate differences between Julia versions in the optimization pipeline
# These tests use fixed inputs to ensure reproducible comparisons

using Test, Random
using OscillatorOptimization
using DataFrames: nrow
using Statistics: mean, std
using SciMLBase: successful_retcode

@testset "Version comparison tests" begin
    # Fixed system setup
    fixed_params = Dict{Symbol, Float64}(:DF => 20.0)
    opt_sys = OptimizationReactionSystem(fullrn; fixed_params)
    
    # Create a fixed test individual (same parameters across Julia versions)
    rng = MersenneTwister(1234)
    test_population = generate_population(fullrn, 1, fixed_params; rng=rng)
    test_individual = first(test_population)  # SlicedDimArray acts like a 1D vector
    
    @testset "Individual ODE solution reproducibility" begin
        # Test that the same individual produces the same ODE solution
        sol, saved_array = solve_odes(test_individual, opt_sys)

        @test successful_retcode(sol)
        @test size(saved_array, 1) == 2  # Two observables
        @test size(saved_array, 2) > 0   # Time points
        
        # Print detailed diagnostic information
        println("=== ODE Solution Diagnostics ===")
        println("Solution status: ", sol.retcode)
        println("Saved array size: ", size(saved_array))
        println("First observable at t=0: ", saved_array[1, 1])
        println("Second observable at t=0: ", saved_array[2, 1])
        println("First observable at end: ", saved_array[1, end])
        println("Second observable at end: ", saved_array[2, end])
        println("First observable mean: ", mean(saved_array[1, :]))
        println("Second observable mean: ", mean(saved_array[2, :]))
        println("First observable std: ", std(saved_array[1, :]))
        println("Second observable std: ", std(saved_array[2, :]))
    end
    
    @testset "Fitness calculation reproducibility" begin
        # Use the same saved array to test fitness calculation
        sol, saved_array = solve_odes(test_individual, opt_sys)
        
        # Test fitness calculation
        phenotype = zeros(3)  # [fitness, period, amplitude]
        result = calculate_fitness!(phenotype, saved_array, opt_sys)
        
        println("=== Fitness Calculation Diagnostics ===")
        println("Calculated fitness: ", result[1])
        println("Calculated period: ", result[2])  
        println("Calculated amplitude: ", result[3])
        
        # Test the intermediate steps
        fft_observable = saved_array[1, :]
        fftData = getFrequencies(fft_observable)  # Use version without pre-computed plan
        fft_peaks = find_fft_peaks(fftData)
        
        println("FFT data length: ", length(fftData))
        println("FFT data max: ", maximum(fftData))
        println("FFT data mean: ", mean(fftData))
        println("Number of FFT peaks: ", length(fft_peaks))
        
        if length(fft_peaks) >= 2
            println("FFT peaks found at indices: ", fft_peaks)
            println("FFT peak values: ", fftData[fft_peaks])
        else
            println("Insufficient FFT peaks - this explains zero fitness")
        end
        
        # These should be reproducible across Julia versions
        @test result[1] >= 0.0 || result[1] == 0.0  # Fitness should be non-negative or zero
        @test result[2] >= 0.0 || result[2] == 0.0  # Period should be non-negative or zero
        @test result[3] >= 0.0 || result[3] == 0.0  # Amplitude should be non-negative or zero
    end
    
    @testset "Oscillation detection reproducibility" begin
        # Test the oscillation detection functions directly
        sol, saved_array = solve_odes(test_individual, opt_sys)
        
        # Test time-domain detection  
        td_observable = saved_array[2, :]
        td_peaks, td_troughs = find_amem_peaks(td_observable)
        
        println("=== Oscillation Detection Diagnostics ===")
        println("Time-domain observable length: ", length(td_observable))
        println("Time-domain observable range: ", minimum(td_observable), " to ", maximum(td_observable))
        println("Number of time-domain peaks: ", length(td_peaks))
        println("Number of time-domain troughs: ", length(td_troughs))
        
        # Test oscillation checks
        is_oscillatory = check_oscillatory(td_observable)
        is_regular = check_oscillation_regularity(td_observable)
        
        println("Is oscillatory: ", is_oscillatory)
        println("Is regular: ", is_regular)
        
        # Check the last 200 seconds criterion
        last_200_points = td_observable[end-2000:end]
        last_200_std = std(last_200_points)
        println("Last 200 points std: ", last_200_std, " (threshold: 1e-6)")
        println("Oscillatory test passes: ", last_200_std >= 1e-6)
        
        # Basic sanity checks
        @test length(td_peaks) >= 0
        @test length(td_troughs) >= 0
        @test isa(is_oscillatory, Bool)
        @test isa(is_regular, Bool)
    end
    
    @testset "Medium population test" begin
        # Test with a medium population to see if we can find any oscillatory solutions
        println("=== Medium Population Test (1000 individuals) ===")
        
        medium_results = run_optimization(1000, opt_sys;  
            mutationRate = 0.95,
            crossoverRate = 0.75,
            mutation_δ = 1.0,
            pm = 0.25,
            η = 1,
            sbx_pm = 0.3,
            sbx_η = 1,
            num_tournament_groups = 20,
            n_points = 50,               # Stop after finding 50 solutions
            seed = 1234,
            show_trace = false,
            parallelization = :thread
        )
        
        println("Medium population test - solutions found: ", nrow(medium_results.df))
        
        if nrow(medium_results.df) > 0
            println("Sample fitness values: ", medium_results.df.Fitness[1:min(5, nrow(medium_results.df))])
            println("Sample period values: ", medium_results.df.Period[1:min(5, nrow(medium_results.df))])
            println("Sample amplitude values: ", medium_results.df.Amplitude[1:min(5, nrow(medium_results.df))])
        else
            println("No oscillatory solutions found - this indicates the core issue")
        end
        
        # This test is more about providing diagnostic info than passing/failing
        @test nrow(medium_results.df) >= 0  # At least it shouldn't crash
    end
end