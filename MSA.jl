import Pkg
Pkg.add("FASTX")
Pkg.add("Distributions")
Pkg.add("BioSequences")
Pkg.add("Plots")

using FASTX
using Distributions
using BioSequences
using Random
using Plots

function codify_chromosome(msa)
    chromosome = [[0 for j in 1:length(msa[1])] for i in 1:length(msa)] 
    for i in 1:length(msa)
        count = 0
        for j in 1:length(msa[i])
            if msa[i][j] == '-'
                chromosome[i][j] = count*(-1) 
            else
                count = count + 1
                chromosome[i][j] = count
            end
        end
    end
    
    return chromosome
end

function remove_columns_with_all_gaps(sequences, codifications)
    columns_with_all_gaps = []
    total_seq = length(sequences)
    
    for (index,nucleotide) in enumerate(sequences[1])
        if nucleotide == '-'
            column_with_all_gaps = true
            
            for i in 2:total_seq
                if sequences[i][index] != '-'
                    column_with_all_gaps = false
                    break
                end
            end
            
            if column_with_all_gaps
                append!(columns_with_all_gaps, index)
            end
        end
    end
    
    for i in 1:total_seq
        sequences[i] = join([c for (i,c) in enumerate(sequences[i]) if ~(i in columns_with_all_gaps)])
        codifications[i] = [c for (i,c) in enumerate(codifications[i]) if ~(i in columns_with_all_gaps)]
    end
                    
    return sequences, codifications
end
                                

function get_cut_point(cut_elements,chromosome_coded)
    cut_point = 0
    for i in 1:length(chromosome_coded)
        if abs(chromosome_coded[i])==abs(cut_elements[1])
            cut_point = i
            break
        end
    end
    return cut_point
end

function crossover(chromosome_string1, chromosome_coded1,chromosome_string2, chromosome_coded2)
    minimum = min(length(chromosome_string1[1]), length(chromosome_string2[1]))
    cut_point = rand(2:minimum-1)

    diff1,diff2 = Inf,-Inf
    
    cut_points_string2 = []
    
    for i in 1:length(chromosome_string1)
        cut_elements =  [chromosome_coded1[i][cut_point],chromosome_coded1[i][cut_point+1]]
        cut_point2 = get_cut_point(cut_elements,chromosome_coded2[i])
        diff1 = min(diff1,cut_point2)
        diff2 = max(diff2,cut_point2)
        push!(cut_points_string2,cut_point2)
    end
    
    crossover_string = [[],[]]
    crossover_coded = [[],[]]
    
    for i in 1:length(chromosome_string1)
        cut_point2 = cut_points_string2[i]
        s1 = chromosome_string1[i][1:cut_point]
        s2 = chromosome_string2[i][cut_point2+1:length(chromosome_string2[i])]
        
        gap_count1 = Int(cut_point2 - diff1)
        gaps1 = "-"^gap_count1
        
        main_string1 = string(s1,gaps1,s2)        
        main_coded1 = codify_chromosome([main_string1])[1]
        
        s1 = chromosome_string2[i][1:cut_point2]
        s2 = chromosome_string1[i][cut_point+1:length(chromosome_string1[i])]

        gap_count2 = Int(diff2 - cut_point2)
        gaps2 = "-"^gap_count2

        main_string2 = string(s1,gaps2,s2)        
        main_coded2 = codify_chromosome([main_string2])[1]
        
        push!(crossover_string[1], main_string1)
        push!(crossover_coded[1], main_coded1)
        
        push!(crossover_string[2], main_string2)
        push!(crossover_coded[2], main_coded2)
    end
    crossover_string[1],crossover_coded[1] = remove_columns_with_all_gaps(crossover_string[1],crossover_coded[1])
    crossover_string[2],crossover_coded[2] = remove_columns_with_all_gaps(crossover_string[2],crossover_coded[2])

    return crossover_string,crossover_coded
    
end

function mutation(parent_sequence, parent_codification)
    child_sequence = []
    child_codification = []
    
    for i in 1:length(parent_codification)
        curr_codification = parent_codification[i]
        curr_sequence = parent_sequence[i]
        
        gaps = []
        gap_start = -1
        for j in 1:length(curr_sequence)
            if curr_sequence[j] == '-'
                if gap_start == -1
                    gap_start = j
                end
            else
                if gap_start != -1
                    push!(gaps,[gap_start,j-1])
                    gap_start = -1
                end
            end
        end
        if gap_start != -1
            push!(gaps,[gap_start,length(curr_sequence)])
        end
        
        if length(gaps) == 0
            mutated_seq = join([c for (i,c) in enumerate(curr_sequence)])
            mutated_codification = [c for (i,c) in enumerate(curr_codification)]
            
            push!(child_sequence, mutated_seq)
            push!(child_codification, mutated_codification)
            
            continue
        end
        
        gap_group = gaps[rand(1:length(gaps))]
        
        mutated_seq = join([c for (i,c) in enumerate(curr_sequence) if ~(i>=gap_group[1] && i<=gap_group[2])])
        mutated_codification = [c for (i,c) in enumerate(curr_codification) if ~(i>=gap_group[1] && i<=gap_group[2])]
                        
        target_position = rand(1:length(mutated_seq)+1)
        
        gap_sequence = join(['-' for i in gap_group[1]:gap_group[2]])
        
        if target_position == length(mutated_seq)+1
            mutated_seq = mutated_seq * gap_sequence
            gap_code = -1*abs(mutated_codification[length(mutated_codification)])
            mutated_codification = cat(mutated_codification,[gap_code for i in gap_group[1]:gap_group[2]], dims=1)
        elseif target_position == 1
            mutated_seq = gap_sequence * mutated_seq
            gap_code = 0
            mutated_codification = cat([gap_code for i in gap_group[1]:gap_group[2]], mutated_codification, dims=1)
        else
            mutated_seq = mutated_seq[1:target_position-1] * gap_sequence * mutated_seq[target_position:length(mutated_seq)]
            gap_code = -1*abs(mutated_codification[target_position-1])
            mutated_codification = cat(mutated_codification[1:target_position-1],[gap_code for i in gap_group[1]:gap_group[2]],mutated_codification[target_position:length(mutated_codification)],dims=1)          
        end
                                    
        push!(child_sequence, mutated_seq)
        push!(child_codification, mutated_codification)
    end
    
    child_sequence,child_codification = remove_columns_with_all_gaps(child_sequence, child_codification)
end


function add_initial_gaps(gap_count,sequence)
    gap_count = floor(gap_count)
    while gap_count!=0
        blank_count = Int(rand(1:gap_count))
        position  = rand(1:length(sequence))
        sequence = string(sequence[1:position],"-"^blank_count,sequence[position+1:length(sequence)])
        gap_count=gap_count-blank_count
    end
    return sequence
end

function generate_single_chromosome(msa)
    
    new_msa = []
    for i in 1:length(msa)
        push!(new_msa,"")
    end
    
    max_sequence_length = length(msa[1])
    max_sequence_index = 1
    for (i,sequence) in enumerate(msa)
        if length(sequence)>max_sequence_length
            max_sequence_length = length(sequence)
            max_sequence_index = i
        end
    end
    
    # Update Max Length Sequence
    
    gap_count = rand(Uniform(0.2,0.4)) * max_sequence_length
    new_msa[max_sequence_index] = add_initial_gaps(gap_count,msa[max_sequence_index])
    
    # Update Rest of the Sequences
    
    for i in 1:length(msa)
        if i==max_sequence_index
            continue
        end
        gap_counts = length(new_msa[max_sequence_index]) - length(msa[i])
        new_msa[i] = add_initial_gaps(gap_counts,msa[i])
    end
    return new_msa
end

function generate_initial_population(msa,population_count)        
    without_crossover = 0.2*population_count
    population_string = []
    population_coded = []
    for i in 1:without_crossover
        chromosome = generate_single_chromosome(msa)
        coded_chromosome = codify_chromosome(chromosome)
        chromosome,coded_chromosome = remove_columns_with_all_gaps(chromosome,coded_chromosome)
        push!(population_coded,coded_chromosome)
        push!(population_string,chromosome)
    end
        
    with_crossover = (0.8*population_count)/2
    for i in 1:with_crossover
        random1 = Int(rand(1:without_crossover))
        random2 = Int(rand(1:without_crossover))
        crossover_string, crossover_coded = crossover(population_string[random1],population_coded[random1],population_string[random2],population_coded[random2])
        population_string = vcat(population_string,crossover_string)
        population_coded = vcat(population_coded,crossover_coded)
    end
    return population_string, population_coded
end    

function compute_tc_score(sequences)
    no_of_aligned_columns = 0
    total_columns = length(sequences[1])
    is_aligned = true
    for i in 1:total_columns
        is_aligned = true
        first_residue= sequences[1][i]
        if  first_residue != '-'
            for j in 1:length(sequences)
                if first_residue != sequences[j][i]
                    is_aligned = false
                    break
                end
            end
            else is_aligned = false
        end
        if is_aligned 
            no_of_aligned_columns += 1
        end
    end
    score = (100.0 * no_of_aligned_columns) / total_columns
    return score 
end

# +1 for match and 0 for no match

function calc_sum_pair(sequences) 
    t = length(sequences)
    k = length(sequences[1])
    score = 0
    for i=1:t
        A = sequences[i]
        for j=i+1:t
            B = sequences[j]
            for idx = 1:k
                if A[idx] == B[idx] && A[idx] != '-'
                    score += 1
                end
            end
        end
    end
    return score
end

function calc_nogap_percentage(sequences)
    t = length(sequences)
    k = length(sequences[1])
    no_gaps = 0
    for i in 1:t
        for j in 1:k
            if sequences[i][j]!='-'
                no_gaps+=1
            end
        end
    end
    total = k*t
    no_gaps_score = (no_gaps/total)*100
    return no_gaps_score
end
    
# +1 match. -1 mismatch and -2 for gaps
function calc_sum_pair2(sequences) 
    t = length(sequences)
    k = length(sequences[1])
    score = 0
    for i=1:t
        A = sequences[i]
        for j=i+1:t
            B = sequences[j]
            for idx = 1:k
                if A[idx]=='-' || B[idx]=='-'
                    score -=2
                elseif A[idx] == B[idx]
                    score += 1
                else
                    score -= 1
                end
            end
        end
    end
    return score
end

function calculate_fitness_score(align)
    return [compute_tc_score(align), calc_sum_pair(align), calc_nogap_percentage(align), calc_sum_pair2(align)]
end

# Basic structure to store the individuals and helper functions

mutable struct Individual{Align, Codification, Fitness}
    align::Align
    codification::Codification
    y::Fitness
    rank::UInt16
    crowding::Float64
    dom_count::UInt16
    dom_list::Vector{UInt16}
    
    Individual(align::Align, codification::Codification ,y::Fitness) where {Align,Codification,Fitness} = new{Align,Codification,Fitness}(align, codification, y,zero(UInt16),0.,zero(UInt16),UInt16[]) 
end

function dominates(a::Individual, b::Individual)
    result = false
    for i in eachindex(a.y)
        a.y[i] < b.y[i] && return false
        a.y[i] > b.y[i] && (result = true)
    end
    return result
end

function isless(a::Individual, b::Individual)
    return a.rank < b.rank || a.rank == b.rank && a.crowding >= b.crowding
end


function calculate_Pareto_Fronts!(population)
    n = length(population)
    for p in population
        empty!(p.dom_list)
        p.dom_count = 0
        p.rank = 0
    end
    
    ranked_individuals = 0
    
    current_rank_index = []
    for i in 1:n
        for j in i+1:n
            if dominates(population[i], population[j])
                push!(population[i].dom_list, j)
                population[j].dom_count += 1
            elseif dominates(population[j], population[i])
                push!(population[j].dom_list, i)
                population[i].dom_count += 1
            end
        end
        if population[i].dom_count == 0
            population[i].rank = 1
            ranked_individuals += 1
            append!(current_rank_index, i)
        end
    end
    
    current_rank = UInt16(2)
    while ranked_individuals != n
        indiv_to_delete = deepcopy(current_rank_index)
        empty!(current_rank_index)
        
        for i in indiv_to_delete
            p = population[i]

            for q in p.dom_list
                population[q].dom_count -= one(UInt16)
                if population[q].dom_count == zero(UInt16)
                    population[q].rank = current_rank
                    ranked_individuals += 1
                    append!(current_rank_index, q)
                end
            end

        end
        current_rank += one(UInt16)        
    end
end

function tournament_selection(P)
    a,b = rand(1:length(P)÷2), rand(1:length(P)÷2)
    if isless(P[a], P[b])
        return P[a]
    else
        return P[b]
    end
end

function crowding_distance_sorting!(population)
    for p in population
        p.crowding = 0.
    end
    for j = 1:length(first(population).y)
        sort!(population, by = x -> x.y[j])
        population[1].crowding = population[end].crowding = Inf
        if population[1].y[j] != population[end].y[j]
            for i in 2:length(population)-1
                population[i].crowding += (population[i+1].y[j] - population[i-1].y[j]) / (population[end].y[j] - population[1].y[j])
            end
        end
    end
end

function NSGA(input_seq, population_size, num_of_gen, pmut)
    alignments, codes = generate_initial_population(input_seq, population_size)
    
    P = Vector{Individual}(undef, 2*population_size)    # Whole population 
    
    for i in 1:population_size                                # Create Initial population
        P[i] = Individual(alignments[i],codes[i], calculate_fitness_score(alignments[i])) 
        P[population_size+i] = deepcopy(P[i])
    end
    
    calculate_Pareto_Fronts!(view(P, 1:population_size))
    
    for gen in 1:num_of_gen
        for i = 1:2:population_size
            
            pa = tournament_selection(P)
            pb = tournament_selection(P)
            
            childs_align, childs_codification = crossover(pa.align, pa.codification, pb.align, pb.codification)
            
            if rand() < pmut
                childs_align[1],childs_codification[1] = mutation(childs_align[1],childs_codification[1])
            end
                
            if rand() < pmut
                childs_align[2], childs_codification[2] = mutation(childs_align[2],childs_codification[2])
            end
            
            P[population_size+i] = Individual(childs_align[1],childs_codification[1], calculate_fitness_score(childs_align[1])) 
            P[population_size+i+1] = Individual(childs_align[2],childs_codification[2], calculate_fitness_score(childs_align[2])) 
            
        end
        calculate_Pareto_Fronts!(P)
        
        sort!(P, by = x -> x.rank, alg = Base.Sort.QuickSort)
        
        let f::Int = 1
            ind = 0
            indnext = findlast(x -> x.rank == f, P)
            while 0 < indnext <= population_size
                ind = indnext
                f += 1
                indnext = findlast(x -> x.rank == f, P)
            end
            indnext == 0 && (indnext = length(P))
            crowding_distance_sorting!(view(P, ind+1:indnext))
            sort!(view(P, (ind + 1):indnext), by = x -> x.crowding, rev = true, alg = PartialQuickSort(population_size - ind))
            
        end
        
    end
    optimal_sols = filter(x -> x.rank == 1, view(P, 1:population_size))
    
    best_sol = optimal_sols[1]
    for i in 2:length(optimal_sols)
        if optimal_sols[i].y[2] > best_sol.y[2]
             best_sol = optimal_sols[i]
        end 
    end
    
    return best_sol
    
end


# Read sequences from the input file
# Reference: https://biojulia.net/FASTX.jl/dev/manual/fasta/

function read_input(filename)
    reader = open(FASTA.Reader, filename)
    
    identifiers = []
    sequences = []
    
    for record in reader
        push!(identifiers, string(FASTA.identifier(record)))
        push!(sequences, string(FASTA.sequence(record)))
    end
    close(reader)
    
    return identifiers, sequences
end
#read_input("input1.txt")


function generate_alignments(file_name)
    identifiers, sequences = read_input(file_name)
    population_size = 200
    num_of_generations = 500
    mutation_prob = 0.8
    
    println("Input sequences:")
    for each in sequences
        println(each)
    end
    println()
    
    println("Alignment process started:")
    println()
    
    optimal_individual = @timed NSGA(sequences, population_size, num_of_generations, mutation_prob)
    optimal_align = optimal_individual.value.align
    score = optimal_individual.value.y[2]
    
    
    println("Generated alignments:")
    for each in optimal_align
        println(each)
    end
    println()
    println("Sum of pairs scor[+1 for match, 0 otherwise]: ", score) 
    println()
    println("Total time taken to align the sequence (in seconds): ", optimal_individual.time)
    println()
    
    open("output1.txt", "w") do file
        for i in 1:length(optimal_align)
            write(file, ">"*identifiers[i]*"\n")
            write(file, optimal_align[i]*"\n")
        end
    end
    
    println("This alignment is stored in the output1.txt file")
    
end
#generate_alignments("input1.txt")

