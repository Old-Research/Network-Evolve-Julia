using DifferentialEquations
using NLsolve
using Random
using DelimitedFiles
using Dates

#--------STRUCTURES-------------------------
mutable struct Reaction
    substrate::Array{String,1}
    product::Array{String,1}
    reg::Array{String,1}
    rxnType
    enzyme::Float64
end

mutable struct Model
    boundary::Array{String,1}
    floating::Array{String,1}
    rxns::Array{Reaction,1}
    kVals::Array{Float64,1}
    qVals::Array{Float64,1}
    by::Array{Float64,1}
    #fitness::Float64
end
#--------READ/WRITE FUNCTIONS-------------------------

function makeLogFile(name::String)
    # Input: a name (string) for the log file
    # Output: path to log file with timestamp (string)
    timeStamp = now()
    dt = Dates.format(timeStamp,"yyyy-mm-dd_HH.MM.SS" )
    fileName = string(pwd(),"\\",name,"_",dt,".txt")
    return fileName
end

function writeHeader(f)
    # Record start time, genetic algorithm info, and notes to file f
    open(f, "a") do io
        writedlm(io, [string("START TIME: ", Dates.format(now(),"mm-dd HH:MM")),
        string("Individuals per Generation: ", nIndividuals), string("Elite Selection: ", nElite), string("Randomly Generated: ", nRand), notes])
    end
end

function writeGen(f,j,generation)
    # Record fitness and k values of generation j in file f
    open(f, "a") do io
        writedlm(io, [string("GENERATION ", j)])
        writedlm(io, generation)
    end
end

function writeFooter(f)
    # Record stop time and list of fitness values to file f
    open(f, "a") do io
        writedlm(io, [string("END TIME: ", Dates.format(now(),"HH:MM")),"K parameters varied: ", positions, "Fitness Over Time: ", fitnessOverTime])
    end
end

#--------BASIC FUNCTIONS----------------------
function getNumberOfFloatingSpecies(m::Model)
    return length(m.floating)
end

function getNumberOfReactions(m::Model)
    return length(m.rxns)
end

function simulate(m, u0, tspan)
    # Produces time series data for each species, useful for plotting
    prob = ODEProblem(odefunc!, u0, tspan, m)
    sol = solve(prob)
    return sol
end

function getSteadyState(du, fy, m,)
    t=0
    SS = nlsolve((du,fy) ->odefunc!(du,fy,m,t), u0)
    return SS.zero
end

function odefunc!(du,fy, m, t)
    #fy = list of floating species initial concentration
    by = m.by #by = list of boundary species initial concentration
    nReactions = length(m.rxns)
    nFloatingSpecies = length(m.floating)
    floatingSpecies = m.floating
    boundarySpecies = m.boundary
    reactions = m.rxns
    for i = 1:nFloatingSpecies
      du[i] = 0.0
    end
    #du = zeros(nFloatingSpecies)
    dv = zeros(nReactions) # reaction rates

    for (rindex, reaction) in enumerate(reactions)

        sp = m.rxns[rindex].substrate # ALL subtrates (strings)
        pp = m.rxns[rindex].product # ALL products (strings)

        #Get substrate (sp) concentration
        s = zeros(length(sp)) # Empty substrate concentration array
        for i = 1:length(m.rxns[rindex].substrate) #for every substrate in the reaction
            bIndex = findall(boundarySpecies -> boundarySpecies==sp[i], boundarySpecies)
            #first check if it's a boundary species
            if bIndex != []
                s[i] = by[bIndex[1]] #if it's a boundary species, assign it to whatever is in boundary species list at that index
            else #If it's not a boundary species, find it in the floating species list
                sIndex = findall(floatingSpecies -> floatingSpecies==sp[i],floatingSpecies)
                if sIndex != []
                    s[i] = fy[sIndex[1]]
                else
                println("Substrate species not found:", sp[i])
                end
            end
        end

        #Get product (pp) concentration
        p = zeros(length(pp))
        for i = 1:length(m.rxns[rindex].product) # for every product in the reaction
            bIndex = findall(boundarySpecies -> boundarySpecies==pp[i], boundarySpecies)
            #If it is a boundary species, get it from the boundary species list
            if bIndex != []
                p[i] = by[bIndex[1]]
            else #If it's not a boundary species, look in floating species for product
                sIndex = findall(floatingSpecies -> floatingSpecies==pp[i],floatingSpecies)
                if sIndex != []
                    p[i] = fy[sIndex[1]]
                else
                println("Product species not found: ", pp[i])

                end
            end
        end


        if m.rxns[rindex].rxnType == uniuni
            dv[rindex]=m.rxns[rindex].rxnType(s[1], p[1], m.kVals[rindex][1], m.kVals[rindex][2])
        elseif m.rxns[rindex].rxnType == biuni
            dv[rindex]=m.rxns[rindex].rxnType(s[1], s[2], p[1], m.kVals[rindex][1], m.kVals[rindex][2])
        elseif m.rxns[rindex].rxnType == unibi
            dv[rindex]=m.rxns[rindex].rxnType(s[1], p[1], p[2], m.kVals[rindex][1], m.kVals[rindex][2])
        elseif m.rxns[rindex].rxnType == bibi
            dv[rindex]=m.rxns[rindex].rxnType(s[1], s[2], p[1], p[2], m.kVals[rindex][1], m.kVals[rindex][2])
        elseif m.rxns[rindex].rxnType == basicEnzyme
            dv[rindex]=m.rxns[rindex].rxnType(s[1], p[1], m.rxns[rindex].enzyme, m.kVals[rindex], m.qVals[rindex])
        end

    end

    for (sindex, fsp) in enumerate(floatingSpecies)
        for (rindex, reaction) in enumerate(reactions)
            if fsp in reaction.substrate
                du[sindex] = du[sindex] - dv[rindex]
            end
            if fsp in reaction.product
                du[sindex] = du[sindex] + dv[rindex]
            end
        end
    end
    return du
end
#--------REACTION TYPES----------------------
function basicEnzyme(substrate::Float64, product::Float64, enzyme::Float64,
    k::Float64, q::Float64)
    return enzyme*k*(substrate - (product/q))
end

function uniuni(substrate::Float64,product::Float64,k1::Float64,k2::Float64)
    return k1*substrate-product*k2
end

function biuni(substrate1::Float64, substrate2::Float64, product::Float64,
    k1::Float64, k2::Float64)
    return k1*(substrate1+substrate2) - k2*product
end

function unibi(substrate::Float64, product1::Float64, product2::Float64,
    k1::Float64, k2::Float64)
    k1*substrate - k2*(product1+product2)
end

function bibi(substrate1::Float64, substrate2::Float64,product1::Float64,
    product2::Float64, k1::Float64, k2::Float64)
    return k1*(substrate1+substrate2) - k2*(product1+product2)
end

#--------"GENETIC ALGORITHM" FUNCTIONS----------

function randomK!(m::Model; positions=collect(1:length(m.rxns)))
    # Input model
        # Default: randomize ALL k values
        # kawgs: positions can be set to a custom integer array to randomize specific k values
        # or to randomize a certain number of k values in any position

    # Output model with random k in ith position
    if positions==collect(1:length(m.rxns))
        m.kVals = rand(0.01:0.01:100, (length(m.rxns)))
    else
        for i in positions
            m.kVals[i] = rand(0.01:0.01:100)
        end
    end
    return m
end

function insertRandModel(generation, start_index, stop_index; nK=Int64[])
    # input: Array with fitness in first column and k values in second column
        # If the optional argument nK is set to an integer, then that number of
        # k values will be randomly selected and varied
        # and the rest of the k values from the previous generation will remain
    # output: the SAME array with new models with random k values and their
    # fitness scores in rows start_index through stop_index
    # Places random models and their fitness scores into the generation array
    for i = start_index:stop_index
        randM = createModel()
        if nK != Int64[] # If the optional argument is used
            randM.kVals = generation[i,2] # carry over previous k values
            for j = 1:nK[1]
                # vary a single random k value, repeat nK times,
                # (could end up being the same k value more than once)
                randM = randomK!(randM, positions =
                rand(1:1:length(randM.rxns))) # randomly vary a single k value
            end
        end
        fitScore = fitness(truePerturbMat, randM)
        if fitScore > generation[i,1]
            β = rand(1)
            if β[1] > 0.8 # accept 20%
                generation[i,1] = fitScore
                generation[i,2] = randM.kVals
            end
            break # reject new model
        else # If new fitness is better than old
            generation[i,1] = fitScore
            generation[i,2] = randM.kVals
        end
        randM = []
    end
    return generation
end

function insertMutatedModel(generation, start_index)
    # input: Array with fitness in first column and k values in second column
    # output: the same array with random models and their fitness scores in
    # rows start_index through start_index + nElite
    # Places random models and their fitness scores into the generation array
    # with 80/20 acceptance/rejection if the new model has worse fitness
    n = 1
    for i = start_index:start_index+nElite
        newM = createModel(bsp, fsp, rxns, Q, b0, generation[n,2])
        newM = mutateK!(newM,[all])
        fitScore = fitness(truePerturbMat, newM)
        if fitScore > generation[n,1]
            β = rand(1)
            if β[1] > 0.8 # accept 20%
                generation[i,1] = fitScore
                generation[i,2] = newM.kVals
            end
            break # reject new model
        else # If new fitness is better than old
            generation[i,1] = fitScore
            generation[i,2] = newM.kVals
        end
        newM = []
        n += 1
    end
    return generation
end

function replaceK!(m,k)
    # Takes an array of k values and places them into model m
    # Returns new model m
    m.kVals = k
    return m
end

function mutateK!(m; positions = collect(1:length(m.rxns)), improve = true )
    # IN: model to be mutated by +/- 60%
        # Default: all k values are mutated, new model is (most likely ) more
        #fit than old
        # If positions is set to a custom integer array, then only k's in those
        # positions will vary
        # If improve = false, then the new model accepted regardless of fitness
    # Return model m and its fitness
    oldM = deepcopy(m)
    oldK = oldM.kVals
    check = false
    while check == false
        if positions == collect(1:length(m.rxns))
            newK = m.kVals
            up_down = rand(0.4:1.2:1.6, length(m.rxns))
            newK = newK.*up_down
            bool = newK.>0.01
            if (0 in bool) == true
                loc = findall(bool->bool==0,bool)
                for j in loc
                    newK[loc[j]] = 0.01
                end
            end
            m.kVals = newK
        else
            for i in positions
                m.kVals[i] = m.kVals[i]*rand(0.4:1.2:1.6)
                if m.kVals[i] < 0.01
                    m.kVals[i] = 0.01
                end
            end
        end

        if improve == true
            if fitness(truePerturbMat,oldM) > fitness(truePerturbMat,newM)
                check == false
            else
                check == true
            end
        else
            check == true
        end
    end
end

function crossover(generation, nElite)
    # Takes a sorted generation and randomly pairs the elite chromosomes
    # Outputs the modified UNSORTED generation array with the offspring k values
    # Pairing is random
    # nElite should be even, I'll put a check for this later unless we use a
    # different pairing strategy
    # 2 parents make 2 offspring, all 4 included in subsequent generation
    # Subsequent generation is the input array with modified values, no new
    # generation array created
    matingPop = randperm(nElite::Int64)
    nPair = Int64(nElite/2) # number of pairs
    nK = getNumberOfReactions(trueM) # number of k values
    mom = matingPop[1:nPair] # Array of row numbers of "mother" parameters
    dad = matingPop[nPair+1:end] # Array of row numbers of "father" parameters
    place = 1
    for i = 1:nPair
        crossPt1 = rand(1:1:nK)
        crossPt2 = rand(1:1:nK)
        while crossPt1 == crossPt2
            # prevent cross points from being the same, which can cause
            # duplication
            crossPt2 = rand(1:1:nK)
        end
        crossPt1, crossPt2 = min(crossPt1,crossPt2), max(crossPt1,crossPt2)
        kNew1 = deepcopy(generation[mom[i],2]) # deepcopy mother k values
        kNew2 = deepcopy(generation[dad[i],2]) # deepcopy father k values
        β = rand(1) # random number between 0 and 1
        for j = crossPt1:crossPt2
            kNew1[j] = β[1]*generation[mom[i],2][j] + (1-β[1])*
            generation[dad[i],2][j]
            kNew2[j] = (1-β[1])*generation[mom[i],2][j] + β[1]*
            generation[dad[i],2][j]
        end
        m1 = createModel(bsp, fsp, rxns, Q, b0, kNew1)
        m2 = createModel(bsp, fsp, rxns, Q, b0, kNew2)
        generation[nElite+place,2] = kNew1
        generation[nElite+place+1,2] = kNew2
        generation[nElite+place,1] = fitness(truePerturbMat,m1)
        generation[nElite+place+1,1] = fitness(truePerturbMat,m2)
        place += 2
    end
    return generation
end

function createPerturbationMatrix(m)
    # Create a perturbation matrix for the model m
    # Assumes that fy and u0 are already defined.
    # fy is the steady state concentration of floating species
    # u0 is the initial concentration of floating species
    simPerturbMat = Array{Any}(undef, getNumberOfFloatingSpecies(m),
    getNumberOfReactions(m))
    for i=1:getNumberOfReactions(m)
        # Get reference enzyme concentration
        refsol = nlsolve((du,fy) ->odefunc!(du,fy,m,t), u0)
        fy = refsol.zero # Steady state species concentrations
        E_ref = m.rxns[i].enzyme
        E_up = 1.05*E_ref
        E_down = 0.95*E_ref
        # Upregulate enzyme by 5%, get new steady state
        m.rxns[i].enzyme = E_up
        upsol = nlsolve((du,fy) ->odefunc!(du,fy,m,t), u0)
        S_up = upsol.zero
        # Downregulate enzyme by 5%, get new steady state
        m.rxns[i].enzyme = E_down
        downsol = nlsolve((du,fy) ->odefunc!(du,fy,m,t), u0)
        S_down = downsol.zero
        # Put scaled sensitivity values in Simulated Perturbation Matrix
        for j = 1:getNumberOfFloatingSpecies(m)
            simPerturbMat[j,i] = ((S_up[j]-S_down[j])/.1)*(E_ref/fy[j])
        end
        # Reset enzyme
        m.rxns[i].enzyme = E_ref
    end
    return simPerturbMat
end

function fitness(perturbationMatrix, m)
    # Input: true perturbation matrix and a model, m
    # Output: fitness score for model m
    # Quantifies the difference between an experimental and simulated
    # perturbation matrix
    # Matrices must be the same size
    simPerturbationMatrix = createPerturbationMatrix(m)
    return fitness  = sqrt(sum((simPerturbationMatrix.-perturbationMatrix).^2))
end

function createEmptyRxnArray(num)
    # In: how many reactions you want in the array
    # Out: array with num empty reactions
    # For now: 1 substrate and 1 product per reaction, all basicEnzyme type
    reactions = Array{Reaction,1}(undef,num)
    for i = 1:num
        reactions[i] = Reaction(["x"],["x"],["x"],basicEnzyme,1.)
    end
    return reactions
end

function createModel(;k=Float64[]::Array{Float64,1}, rxns=true::Bool,
    Q=trueQ::Array{Float64,1})::Model
    # This function creates a new instance of the type Model with random k vals
    # Input: Array of q values, array of k values, boolean
        # If createModel() is called with no arguments, output model will have
        # random k values,
        # reaction topology and q values identical to true model
    # Output: new model
        # k = [] (default) -- generate random k values
        # k = Float64 Array -- use input k values in model
        # rxns = true (default) -- use true reaction topology in new model
        # rxns = false -- output model will have empty reaction array
    if rxns == true
        if k == []
            kvals = Array{Float64,1}(undef, length(trueM.rxns))
            m = Model(bsp, fsp, trueM.rxns, kvals, Q, by)
            m = randomK!(m)
        else
        m = Model(bsp, fsp, trueM.rxns, k, Q, by)
        end
    end
    if rxns == false
        rxns = createEmptyRxnArray(length(trueM.rxns))
        if  k == []
            kvals = Array{Float64,1}(undef, length(trueM.rxns))
            m = Model(bsp, fsp, rxns, kvals, Q, by)
            m = randomK!(m)
        else
            m = Model(bsp, fsp, rxns, k, Q, by)
        end
    end
    return m
end

function createGen1(m, nIndividuals)
    # This function requires that the arguments for createModel, bsp, fsp, rxns,
    # Q, b0,
    # are defined elsewhere in the script
    # Input: inital model m, number of individuals in generation
    # Outputs a sorted generation of nIndividuals random models (random k
    # values)
    generation1 = Array{Any}(undef, nIndividuals, 2)
    for i = 1:nIndividuals
        temp = createModel(bsp, fsp, rxns, Q, b0)
        generation1[i,2] = temp.kVals
        generation1[i,1] = fitness(truePerturbMat, temp)
        temp = []
    end
    generation1 = generation1[sortperm(generation1[:,1]),:]
    return generation1
end

function sortByFitness(generation)
    # Input: an unsorted generation array with fitness in first column
    # Output: the same array sorted by fitness (best to worst)
    generation = generation[sortperm(generation[:,1]),:]
    return generation
end

function createNewGenElite(generation1,nIndividuals,nElite)
    # Create array for a new generation with nIndividuals
    # Carry over nElite from generation1
    # Return a NEW array for a new generation that is empty save for elites
    newGeneration = Array{Any}(undef, nIndividuals, 2)
    newGeneration[1:nElite,:] = generation1[1:nElite,:]
    return newGeneration
end

function createNetworkModel(bsp,fsp,k=trueK,q=trueQ,enzyme=enzyme)
    # ATTN: 2020-06-22 there's a bug that allows Xo to be selected as a product
    # Input: string arrays for boundary and floating species
        # Optional arguments: k values, q values, enzymes.
        # Set to true values by default for now.
    # Output: model with random toplogy
    # Puts in true k and q values, basicEnzyme as reaction type, and ignores
    # reg for now.

    reactions = fill(Reaction(["s"], ["p"], ["reg"], basicEnzyme, 1.0),
    (length(enzyme)))
    sub_list = Array{String,1}(undef,0)
    enzyme_list = shuffle!(deepcopy(enzyme))
    sp_sets = fill(("s","p"), (length(enzyme)))
    reg = ["reg"]

    # First substrate: Xo
    # First product: random selection from floating species
    s1 = bsp[1]
    α = rand(1:1:length(fsp))
    p1 = fsp[α]
    sp_sets[1] = (s1,p1)
    reactions[1].substrate = [s1]
    reactions[1].product = [p1]

    # Second substrate: random selection from floating species
    # Second product: X1
    α = rand(1:1:length(fsp))
    s2 = fsp[α]
    p2 = bsp[2]
    sp_sets[2] = (s2,p2)
    reactions[2].substrate = [s2]
    reactions[2].product = [p2]

    # Substrate selection for reactions 3:end
    for i = 3:length(enzyme)
        # substrate selection
        if (i-2) <= length(fsp)
            s = fsp[i-2]
        else # once every floating species is used, pick randomly
            β = rand(1)
            if β[1] < 0.05
                s = bsp[1]
            else
                α = rand(1:1:length(fsp))
                s = fsp[α]
            end
        end
        push!(sub_list,s)
    end

    # product selection for reactions 3:end
    @label product_selection
    sp_sets[3:end] = fill(("s","p"),(length(sp_sets)-2,1))
    prod_list = deepcopy(fsp)
    for i = 3:length(rxns)
        s = sub_list[i-2]
        p = s
        while p == s
            set = (s, p)
            sp_sets[i] = set
            let set = set
                if length(prod_list) != 0
                    # If there are still species that have not been products,
                    # use those
                    if prod_list == [s]
                        # If the only product left to choose from is the same as
                        # the substrate, re-do selection
                        @goto product_selection
                    end
                    let counter = 1 # To prevent infinite loops
                        while isPairIn(set, sp_sets) == true
                            if counter > 10
                                @goto product_selection
                            end
                            α = rand(1:1:length(prod_list))
                            p = prod_list[α]
                            set = (s,p)
                            counter += 1
                        end
                    end #of let counter...
                else
                    #if every species has been a product once, choose randomly
                    # from all of them for remaining reactions
                    while isPairIn(set, sp_sets) == true
                        β = rand(1)
                        if β[1] < 0.05
                                p = bsp[2]
                        else
                            println(length(fsp))
                            α = rand(1:1:length(fsp))
                            p = fsp[α]
                        end
                        set = (s,p)
                    end #of while isPairIn while-loop
                end # if prod_list list != 0o
            end # of "let" block
        end #of while p==s loop
        filter!(prod_list->prod_list≠p,prod_list)
        sp_sets[i] = (s,p)
        reactions[i] = Reaction([s],[p],reg,basicEnzyme,enzyme_list[i])
    end # of reaction for-loop
    m = Model(bsp,fsp,reactions,trueK,trueQ,by)
    return m
end

function isPairIn(sp, spArray)
    # Checks if a substrate-product species pair is in another array without
    # respect to the order of the tuple
    # Meaning ("S2","S4") is equivalent to ("S4","S2") for example
    # IN: substrate-product tuple and an array of substrate-product tuples
    # OUT: boolean - true/false sp is in spArray
    check = false
    if (sp in spArray) == true
        check = true
    elseif ((sp[2],sp[1]) in spArray) == true
        check = true
    end
    return check
end

function isSetUnique(product,enzyme)
    # For use in network crossover function
    # Returns true if neither the product nor the enzyme has been used
    # previously
    if (product in pList) == false && (enzyme in eList) == false
        return true
    else
        return false
    end
end

function isUnique(item, list)
    # Checks if a substrate/product/enzyme has been used before.
    # In: a specific substrate, enzyme, or product and the list to check
        # sList, pList, or eList
    # Out: boolean, true if the item has NOT been previsouly used.
    check = item in list
    if check == false # if the item is NOT on the list, it IS unique
        return true
    else
        return false
    end
end


function testTopologyModel(m::Model)
    # This is a function to test if generated models have "valid" networks
    substrates = Array{String}(undef, length(m.rxns))
    products = Array{String}(undef, length(m.rxns))
    for i = 1:length(m.rxns)
        substrates[i] = m.rxns[i].substrate[1]
        products[i] = m.rxns[i].product[1]
    end
    # Make sure every floating species is used at least once as both a
    # product and reactant
    for sp in m.floating
        if (sp in substrates) == false
            error("Species ",sp," not found amongst substrates")
        elseif (sp in products) == false
            result = "fail"
            error("Species ",sp," not found amongst products")
        else
            result = "PASS"
        end
    end
    # Make sure the boundary species are properly connected (eg. Xo is
    # never a product)
    if (m.boundary[1] in products) == true
        error(m.boundary[1], "found amongst products")
        result = "FAIL"
    elseif (m.boundary[2] in substrates) == true
        error(m.boundary[2], "found amongst substrates")
        results = "FAIL"
    end
    return result
end

#--------MODEL---------------------------------


s1 = ["Xo"]; p1 = ["S1"]; reg1 = []; q1 = 1.5;
s2 = ["S1"]; p2 = ["S2"]; reg2 = []; q2 = 3.9;
s3 = ["S1"]; p3 = ["S3"]; reg3 = []; q3 = 1.1;
s4 = ["S2"]; p4 = ["S4"]; reg4 = []; q4 = 0.7;
s5 = ["S4"]; p5 = ["S5"]; reg5 = []; q5 = 2.3;
s6 = ["S3"]; p6 = ["S5"]; reg6 = []; q6 = 3.1;
s7 = ["S5"]; p7 = ["X1"]; reg7 = []; q7 = 4.0;


enzyme = [1; .5; 3; 1.7; 4.2; .7; 3.8]
trueQ = [1.5; 3.9; 1.1; 0.7; 2.3; 3.1; 4.]

trueK = [1.; 0.5; 3.;2.;1.;1.5;2.5]

r1 = Reaction(s1,p1,reg1,basicEnzyme,enzyme[1])
r2 = Reaction(s2,p2,reg2,basicEnzyme,enzyme[2])
r3 = Reaction(s3,p3,reg3,basicEnzyme,enzyme[3])
r4 = Reaction(s4,p4,reg4,basicEnzyme,enzyme[4])
r5 = Reaction(s5,p5,reg5,basicEnzyme,enzyme[5])
r6 = Reaction(s6,p6,reg6,basicEnzyme,enzyme[6])
r7 = Reaction(s7,p7,reg7,basicEnzyme,enzyme[7])

trueM = Model(["Xo", "X1"], ["S1", "S2","S3","S4","S5"], [r1, r2, r3, r4,
r5, r6, r7], trueK, trueQ, [10.,0.])

#--------ODE INFO------------------------
const u0 = zeros(length(trueM.floating))
fy = zeros(length(trueM.floating))
const t=0
du = zeros(length(trueM.floating))
const tspan = (0.,40.)
#--------TRUE PERTURBATION MATRIX
truePerturbMat = createPerturbationMatrix(trueM)
#--------MODEL CREATION INFO------------------------
rxns = trueM.rxns
fsp = trueM.floating
bsp = trueM.boundary
by = trueM.by



#-----------DEVELOPING TOPOLOGY CROSSOVER-----------------
#
# # Parents
# modelA = createNetworkModel(bsp,fsp)
# modelB = createNetworkModel(bsp,fsp)
#
# # Empty offspring
# newM = createModel(rxns = false) #has trueK and trueQ and true enzymes
# newM.rxns[1].substrate = [bsp[1]]
# newM.rxns[2].product = [bsp[2]]
#
# # Make arrays of (substrate, product) tuples for each model
# setA = Array{Tuple{String,String}}(undef, length(modelA.rxns),1)
# setB = Array{Tuple{String,String}}(undef, length(modelA.rxns),1)
# for i=1:length(modelA.rxns)
#     setA[i] = (modelA.rxns[i].substrate[1], modelA.rxns[i].product[1])
#     setB[i] = (modelB.rxns[i].substrate[1], modelB.rxns[i].product[1])
# end
#
# # Keep track of what substrates, products, enzymes have already been used
# tempfsp = deepcopy(fsp)
# sList = Array{String,1}(undef,0)
# pList = Array{String,1}(undef,0)
# eList = Array{Float64,1}(undef,0)
# SPset = fill(("s","p"), (length(modelA.rxns), 1)) # Keep track of offspring
# # s-p pairs
#
#
# # Make sure that boundary species are in offspring model:
# let check = false
#     while check == false
#         β = rand(1)
#         if β[1] > 0.5 # Randomly select product, enzyme from one parent for Xo substrate rxn
#             newM.rxns[1].product = modelA.rxns[1].product
#             newM.rxns[1].enzyme = modelA.rxns[1].enzyme
#         else
#             newM.rxns[1].product = modelB.rxns[1].product
#             newM.rxns[1].enzyme = modelB.rxns[1].enzyme
#         end
#
#         # Make sure that boundary product is in offspring
#         β = rand(1)
#         if β[1] > 0.5 # Randomly select substrate, enzyme from one parent for X1 product rxn
#             newM.rxns[2].substrate = modelA.rxns[2].substrate
#             newM.rxns[2].enzyme = modelA.rxns[2].enzyme
#             if newM.rxns[2].substrate == newM.rxns[1].product[1] # If selected substrate 2 = product 1, then chose the other substrate for 2
#                 newM.rxns[2].substrate = modelB.rxns[2].substrate
#                 newM.rxns[2].enzyme = modelB.rxns[2].enzyme
#             end
#         else
#             newM.rxns[2].substrate = modelB.rxns[2].substrate
#             newM.rxns[2].enzyme = modelB.rxns[2].enzyme
#             if newM.rxns[2].substrate == newM.rxns[1].product[1]
#                 newM.rxns[2].substrate = modelA.rxns[2].substrate
#                 newM.rxns[2].enzyme = modelA.rxns[2].enzyme
#             end
#         end
#         check = newM.rxns[1].enzyme != newM.rxns[2].enzyme # Make sure rxn 1 and 2 use different enzymes, otherwise re-do
#     end
# end
#
# # Put reaction 1 & 2 info in offspring model and store used substrates, products, and eznymes in lists
# newM.rxns[1].substrate = [bsp[1]]
# push!(pList, newM.rxns[1].product[1])
# push!(eList, newM.rxns[1].enzyme)
# SPset[1] = (newM.rxns[1].substrate[1], newM.rxns[1].product[1])
# newM.rxns[2].product = [bsp[2]]
# push!(sList, newM.rxns[2].substrate[1])
# push!(eList, newM.rxns[2].enzyme)
# SPset[2] = (newM.rxns[2].substrate[1], newM.rxns[2].product[1])
#
# # Check if any reactions are in both parents, pass down any shared reactions
# rxnNum = 3# Index of next offspring reaction to be filled
# for i = 3:length(modelA.rxns)
#     global rxnNum
#     if (setA[i] in setB) == true
#         #If the s-p pair is in both sets, pass it down to offspring, add substrate and product to tracking lists
#         newM.rxns[rxnNum].substrate = [setA[i][1]]
#         push!(sList, setA[i][1])
#         newM.rxns[rxnNum].product = [setA[i][2]]
#         push!(pList, setA[i][2])
#         SPset[rxnNum] = (newM.rxns[rxnNum].substrate[1], newM.rxns[rxnNum].product[1])
#         # Decide which enzyme to pass on
#         eA = modelA.rxns[i].enzyme
#         eIndex = findall(setB -> setB==setA[i], setB)
#         eB = modelB.rxns[eIndex[1][1]].enzyme
#         if (eA in eList) == false && (eB in eList) == false # If neither enzyme has been used, choose randomly
#             β = rand(1)
#             if β[1] > 0.5
#                 newM.rxns[rxnNum].enzyme = eA
#             else
#                 newM.rxns[rxnNum].enzyme = eB
#             end
#         elseif (eA in eList) == true && (eB in eList) == false # If the enzyme from A has been used, but not B, use B
#             newM.rxns[rxnNum].enzyme = eB
#         elseif (eA in eList) == false && (eB in eList) == true # If the enzyme from B has been used, but not A, use A
#             newM.rxns[rxnNum].enzyme = eA
#         else # If both enzymes have been previously used, model fails
#             println("FAIL: enzyme redundancy in reaction ",rxnNum)
#             break
#         end
#         push!(eList, newM.rxns[rxnNum].enzyme)
#         global rxnNum += 1
#     end
# end
#
# # Reactions that were NOT shared between parents
# for (fIndex, f) in enumerate(tempfsp) # For each floating species in the system
#     global rxnNum
#     global pA
#     global pB
#     global eA
#     global eB
#
#     if (f in sList) == false # If the species did not appear in the list of substrates that are alrady in the offspring model
#         newM.rxns[rxnNum].substrate = [f]
#         push!(sList, f)
#         # Find candidate products & enzymes: those affiliated with the species in question in each parent
#         for i = 1:length(modelA.rxns) # Looking through tuples lists to find the substrate
#             if f == setA[i][1]
#                 pA = setA[i][2] # p1 is the product affiliated with substrate f in setA
#                 eA = modelA.rxns[i].enzyme
#             end
#             if f == setB[i][1] # Now lookg through setB for f
#                 pB = setB[i][2]
#                 eB = modelB.rxns[i].enzyme
#             end
#         end
#
#         if isUnique(eA, eList) == false && isUnique(eB, eList) == false
#             # If both enzymes have been previously used in another reaction, no valid model
#             error("Enzyme redundancy in reaction ", rxnNum)
#         end
#
#         if (isSetUnique(pA,eA) == true && isSetUnique(pB,eB) == true) ||
#             (isUnique(eA,eList) == true && isUnique(eB,eList) == true &&
#             isUnique(pA,pList) == false && isUnique(pB,pList == false))
#             # Both enzymes and both products are unique
#             # OR both enzymes but neither products are unique
#             # -> randomly pick a set
#             β = rand(1)
#             if β[1] > 0.5
#                 newM.rxns[rxnNum].enzyme = eB
#                 newM.rxns[rxnNum].product = [pB]
#             else
#                 newM.rxns[rxnNum].enzyme = eA
#                 newM.rxns[rxnNum].product = [pA]
#             end
#         elseif isSetUnique(pA,eA) == true && isSetUnique(pB,eB) == false # Complete unique set for A but not B -> choose set A
#             newM.rxns[rxnNum].enzyme = eA
#             newM.rxns[rxnNum].product = [pA]
#         elseif isSetUnique(pA,eA) == false && isSetUnique(pB,eB) == true  # Complete unique set for B but not A -> choose set B
#             newM.rxns[rxnNum].enzyme = eB
#             newM.rxns[rxnNum].product = [pB]
#         elseif (isUnique(pA,pList) == true && isUnique(eA,eList) == false &&
#             isUnique(pB,pList) == false && isUnique(eB,eList) == true) # Only product A and enzyme B are unqiue -> choose them
#             newM.rxns[rxnNum].enzyme = eB
#             newM.rxns[rxnNum].product = [pA]
#         elseif (isUnique(pA,pList) == false && isUnique(eA,eList) == true &&
#             isUnique(pB,pList) == true && isUnique(eB,eList) == false) # Only product B and enzyme A are unique -> choose them
#             newM.rxns[rxnNum].enzyme = eA
#             newM.rxns[rxnNum].product = [pB]
#         elseif (isUnique(pA,pList) == false && isUnique(eA,eList) == true &&
#             isUnique(pB,pList) == false && isUnique(eB,eList) == false) # Only enzyme A is unique -> choose set A
#             newM.rxns[rxnNum].enzyme = eA
#             newM.rxns[rxnNum].product = [pA]
#         elseif (isUnique(pA,pList) == false && isUnique(eA,eList) == false &&
#             isUnique(pB,pList) == true && isUnique(eB,eList) == false) # Only enzyme B is unique -> choose set B
#             newM.rxns[rxnNum].enzyme = eB
#             newM.rxns[rxnNum].product = [pB]
#         else
#             error("Conditions not met for reaction ", rxnNum)
#         end
#     end
# end
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # Let's say we have an empty reaction spot after all this
#
# # Check which enzyme(s) are unaccounted for. Make an array of them
# # unusedE = Array{Float64, 1}(undef, 0)
# # for e in enzyme
# #     if (e in eList) == false
# #         push!(unusedE, e)
# #     end
# # end
# #
# #
# # for e in unusedE
# #     global rxnNum
# #     global modelA
# #     global modelB
# #     global eIndexA
# #     global eIndexB
# #     for i = 1:length(modelA.rxns)
# #         if e == modelA.rxns[i].enzyme
# #             eIndexA = i
# #         end
# #         if e == modelB.rxns[i].enzyme
# #             eIndexB = i
# #         end
# #     end
# #
# #     A_pair = setA[eIndexA] # S-P tuples affiliated with the missing enzyme in model A and B
# #     B_pair = setB[eIndexB]
# #
# #     # Check to make sure A_pair or B_pair are not making loops in new model
# #     if pairIsIn(A_pair, SPset) == true && pairIsIn(B_pair, SPset) == true # If both reactions form loop, reject model
# #         println("New Model Failed")
# #         break
# #     elseif pairIsIn(A_pair, SPset) == true && pairIsIn(B_pair, SPset) == false
# #         newM.rxns[rxnNum].substrate = [B_pair[1]]
# #         newM.rxns[rxnNum].product = [B_pair[2]]
# #     elseif pairIsIn(A_pair, SPset) == false && pairIsIn(B_pair, SPset) == true
# #         newM.rxns[rxnNum].substrate = [A_pair[1]]
# #         newM.rxns[rxnNum].product = [A_pair[2]]
# #     else # If neither reaction forms a loop, randomly pick one
# #         β = rand(1)
# #         if β[1] < 0.5
# #             newM.rxns[rxnNum].substrate = [A_pair[1]]
# #             newM.rxns[rxnNum].product = [A_pair[2]]
# #         else
# #             newM.rxns[rxnNum].substrate = [B_pair[1]]
# #             newM.rxns[rxnNum].product = [B_pair[2]]
# #         end
# #     end
# #     newM.rxns[rxnNum].enzyme = e
# #     push!(eList, e)
# #     push!(sList, newM.rxns[rxnNum].substrate[1])
# #     push!(pList, newM.rxns[rxnNum].product[1])
# #     rxnNum +=1
# # end
# #
