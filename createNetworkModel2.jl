function createNetworkModel(bsp,fsp,k=trueK,q=trueQ,enzyme=enzyme)
    # ATTN: 2020-06-22 there's a bug that allows Xo to be selected as a product
    # Input: string arrays for boundary and floating species
        # Optional arguments: k values, q values, enzymes.
        # Set to true values by default for now.
    # Output: model with random toplogy
    # Puts in true k and q values, basicEnzyme as reaction type, and ignores reg for now.

    # Let's allocate some arrays!
    reactions = fill(Reaction(["s"], ["p"], ["reg"], basicEnzyme, 1.0),(length(enzyme)))
    prod_list = deepcopy(fsp)
    sub_list = Array{String,1}(undef,0)
    enzyme_list = deepcopy(enzyme)
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
        if (i-2) <= length(fsp) # add each floating species as a substrate first
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
    println("we visited the product_selection label")
    sp_sets[3:end] = fill(("s","p"),(length(sp_sets)-2,1))
    for i = 3:length(rxns)
        s = sub_list[i-2] # Set s to first item on sub_list, in this case S1 (randomly picked as sub for X1)
        p = s
        while p == s # repeat product selection until a product that is different from the substrate is selected
            set = (s, p) # (S1, S1)
            sp_sets[i] = set
            let set = set
                if length(prod_list) != 0 # If there are still species that have not been products, use those
                    if prod_list == [s]
                        println("product list exception block reached") #If the only product left to choose from is the same as the substrate, back track 1 reaction and re-do selection
                        @goto product_selection
                        #@goto endfunction # Go back to beginning of product selection block
                    end
                    let counter = 1 # To prevent infinite loops
                        while isPairIn(set, sp_sets) == true #we got stuck in this loop, prod_list != [s]
                            if counter > 10
                                println("Invalid Model")
                                @goto product_selection
                                #@goto endfunction
                            end
                            α = rand(1:1:length(prod_list))
                            p = prod_list[α]
                            set = (s,p)
                            counter += 1
                        end
                    end #of let counter...
                else #if every species has been a product once, choose randomly from all of them for remaining reactions
                    while isPairIn(set, sp_sets) == true
                        println("we are in the while loop on line 84")
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
        α = rand(1:1:length(enzyme_list)) # randomly select an enzyme that has yet to be used
        e = enzyme_list[α]
        filter!(enzyme_list->enzyme_list≠e,enzyme_list)
        reactions[i] = Reaction([s],[p],reg,basicEnzyme,e)
    end # of reaction for-loop
    m = Model(bsp,fsp,reactions,trueK,trueQ,by)
    #return m
    @label endfunction
    return sp_sets
end
