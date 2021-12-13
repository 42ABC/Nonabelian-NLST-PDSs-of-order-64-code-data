#Find PDSs with k=18 in grps of order 64 with an EA(8) image\

#Acknowledgements
#The procedure that iterates over normal subgroups to find the groups we are interested in is based on code written by Dr. Smith

#Run with
#Read("PDS_18_EA8img2.g");

#CHANGEPATH
Read("helper_functions.g");
Read("11IncidenceMatrices.txt");
Read("1ConvolutionTable.txt");

#Input: a group g
#Output: Precomputes group information for the accompanying Java file PDS18_looprun.java
givePDS18_EA8img := function(g)
    local ti, parmset, counter, failedcharmap, failedcharmap2, failedcharmap3, temp, PDS_list, PDS, \
    isPDS, sss, sss2, sss3, chrmap, pos, charmaps, charmap1, charmap2, charmap3, k1, k2, k3, k4, k5, k6, \
    k7, k8, a1, a2, a3, a4, a5, a6, a7, a8, actual_choices_by_level, possible_parmsets, \
    posschars, iii, tcharmap, gens, e, k, j, i, charcombo, MT, positions_notorder2elts_inverse, positions_notorder2elts, \
    order2elts, positions_order2elts, level, positions_elts_by_sg, positions_elts_by_sg_inverse, combos_list, \
    notorder2elts, c1, numparms, pos_parmsetsindex_temp, addit, numorder2elts, elt, elts_by_sg, lvl, rightcosets, \
    subgroup1, possible_sols, e_subgroup, ng, h8,failedconvcheck, listy, listy2, ppos, allsixes, i1, i2, i3, conprod, \
    efacgroup, facgroup, testmat, eltpos, gensfg,goon,conprodmid,ii1,ii2,ii3,ii4,intssofar, gid, hom, fg, efg, mat, matinv, \
    rowwith18, row;


    gid := IdGroup(g)[2];
    Print("Group: ", gid, "\n");

    e := Elements(g);
    gens := MinimalGeneratingSet(g);
    ng := NormalSubgroups(g);


    for k in ng do
        if IdGroup(FactorGroup(g,k)) = [8,5] then   # SmallGroup(8, 5) = C2 x C2 x C2. Change this line for other images
            subgroup1 := k;
            break;

        fi; 
    od; 

    e_subgroup := Elements(subgroup1);


    ###################################################################
    ###################################################################


    #STEP 1
    #GET PARM SETS


    ###################################################################
    ###################################################################


    #Getting the characters in the smaller group in the right order
    #Promised that the factor group will be like (C_2)^3
    #Note this doesn't need to be efficient because of how small it is
    hom := NaturalHomomorphismByNormalSubgroup(g, subgroup1);
    fg := Image(hom, g); #using a homomorphism command instead to get the factor group
    gensfg := MinimalGeneratingSet(fg); #assumes gensfg will have the three order 2 generators we want at the front, why we're using MinimalGeneratingSet
    
    #Print("gensfg: ", gensfg, "\n");
    #Print(Order(gensfg[1]), " ", Order(gensfg[2]), "\n");

    efg := Elements(fg); 
    #Print("efg: ", efg,"\n");

    posschars := dfslist(3,List([1..2],x -> (-1)^x)); #match up to generators
        
    #build the character matrix in the order of efg
    mat := [];
    for row in [1..Size(posschars)] do
        Add(mat, ListWithIdenticalEntries(8,-500)); #to make sure correctly loading in (we would know if a spot was missed)
        for i in [1..2] do
            for j in [1..2] do
                for k in [1..2] do
                    eltpos := Position(efg, gensfg[1]^i * gensfg[2]^j *gensfg[3]^k); #find the position in the list we want
                    if eltpos=fail then continue; fi; #skip if the inverse has already been included #redundant for c2^3

                    mat[row][eltpos] := posschars[row][1]^i * posschars[row][2]^j * posschars[row][3]^k;
       
                od;
            od;
        od;
    od;

    rowwith18 := -1;
    for i in [1..Size(mat)] do
        if ForAll(mat[i], x->x >= 0) then
            rowwith18 := i;
            break;
        fi;
    od;

    matinv := mat^-1;


    possible_sols := List(dfslist(7,[-6,2]), function(x) Add(x,18,rowwith18); return x; end);

    possible_parmsets := List(possible_sols, combo -> matinv * combo);


    possible_parmsets := Set(possible_parmsets);

    pos_parmsetsindex_temp := [];

    for i in [1..Size(possible_parmsets)] do
        if Minimum(possible_parmsets[i]) >= 0 then
            Add(pos_parmsetsindex_temp, i);
        fi;
        
    od;

    possible_parmsets := List(pos_parmsetsindex_temp, i -> possible_parmsets[i]);

    #
    ###########
    #Actual choices by level piece
    ###########
    #
    rightcosets := [];
    for elt in efg do #note that this came from a hom command! (NOT a FactorGroup call)
        Add(rightcosets, PreImagesElm(hom, elt));
    od;

    elts_by_sg := List(rightcosets, x -> Elements(x));
    positions_elts_by_sg := List(elts_by_sg, x -> List(x,y -> Position(e,y)));
    positions_elts_by_sg_inverse := List(elts_by_sg, x -> List(x, y -> Position(e,y^-1)));

    

    #Note that we need to add 1 to the number we seek to get to the correct index
    combos_list := List([0..8], x -> Combinations([1..8], x));


    #count the order 2 elements in each coset:
    order2elts := [];
    positions_order2elts := [];
    numorder2elts := 0*[1..8]; #by level
    #and the not order 2 elts
    notorder2elts := [];
    positions_notorder2elts := [];
    positions_notorder2elts_inverse := [];

    for lvl in [1..8] do
        Add(order2elts,[]);
        Add(positions_order2elts, []);
        Add(notorder2elts, []);
        Add(positions_notorder2elts, []);
        Add(positions_notorder2elts_inverse,[]);
        for i in [1..8] do
            elt := elts_by_sg[lvl][i];
            if Order(elt)=2 or Order(elt)=1 then #order 2 elts and the identity (anything that does not need to be paired)
                Add(order2elts[lvl], elt);
                Add(positions_order2elts[lvl], Position(e,elt));
            else
                Add(notorder2elts[lvl], elt);
                Add(positions_notorder2elts[lvl], positions_elts_by_sg[lvl][i]);
                Add(positions_notorder2elts_inverse[lvl], positions_elts_by_sg_inverse[lvl][i]);
            fi;
        od;
        numorder2elts[lvl] := Size(positions_order2elts[lvl]);
    od;

    #then whittle down parm sets (only examine those that are possible)
    pos_parmsetsindex_temp := []; #add the good ones here
    for i in [1..Size(possible_parmsets)] do
        addit := true;
        for j in [1..8]  do
            if RemInt(possible_parmsets[i][j],2)=1 and numorder2elts[j]=0 then #an odd elt is needed but no odd elt is present
                addit := false;
                break;
            fi;
        od;
        if addit then
            Add(pos_parmsetsindex_temp, i);
        fi;
    od;

    possible_parmsets := List(pos_parmsetsindex_temp, i -> possible_parmsets[i]);



    #first dimension is level
    #second dimension is the # of C_2 x C_2 x C_2 subgroup elts being chosen: 0,1,2,...,8. Assuming each time that this is a valid number
    #third dimension is a list of all the possible positions that could be added
    actual_choices_by_level := [];

    #For ALL levels
    #must pick the inverse as well, meaning that numbers here on out must all be even
    #there are order 2 elts as well though
    for level in [1..8] do
        Add(actual_choices_by_level,[]);


        #this is choosing 0 
        Add(actual_choices_by_level[level], [[]]); 

        #this is one coset elt chosen
        Add(actual_choices_by_level[level], []); 
       
        for pos in positions_order2elts[level] do
            Add(actual_choices_by_level[level][2], [pos]); #note adding by one to get the indexing right (this is choosing one parm but the index is 2)
        od;

        if Size(actual_choices_by_level[level][2]) = 0 then Add(actual_choices_by_level[level][2], []); fi;

        #picking two elts: two cases (two odd or elt & inverse)


        #pick two odd elts
        Add(actual_choices_by_level[level],[]);#this is choosing two coset elts (param=2)
        for i in [1..Size(order2elts[level])] do
            for j in [i+1..Size(order2elts[level])] do
                temp := [positions_order2elts[level][i], positions_order2elts[level][j]];
                if Size(temp) = Size(Set(temp)) and Size(temp)=2 then 
                    Add(actual_choices_by_level[level][3], Set(temp));
                fi;
            od;
        od;

        #or can pick an elt and its inverse
        for i in [1..Size(notorder2elts[level])] do
            temp := [positions_notorder2elts[level][i], positions_notorder2elts_inverse[level][i]];
            if Size(temp)=Size(Set(temp)) and Size(temp)=2 then
                Add(actual_choices_by_level[level][3], Set(temp));
            fi;
        od;

        #Error("pause\n");

        actual_choices_by_level[level][3] := Set(actual_choices_by_level[level][3]); #will take out redundant choices
                
                
        if Size(actual_choices_by_level[level][3]) = 0 then Print("surprise!\n"); Add(actual_choices_by_level[level][3], []); fi; #should never happen


        #3 parms
        #pick from parm 2 and parm 1
        Add(actual_choices_by_level[level], []); 
        for i in [1..Size(actual_choices_by_level[level][2])] do
            for j in [1..Size(actual_choices_by_level[level][3])] do
                temp := Concatenation(actual_choices_by_level[level][2][i], actual_choices_by_level[level][3][j]);
                if Size(temp) = Size(Set(temp)) and Size(temp)=3 then
                    Add(actual_choices_by_level[level][4], Set(temp));
                fi;
            od;
        od;
        actual_choices_by_level[level][4] := Set(actual_choices_by_level[level][4]);


        if Size(actual_choices_by_level[level][4]) = 0 then Add(actual_choices_by_level[level][4], []); fi;


        #4 parms
        Add(actual_choices_by_level[level], []); 
        for i in [1..Size(actual_choices_by_level[level][3])] do
            for j in [i+1..Size(actual_choices_by_level[level][3])] do
                temp := Concatenation(actual_choices_by_level[level][3][i], actual_choices_by_level[level][3][j]);
                if Size(temp) = Size(Set(temp)) and Size(temp)=4 then
                    Add(actual_choices_by_level[level][5], Set(temp));
                fi;
            od;
        od;
        actual_choices_by_level[level][5] := Set(actual_choices_by_level[level][5]);
        
        if Size(actual_choices_by_level[level][5]) = 0 then Print("surprise!\n"); Add(actual_choices_by_level[level][5], []); fi;


        #5 parms: 4 parm & 1 parm combo
        Add(actual_choices_by_level[level], []); 
        for i in [1..Size(actual_choices_by_level[level][2])] do
            for j in [1..Size(actual_choices_by_level[level][5])] do
                temp := Concatenation(actual_choices_by_level[level][2][i], actual_choices_by_level[level][5][j]);
                if Size(temp) = Size(Set(temp)) and Size(temp)=5 then
                    Add(actual_choices_by_level[level][6], Set(temp));
                fi;
            od;
        od;
        actual_choices_by_level[level][6] := Set(actual_choices_by_level[level][6]);

        if Size(actual_choices_by_level[level][6]) = 0 then Add(actual_choices_by_level[level][6], []); fi;


        #6 parms: same as 4 parm & 2 parm
        Add(actual_choices_by_level[level], []); 
        for i in [1..Size(actual_choices_by_level[level][3])] do
            for j in [1..Size(actual_choices_by_level[level][5])] do
                temp := Concatenation(actual_choices_by_level[level][3][i], actual_choices_by_level[level][5][j]);
                if Size(temp) = Size(Set(temp)) and Size(temp)=6 then
                    Add(actual_choices_by_level[level][7], Set(temp));
                fi;
            od;
        od;
        actual_choices_by_level[level][7] := Set(actual_choices_by_level[level][7]);
        if Size(actual_choices_by_level[level][7]) = 0 then Print("surprise!\n"); Add(actual_choices_by_level[level][7], []); fi;


        #7 parms: same as 6 parm & 1 parm
        Add(actual_choices_by_level[level], []); 
        for i in [1..Size(actual_choices_by_level[level][2])] do
            for j in [1..Size(actual_choices_by_level[level][7])] do
                temp := Concatenation(actual_choices_by_level[level][2][i], actual_choices_by_level[level][7][j]);
                if Size(temp) = Size(Set(temp))and Size(temp)=7 then
                    Add(actual_choices_by_level[level][8], Set(temp));
                fi;
            od;
        od;
        actual_choices_by_level[level][8] := Set(actual_choices_by_level[level][8]);

        if Size(actual_choices_by_level[level][8]) = 0 then Add(actual_choices_by_level[level][8], []); fi;


        #8 parms: well this means you're choosing everything
        #for fun let's say it's 6 parms and 2 parms and see if it works
        #less efficient but not that long anyway :) CBT
        Add(actual_choices_by_level[level], []); 
        for i in [1..Size(actual_choices_by_level[level][3])] do
            for j in [1..Size(actual_choices_by_level[level][7])] do
                temp := Concatenation(actual_choices_by_level[level][3][i], actual_choices_by_level[level][7][j]);
                if Size(temp) = Size(Set(temp)) and Size(temp)=8 then
                    Add(actual_choices_by_level[level][9], Set(temp));
                fi;
            od;
        od;
        actual_choices_by_level[level][9] := Set(actual_choices_by_level[level][9]);
        
        if Size(actual_choices_by_level[level][9]) = 0 then Print("big surprise!\n"); Add(actual_choices_by_level[level][9], []); fi;


    od;

    #
    #######
    #PDS testing framework
    #######
    #

    #CHANGEPATH all throughout here

    #Convolution table method
    MT := ConvolutionTable_f(g,e,x->x);
    allsixes := ListWithIdenticalEntries(64,6);

    #Print out the info to text files to use in Java
    PrintTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"),"//Table for grp ", gid, ":\n");

    #Print the possible parmsets:
    AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), Size(possible_parmsets), "\n"); 
    for i in [1..Size(possible_parmsets)] do
        for k in [1..8] do
            AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), possible_parmsets[i][k], " "); 
        od;
        AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), "\n"); #
    od;



    #print the convolution table first
    for i in [1..64] do
        for j in [1..64] do
            AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), MT[i][j], " "); #print the convolution table first
        od;
        AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), "\n"); #print the convolution table first
    od;

    #start printing actual_choices_by_level
    for i in [1..8] do #the level
        AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), i, "\n"); #denote the level
        for j in [1..9] do #the number of parameters chosen
            AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), j-1," ", Size(actual_choices_by_level[i][j]), " ", "\n"); #denote the level
            for choice in actual_choices_by_level[i][j] do #each choice gets its own line
                for k in choice do
                    AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), k, " "); 
                od;
                AppendTo(Concatenation("PDS18_gap_tables/gaptables64.",String(gid),".txt"), "\n"); 
            od;
        od;
    od;
    return -1;
        
    end;


grpswithc2c2c2img := [55..266]; #all of them

for id in grpswithc2c2c2img do
    Print(id,"\n");
    givePDS18_EA8img(SmallGroup(64,id));
od;
