#Library of nice functions to help with Gap

#Use this to run file
#Read("helper_functions.g");

print_2D_list := function(list)
    local k,l;
    Print("##############################################\n");

    for l in list do
        for k in l do
            Print(k,"\n");
        od;
        Print("##############################################\n");
    od;
    end;


length_2D_list := function(list)
    return Sum(List(list, x -> Sum(x)));
    end;


#Inefficiently creates a list of all possible lists of length len, where each elm in the list is one of nums
#@param len the length of the lists of interest
#@param nums the numbers one is choosing from (with repetition)
#@return a list of lists with length len and entries all possible combos of numbers in nums
dfslist := function(len,nums)
    local lis;
    lis := [];
    dfslistrec(len,nums,1,[],lis);
    return lis;
    end;


#Better way to implement this?
dfslistrec := function(len,nums,ind,temp,lis)
    local i;

    #Print(ind, " ", ShallowCopy(temp), "\n");
    #Print(ind,"\n");

    if ind = len+1 then 
        Add(lis,List(temp,x->x)); #A very sad way to copy a list but it worked
    else 
        for i in nums do
            Add(temp,i);
            dfslistrec(len,nums,ind+1,temp,lis);
            Remove(temp); #removes the last elt
        od;
    fi;
    end;

#Returns if pds1 and pds2 are equivalent (by automorphisms)
isEquivalentPDS := function(g, pds1, pds2)
    local eq, aut, eaut, pds1_eform, pds2_eform, k, e;
    
    eq := false;
    e := Elements(g);
    aut := AutomorphismGroup(g);
    eaut := Elements(aut);

    #element form
    pds1_eform := List(Set(pds1), x -> e[x]);
    pds2_eform := List(Set(pds2), x -> e[x]);

    for k in [1..Size(eaut)] do
        if pds1_eform=ImagesSet(eaut[k], pds2_eform) then
            eq:=true;
            break;
        fi;
    od;

    return eq;
    end;

