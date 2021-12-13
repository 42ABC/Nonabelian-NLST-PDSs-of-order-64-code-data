#Determine the possible coset distributions of PDSs in groups of order 64 with k=18
#Read("PDS18_cosetdistribution.g");
LoadPackage("difsets");

#manual dset check
isdsetish := function(g,dset,checkpos)
    local e,col,num2,i,j,num6;
    e := Elements(g);
    col := ListWithIdenticalEntries(64,0);
    for i in dset do
        for j in dset do
            
            col[Position(e,i * j^-1)] := col[Position(e,i * j^-1)] + 1;
            
        od;
    od;
    num2 := 0;
    num6 := 0;
    for i in checkpos do
        if col[i]=2 then
            num2 := num2 +1;
        elif col[i]=6 then
            num6 := num6 + 1;
        fi;
    od;
    #Print(checkpos,"\n");

    #Print(Collected(col),"\n");
    #Print(num2,"\n");
    if col[1]=6 and num2=15 then
        return true;
    else 
        return false;
    fi;
    end;


revable := [];
Read("PDS18Ineq_table.txt"); #CHANGEPATH
for i in [1..Size(ineq_PDS_table)] do
    rev_i := [];
    Read(Concatenation("PDS18_output/pds18_64.",String(ineq_PDS_table[i][1]),".txt")); #CHANGEPATH
    g := SmallGroup(64,ineq_PDS_table[i][1]);
    e := Elements(g);
    for j in [1..Size(ineq_PDS_table[i][2])] do

        pds := pds_list[ineq_PDS_table[i][2][j]][2];

        ng := NormalSubgroups(g);
        for k in ng do
            if IdGroup(FactorGroup(g,k))=[4,2] then
                #Print("yay\n");
                sg := k;       
                e_subgroup1 := Elements(sg);
                #Print("normal?", IsNormal(g,subgroup1),"\n");

                hom := NaturalHomomorphismByNormalSubgroup(g, sg);
                fg := Image(hom, g); #FactorGroup(g,subgroup1); #using a homomorphism command instead to get the factor group
                efg := Elements(fg);
                pdsimg := List(pds,x->Image(hom,e[x]));
                gensfg := MinimalGeneratingSet(fg); #assumes gensfg will have the three order 2 generators we want at the front, why we're using MinimalGeneratingSet

                #Print(pdsimg,"\n");
                if Size(Collected(pdsimg))=3 then  
                    #Error("pause"); 

                    rightcosets := List(efg,x->PreImagesElm(hom,x));
                    reps := List(rightcosets,x->Representative(x));
                    #Print("reps: ", reps,"\n");

                    dset1 := [];
                    dset2 := [];
                    dset3 := [];
                    checkpos := [];
            
                    for pos in [1..64] do
                        if Image(hom,e[pos])=efg[1] then
                            Add(checkpos,pos);
                    
                        fi;
                    od;

                    for pos in pds do 
                        if Image(hom,e[pos])=efg[2] then
                            Add(dset1,e[pos]); fi;
                        if Image(hom,e[pos])=efg[3] then
                            Add(dset2,e[pos]); fi;
                        if Image(hom,e[pos])=efg[4] then
                            Add(dset3,e[pos]);  fi;
                    od;
                    #Print(dset1, " ", dset2, " ",dset3, "\n");

                    if isdsetish(g,dset1,checkpos) and isdsetish(g,dset2,checkpos) and isdsetish(g,dset3,checkpos) then
                        #Print("Three dsets!\n");

                        inv_dset1 := List(dset1, x -> x^-1);
                        int1 := Size(Intersection(dset1, inv_dset1));
        
                        inv_dset2 := List(dset2, x -> x^-1);
                        int2 := Size(Intersection(dset2, inv_dset2));
        
                        inv_dset3 := List(dset3, x -> x^-1);
                        int3 := Size(Intersection(dset3, inv_dset3));
        
                        if int1=6 and int2=6 and int3=6 then #note that the 6 is specialized for order 16
                            Print("3rev: ", ineq_PDS_table[i][1], " ",ineq_PDS_table[i][2][j],"\n");
                            Add(rev_i,ineq_PDS_table[i][2][j]);
                            #Print(dset1,"\n");
                            #Print(dset2,"\n");
                            #Print(dset3,"\n");
                            #Error("pause");

                            break;
                        fi;
                    fi;

                fi;
            fi;
        od;
    od;

    Add(revable,[ineq_PDS_table[i][1],rev_i]);

od;

Print("printing revable:\n");
for i in [1..Size(revable)] do
    Print(revable[i],"\n");
od;

PrintTo("PDS18_revablePDSs.txt","revable := ", revable, ";\n"); #CHANGEPATH