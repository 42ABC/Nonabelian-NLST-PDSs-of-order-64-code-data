#Check PDSs found to see if they are equivalent

#Input: The lists of PDSs found in the PDS18_output folder
#Output: A list of lists. Each element of the list is of the form [gid, [x1, x2, ..., xn]]. gid is a group number, and x1, x2, ..., xn are the inequivalent PDSs in the group SmallGroup(64,gid)

#run code with
#Read("checkPDSInequivalence3.g");

#CHANGEPATH
Read("helper_functions.g");
Read("11IncidenceMatrices.txt");
Read("1ConvolutionTable.txt");

LoadPackage("design"); #just to make sure its there (probably have already loaded it in one of the reads)


all_ineq_pds := [];

for i in [55..266] do
    Read(Concatenation("PDS18_output/pds18_64.",String(i),".txt")); #CHANGEPATH
    if Size(pds_list)=0 then
        continue;
    fi;

    g := SmallGroup(64,i);
    Print("group: ", i, "\n");

    e := Elements(g);
    aut := AutomorphismGroup(g);
    eaut := Elements(aut);
    #Print("Size of aut: ", Size(aut),"\n");

    #element form   
    PDS_list_eform := [];

    for pd in pds_list do
        Add(PDS_list_eform, List(Set(pd[2]), x-> e[x]));
    od;

    #Print("# of PDSs to check: ", Size(PDS_list_eform), "\n");

    ineq_pds := [];

    eq_pds := [1];
    pds_not_in := 1;
    while pds_not_in >= 0 do
        #Print("Round of auto!\n");
        Print(pds_not_in,"\n");
        Add(ineq_pds,pds_not_in);
        for k in [1..Size(eaut)] do
            pos := Position(PDS_list_eform, ImagesSet(eaut[k],PDS_list_eform[pds_not_in]));
            if pos=fail then 
                continue; 
            else
                Add(eq_pds, pos);
            fi;
        od;
        eq_pds := Set(eq_pds);
        pds_not_in := -1;
        for k in [1..Size(PDS_list_eform)] do
            if not (k in eq_pds) then
                pds_not_in := k;
                break;
            fi;
        od;
    od;

    Add(all_ineq_pds,[i,ineq_pds]);

    #Print("By auts: ", Size(Set(eq_pds)),"\n");
od;


PrintTo("PDS18Ineq_table.txt","ineq_PDS_table := ", all_ineq_pds, ";\n"); #CHANGEPATH
