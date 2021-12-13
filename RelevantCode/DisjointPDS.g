#Check PDSs found to see which groups contain a pair of disjoint PDSs
#Input: The PDSs found by exhaustive search
#Output: A list of lists. Each element of the list is of the form [gid, pds1, pds2]. 
#gid is a group number (representing the group SmallGroup(64,gid)
#pds1 and pds2 are PDSs in position form, corresponding to the list of group elements found by Elements(SmallGroup(64,gid))

#Acknowledgements:
#This code takes inspiration from Dr. Ken Smith's disjoint difference set program

#Run with
#Read("DisjointPDS.g");

#CHANGEPATH 
Read("helper_functions.g");
Read("11IncidenceMatrices.txt");
Read("1ConvolutionTable.txt");

LoadPackage("design"); #just to make sure its there (probably have already loaded it in one of the reads)
LoadPackage("difsets");

grps_with_dispds := [];

for gid in [55..266] do
    Read(Concatenation("PDS18_output/pds18_64.",String(gid),".txt")); #CHANGEPATH
    if Size(pds_list)=0 then
        continue;
    fi;
    g := SmallGroup(64,gid);
    Print("group: ", gid,"\n");
    e := Elements(g);
  
    found_one := false;

    for i in [1..Size(pds_list)] do
        for j in [i+1..Size(pds_list)] do
    
            #Print(Size(Intersection(pds1,pds2)),"\n");
            if Size(Intersection(pds_list[i][2],pds_list[j][2]))=0 then
                #Print("found two disjoint PDSs!\n");
                #Print(pds_list[i][2], " " , pds_list[j][2],"\n");
                Add(grps_with_dispds,[gid,pds_list[i][2],pds_list[j][2]]);

                found_one:=true;
                break;
                
            fi;
        od;
        if found_one then break; fi;
    od;

    

od;


Print(Size(grps_with_dispds),"\n");

PrintTo("PDS18_disjointPDSs_list.txt","grps_with_dispds :=", grps_with_dispds,";\n"); #CHANGEPATH