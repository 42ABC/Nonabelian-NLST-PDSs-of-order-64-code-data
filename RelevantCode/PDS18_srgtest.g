#Goal: Determine how many of the PDSs we found have the same strongly regular graph
#list of PDSs and corresponding groups
#Index 1: group ID
#Index 2: the PDS
#Index 3: the graph
#Index 4: the isomorphism ID (which equivalence class it belongs to)

#Read("PDS18_srgtest.g");

#CHANGEPATH
LoadPackage("grape");
Read("PDS18Ineq_table.txt");
Read("11IncidenceMatrices.txt");

#create the graphs
PDS_graphs := [];
counter := 1;
for i in [1..Size(ineq_PDS_table)] do
    Read(Concatenation("PDS18_output/pds18_64.",String(ineq_PDS_table[i][1]),".txt")); #CHANGEPATH
    g := SmallGroup(64,ineq_PDS_table[i][1]);
    for j in [1..Size(ineq_PDS_table[i][2])] do

        pds := pds_list[ineq_PDS_table[i][2][j]][2];
        
        adj := IncidenceMatrix(LeftTranslates_Nums(g,pds)); #adjacency matrix
        Add(PDS_graphs,[ineq_PDS_table[i][1],ineq_PDS_table[i][2][j],Graph( Group(()), [1..64], OnPoints, function(x,y) return adj[x][y]=1; end, true),counter, -1]); #from GRAPE
        counter := counter + 1;
    od;
od;

# Print("before comparisons: \n");
# for i in [1..Size(PDS_graphs)] do
#     Print(PDS_graphs[i][1]," " , PDS_graphs[i][4],"\n");
# od;


#compare the PDSs to each other
for i in [1..Size(PDS_graphs)] do
    for j in [i+1..Size(PDS_graphs)] do
        if IsIsomorphicGraph(PDS_graphs[i][3],PDS_graphs[j][3]) then
            PDS_graphs[j][4] := PDS_graphs[i][4]; #set j's iso ID to i
        fi;
    od;
od;



PDS_graphs := SortedList(PDS_graphs, function(x,y) return x[4] < y[4]; end); #sort by the id

Print("after comparisons: \n");
for i in [1..Size(PDS_graphs)] do
    Print(PDS_graphs[i][1]," " , PDS_graphs[i][2], " " , PDS_graphs[i][4],"\n");
od;

#Store the graphs in a smart way
graph_iso_classes := List([1..10],x->[x,[]]); #note that the graphs line up with Spence's first 10 graphs

#key: iso class
#value : dictionary with key of grp, value of list of PDSs in that grp with that iso class
iso_dict := NewDictionary(0,true);
for i in [1..10] do
    AddDictionary(iso_dict,i,NewDictionary(0,true));
    for j in [1..49] do
        AddDictionary(LookupDictionary(iso_dict,i),ineq_PDS_table[j][1],[]);
    od;
od;

#create Spence's graphs
Read("Spence_srg64_cleaned.txt"); #CHANGEPATH
spence_srg_graphform := [];
for i in [1..Size(spence_srg)] do
    Add(spence_srg_graphform,Graph( Group(()), [1..64], OnPoints, function(x,y) return spence_srg[i][x][y]=1; end, true));
od;
for index in [1..Size(PDS_graphs)] do
    for j in [1..Size(spence_srg)] do
        if (IsIsomorphicGraph(PDS_graphs[index][3],spence_srg_graphform[j])) then
            Print(PDS_graphs[index][1], " is ismorphic to graph# ", j, "\n");
            PDS_graphs[index][5] := j;
            #Add(graph_iso_classes[j],[PDS_graphs[i][1],PDS_graphs[i][2]]);
            Add(LookupDictionary(LookupDictionary(iso_dict,j),PDS_graphs[index][1]),PDS_graphs[index][2]);

            break; #can't be isomorphic to more than one graph

        fi;
    od;
od;

for i in [1..10] do
    for j in [1..49] do
        if Size(LookupDictionary(LookupDictionary(iso_dict,i),ineq_PDS_table[j][1])) > 0 then
            Add(graph_iso_classes[i][2],[ineq_PDS_table[j][1],LookupDictionary(LookupDictionary(iso_dict,i),ineq_PDS_table[j][1])]);
        fi;
    od;
od;

PrintTo("PDS18_fullsrgWithSpence_graphiso.txt","graph_iso_classes := ", graph_iso_classes, ";\n"); #CHANGEPATH