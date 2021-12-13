#Read("NumberOfPDSs.g");

cc := 0;
for i in [55..266] do
    Read(Concatenation("PDS18_output/pds18_64.",String(i),".txt")); #CHANGEPATH
    cc := cc + Size(pds_list);
od;
Print(cc,"\n");