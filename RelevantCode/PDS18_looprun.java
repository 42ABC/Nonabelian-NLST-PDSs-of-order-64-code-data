//Run the main loop in a PDS 18 search code in Java, for possible speed boost

//Acknowledgements:

//Thanks to Nikita Morozov & Sreya Aluri for advice on how to make this code more efficient: \
//Specifically, using index assignment to avoid a remove call and not using a recursive method (iterative or for loops instead)
//Thanks to Dr. Prateek Bhakta for suggesting checks in the middle of the nested loops to reduce runtime, this helps a lot!
//The convolution check for PDSs was adapted from Drs. Smith and Kaliszewski's convolution table procedures

import java.util.*;
import java.io.*;

public class PDS18_looprun {

    static class PDS {
        
        
        int gid;
        int[] pds;
        boolean full;
        public PDS(int id,int[] temp) {
            this.gid=id;

            pds = new int[18];
            for (int i = 0; i < 18; i++) {
                pds[i]=temp[i];
            }
            full=true;
        }

        public PDS(int id) {
            pds = new int[18];
            full=false;

        }

        public String toString() {
            if (full) {
                return "[" + gid + ", " + Arrays.toString(pds) + "]";
            }
            else {
                return gid + " []";
            }
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {


        System.out.println(java.util.Calendar.getInstance().getTime());

        PrintWriter pw;

        ArrayList<PDS> PDS_list = new ArrayList<PDS>();
        for (int i =55;i <267; i++) {
            
            PDS_list = findPDS(i);
            System.out.println("num PDSs: " + PDS_list.size());

            pw = new PrintWriter(new FileWriter(new File("PDS18_output/pds18_64." + i + ".txt"))); //CHANGEPATH
            pw.write("pds_list := [");
            for (int j = 0; j < PDS_list.size()-1; j++) {
                pw.write(PDS_list.get(j).toString());
                pw.write(",\n");

            }
            if (PDS_list.size()==0) {
                pw.write("];\n");
            }
            else {
                pw.write(PDS_list.get(PDS_list.size()-1).toString());
                pw.write("];\n");

            }
            
            pw.close();
            
            System.out.println(java.util.Calendar.getInstance().getTime());
        }

        for (PDS p : PDS_list) {
            if (p.full) {
                System.out.println(p);
            }
        }
    }


    static void print2Darr(int[][] arr){
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr[0].length; j++) {
                System.out.print("" + arr[i][j] + " ");
            }
            System.out.println();
        }
    }

    static ArrayList<PDS> findPDS(int gid) throws FileNotFoundException {


        int[][] MT = new int[64][64];

        Scanner scan = new Scanner(new File("PDS18_gap_tables/gaptables64." + gid + ".txt"));
        StringTokenizer st = new StringTokenizer(scan.nextLine()); //first line is junk

        //System.out.println(st.nextToken());
        System.out.println(gid);


        //fill in possible parmsets
        st = new StringTokenizer(scan.nextLine());
        int numpossparmsets = Integer.parseInt(st.nextToken());

    
        int[][] possible_parmsets = new int[numpossparmsets][8];
        for (int i = 0; i < numpossparmsets; i++) {
            st = new StringTokenizer(scan.nextLine());
            for (int j = 0; j < 8; j++) {
                possible_parmsets[i][j] = Integer.parseInt(st.nextToken());
            }

        }

        //fill in Convolution table
        for (int i = 0; i < 64; i++) {
            st = new StringTokenizer(scan.nextLine());
            for (int j = 0; j < 64; j++) {
                MT[i][j] = Integer.parseInt(st.nextToken())-1;
            }
        }

        
        //fill in actual_choices_by_level table
        int[][][][] actual_choices_by_level = new int[8][9][70][8];
        int[][] actual_choices_num_choices = new int[8][9]; //store the number of choices at each level and parmset
        for (int i = 0; i < 8; i++) {

            st = new StringTokenizer(scan.nextLine());

            //System.out.println("i: " + i);

            for (int jj = 0; jj < 9; jj++) {


                //System.out.println("\tjj: " + jj);

                st = new StringTokenizer(scan.nextLine());

                int j = Integer.parseInt(st.nextToken());
                int choices = Integer.parseInt(st.nextToken());
                actual_choices_num_choices[i][j] = choices;

                for (int k = 0; k < choices; k++) {
                    //System.out.println("\t\tk: " + k);

                    st = new StringTokenizer(scan.nextLine());
                    if (st.hasMoreTokens()) {
                        for (int l = 0; l < j; l++) {
                            actual_choices_by_level[i][j][k][l] = Integer.parseInt(st.nextToken());
                        }
                    }
                    

                }

            }
            
            
            

        }

        // for (int i =0; i < 8; i++) {
        //     for (int j = 0;j< 9; j++) {
        //         System.out.println("level: " + i + " parm#: " + j);
        //         System.out.println(actual_choices_num_choices[i][j]);
        //         //System.out.println(actual_choices_by_level[i][j]);
        //         //print2Darr(actual_choices_by_level[i][j]);
        //     }
        // }

    
        int[] allsixes = new int[64];
        Arrays.fill(allsixes,6);
            

        //Start iterating over all the possible combinations:
        //counts how many parameter values have been checked
        int counter = 0;
        boolean goon;
        int[] conprod;
        int[] conv;
        int[] conprodmid;
        int intssofar;
        ArrayList<PDS> pdss = new ArrayList<PDS>();

        //Debugging: see if a working example actually works

        // int[] realpds = { 5, 20, 21, 22, 2, 11, 3, 16, 4, 17, 29, 51, 53, 63, 8, 61, 44, 59 };

        //  conprod = new int[64]; //got this from Dr. Smith & Dr. Kaliszewski

        // for (int i1 = 0; i1 < intssofar; i1++) {
        // for (int i2 = 0; i2 < intssofar; i2++) {
        // int i3 = MT[realpds[i1]-1][realpds[i2]-1]; //notice the -1 because of Java
        // conprod[i3] = conprod[i3] + 1;
        // }
        // }

        // int numtwoo = 0;
        // int numsixx = 0;
        // for (int i4 = 1; i4 < 64; i4++) {
        // if (conprod[i4]==2) numtwoo++;
        // if (conprod[i4]==6) numsixx++;
        // }
        // System.out.println(numtwoo);
        // System.out.println(numsixx);
        // System.out.println(Arrays.toString(conprod));

        for (int iii = 0; iii < possible_parmsets.length; iii++) {

            //System.out.println(iii);
            int[] parmset = possible_parmsets[iii];
            int[] temp = new int[18];
            int ti = 0; //ti stands for temp index

            //System.out.println(Arrays.toString(parmset));


            //note we start at 0 here but at 1 in GAP
            for (int a0 = 0; a0 < actual_choices_num_choices[0][parmset[0]]; a0++) {
                for (int k0 = 0; k0 < parmset[0]; k0++) temp[ti++] = actual_choices_by_level[0][parmset[0]][a0][k0];

                for (int a1 = 0; a1 < actual_choices_num_choices[1][parmset[1]]; a1++) {
                    for (int k1 = 0; k1 < parmset[1]; k1++) temp[ti++] = actual_choices_by_level[1][parmset[1]][a1][k1];

                    for (int a2 = 0; a2 < actual_choices_num_choices[2][parmset[2]]; a2++) {
                        for (int k2 = 0; k2 < parmset[2]; k2++) temp[ti++] = actual_choices_by_level[2][parmset[2]][a2][k2];

                        for (int a3 = 0; a3 < actual_choices_num_choices[3][parmset[3]]; a3++) {
                            for (int k3 = 0; k3 < parmset[3]; k3++) temp[ti++] = actual_choices_by_level[3][parmset[3]][a3][k3];

                            for (int a4 = 0; a4 < actual_choices_num_choices[4][parmset[4]]; a4++) {
                                for (int k4 = 0; k4 < parmset[4]; k4++) temp[ti++] = actual_choices_by_level[4][parmset[4]][a4][k4];

                                //System.out.println(Arrays.toString(temp));

                                goon = true;
                                conprodmid = new int[64];
                                intssofar = parmset[0] + parmset[1] + parmset[2] + parmset[3] + parmset[4];
                                for (int ii1 = 0; ii1 < intssofar; ii1++) {
                                    for (int ii2 = 0; ii2 < intssofar; ii2++) {
                                        int ii3 = MT[temp[ii1]-1][temp[ii2]-1]; //notice the -1 because of Java
                                        conprodmid[ii3] = conprodmid[ii3] + 1;
                                    }
                                }
                                for (int ii4 = 1; ii4 < 64; ii4++) {

                                    if (conprodmid[ii4] > 6)  {
                                        goon=false;
                                        break;
                                    }
                                }

                                if (!goon) {
                                    ti -= parmset[4]; //set the ti back
                                    continue;
                                }

                                for (int a5 = 0; a5 < actual_choices_num_choices[5][parmset[5]]; a5++) {
                                    for (int k5 = 0; k5 < parmset[5]; k5++) temp[ti++] = actual_choices_by_level[5][parmset[5]][a5][k5];

                                    goon = true;
                                    conprodmid = new int[64];
                                    intssofar = parmset[0] + parmset[1] + parmset[2] + parmset[3] + parmset[4] + parmset[5];
                                    for (int ii1 = 0; ii1 < intssofar; ii1++) {
                                    for (int ii2 = 0; ii2 < intssofar; ii2++) {
                                    int ii3 = MT[temp[ii1]-1][temp[ii2]-1]; //notice the -1 because of Java
                                    conprodmid[ii3] = conprodmid[ii3] + 1;
                                    }
                                    }
                                    for (int ii4 = 1; ii4 < 64; ii4++) {

                                    if (conprodmid[ii4] > 6)  {
                                    goon=false;
                                    break;
                                    }
                                    }

                                    if (!goon) {
                                    ti = ti-parmset[5]; //set the ti back
                                    continue;
                                    }

                                    for (int a6 = 0; a6 < actual_choices_num_choices[6][parmset[6]]; a6++) {
                                        for (int k6 = 0; k6 < parmset[6]; k6++) temp[ti++] = actual_choices_by_level[6][parmset[6]][a6][k6];

                                        goon = true;
                                        conprodmid = new int[64];
                                        intssofar = parmset[0] + parmset[1] + parmset[2] + parmset[3] + parmset[4] + parmset[5] + parmset[6];
                                        for (int ii1 = 0; ii1 < intssofar; ii1++) {
                                        for (int ii2 = 0; ii2 < intssofar; ii2++) {
                                        int ii3 = MT[temp[ii1]-1][temp[ii2]-1]; //notice the -1 because of Java
                                        conprodmid[ii3] = conprodmid[ii3] + 1;
                                        }
                                        }
                                        for (int ii4 = 1; ii4 < 64; ii4++) {

                                        if (conprodmid[ii4] > 6)  {
                                        goon=false;
                                        break;
                                        }
                                        }

                                        if (!goon) {
                                        ti = ti-parmset[6]; //set the ti back
                                        continue;
                                        }


                                        for (int a7 = 0; a7 < actual_choices_num_choices[7][parmset[7]]; a7++) {
                                            for (int k7 = 0; k7 < parmset[7]; k7++) temp[ti++] = actual_choices_by_level[7][parmset[7]][a7][k7];

                                            conprod = new int[64]; //got this from Dr. Smith & Dr. Kaliszewski

                                            for (int i1 = 0; i1 < 18; i1++) {
                                            for (int i2 = 0; i2 < 18; i2++) {
                                            int i3 = MT[temp[i1]-1][temp[i2]-1]; //notice the -1 because of Java
                                            conprod[i3] = conprod[i3] + 1;
                                            }
                                            }

                                            int numtwo = 0;
                                            int numsix = 0;
                                            for (int i4 = 1; i4 < 64; i4++) {
                                            if (conprod[i4]==2) numtwo++;
                                            if (conprod[i4]==6) numsix++;
                                            }
                                            if (conprod[0]==18 && numtwo==18 && numsix==45) {
                                                //System.out.println(temp);

                                            int[] listy = new int[64];
                                            int[] listy2 = new int[64];

                                            for (int ppos = 0; ppos < 18; ppos++) {
                                                listy[temp[ppos]-1] += 1;
                                                listy2[temp[ppos]-1] += 1;
                                            }

                                            //Add the identity
                                            listy[0] -= 2;
                                            listy2[0] += 6;
                                            conv = new int[64];

                                            for (int i1 = 0; i1 < 64; i1++) {
                                            for (int i2 = 0; i2 < 64; i2++) {
                                            int i3 = MT[i1][i2];
                                            conv[i3] += listy[i1] * listy2[i2];
                                            }
                                            }

                                            if (Arrays.equals(conv,allsixes)) {

                                            //System.out.println("Found a PDS!");
                                            //System.out.println(Arrays.toString(temp));
                                            //System.out.println(counter);
                                            pdss.add(new PDS(gid,temp));
                                            }
                                            }

                                            counter++;

                                            ti -= parmset[7];
                                        }
                                        ti -= parmset[6];

                                    }
                                    ti -=parmset[5];
                                }
                                ti -= parmset[4];
                            }
                            ti -= parmset[3];
                        }

                        ti -= parmset[2];
                    }

                    ti -= parmset[1];
                }
                ti -= parmset[0];
            }
            
        }
        //System.out.println("No PDSs found, # looked through: " + counter);

        scan.close();

        System.out.println("counter: " + counter);

        return pdss;
    }

}