#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>



int main(int argc, char *argv[])
{
    
 char filename[100];
 char string[200];
 int i,n,m,k,l;
 char character_string[50];   
 char *token,*cp;
 char delimiter[] = ",";
 char aminoacid1[3],aminoacid2[3];
 char stringala[3];
 float population[15][22][22][201];
 float counter3[15][22][22];
 
 float dist[15][22][22],ranger=0.0; 
 int number_of_resid=0,boolean=0;
 char sam[5];
 float counter=0.0,counter2=0.0;
 float ranger_l,ranger_k;
 float d_t[22][22],d_t1;
 float maximum;
 int   maximum_n[422],maximum_m[422];
 char char1[1],char10[10]; 
 int   natoms;
 int   res_id[1000],at_id[1000];
 char  res_char[3][1000];
 char  at_char[3][1000];
 float x,y,z; 
 int   r,count=0;
 char  rdf_name[200];
 
 
 FILE *fp = fopen("histogram_matrix_pmf","r");
          
 
 n = 1;
 m = 1;
 k = 1;
 
 while(fgets(string,200,fp) != 0) {
 
 
         count = 1;

         sscanf(string,"%1s",char1);
         
         sscanf(string,"%10s",char10);
         
         if(strcmp(char1,"#") == 0) count = 0;
         if(strcmp(char10,"         ") == 0) count = 0; 
         
          
            if(count == 1) {               
     

                      sscanf(string,"%e%e%e%e%e%e%e%e%e%e%e",
                                                 &ranger,&population[1][n][m][k],&population[2][n][m][k],
                                                 &population[3][n][m][k],&population[4][n][m][k],
                                                 &population[5][n][m][k],&population[6][n][m][k],
                                                 &population[7][n][m][k],&population[8][n][m][k],
                                                 &population[9][n][m][k],&population[10][n][m][k]);
                      
                   /*  printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                                                 ranger,population[1][n][m][k],population[2][n][m][k],
                                                 population[3][n][m][k],population[4][n][m][k],
                                                 population[5][n][m][k],population[6][n][m][k],
                                                 population[7][n][m][k],population[8][n][m][k],
                                                 population[9][n][m][k],population[10][n][m][k]); */ 
                    
//                    printf("%d\t%d\n",k,count); 
                    
                    k++;

                    if(k==201) {
                        
                        k=0;
                        m++;
                    
//                    printf("%d%s\n",m,"m");
                        
                    };
                    
                    if(m>21) {
                        
                        m=1;
                        n++;

//                      printf("%d%s\n",n,"n");
                        
                    };
  
            };      
               
};
 
 fclose(fp);
 
 fp = fopen("start.gro","r");
 
 k = 0;
 l = 1;
 
 while(fgets(string,48,fp) !=0) {
     
     k++;
     
    if(k==2) sscanf(string,"%d",&natoms);
     
    if(k>2 && k <= natoms+2) {
    
        sscanf(string,"%d%s%s%d%f%f%f",&res_id[l],res_char[l],at_char[l],&at_id[l],&x,&y,&z);
   //     printf("%d\t%s\t%s\t%d\t%e\t%e\t%e\n",res_id[l],res_char[l],at_char[l],at_id[l],x,y,z); 
        
        l++;
        
    };
     
 };
 
 fclose(fp);
 
 fp = fopen("atom_index.dat","w");

for(k=1;k<=natoms;k++) {
  
  if(strcmp(at_char[k],"N")==0 || 
     strcmp(at_char[k],"C")==0 ||
     strcmp(at_char[k],"CA")==0 ||
     strcmp(at_char[k],"O")==0 ||
     strcmp(at_char[k],"N")==0) {
         
      
      fprintf(fp,"%d\n",k);
      
    };
    
};
 
 fclose(fp); 
 
for(k=1;k<=natoms;k++) {
  
      for(l=1;l<=natoms;l++) {       
             
            
        strcpy(aminoacid1,res_char[k]);
        strcpy(aminoacid2,res_char[l]);        
        
        strcpy(stringala,"ALA");   
        if(strcmp(stringala,aminoacid1) == 0) n = 1;
        if(strcmp(stringala,aminoacid2) == 0) m = 1;
        strcpy(stringala,"VAL");   
        if(strcmp(stringala,aminoacid1) == 0) n = 2;
        if(strcmp(stringala,aminoacid2) == 0) m = 2;
        strcpy(stringala,"ILE");   
        if(strcmp(stringala,aminoacid1) == 0) n = 3;
        if(strcmp(stringala,aminoacid2) == 0) m = 3;
        strcpy(stringala,"LEU");   
        if(strcmp(stringala,aminoacid1) == 0) n = 4;
        if(strcmp(stringala,aminoacid2) == 0) m = 4;
        strcpy(stringala,"PHE");   
        if(strcmp(stringala,aminoacid1) == 0) n = 5;
        if(strcmp(stringala,aminoacid2) == 0) m = 5;
        strcpy(stringala,"TYR");   
        if(strcmp(stringala,aminoacid1) == 0) n = 6;
        if(strcmp(stringala,aminoacid2) == 0) m = 6;
        strcpy(stringala,"TRP");   
        if(strcmp(stringala,aminoacid1) == 0) n = 7;
        if(strcmp(stringala,aminoacid2) == 0) m = 7;        
        strcpy(stringala,"ASN");   
        if(strcmp(stringala,aminoacid1) == 0) n = 8;
        if(strcmp(stringala,aminoacid2) == 0) m = 8;
        strcpy(stringala,"GLN");   
        if(strcmp(stringala,aminoacid1) == 0) n = 9;
        if(strcmp(stringala,aminoacid2) == 0) m = 9;        
        strcpy(stringala,"ASP");   
        if(strcmp(stringala,aminoacid1) == 0) n = 10;
        if(strcmp(stringala,aminoacid2) == 0) m = 10;
        strcpy(stringala,"GLU");   
        if(strcmp(stringala,aminoacid1) == 0) n = 11;
        if(strcmp(stringala,aminoacid2) == 0) m = 11;         
        strcpy(stringala,"MET");   
        if(strcmp(stringala,aminoacid1) == 0) n = 12;
        if(strcmp(stringala,aminoacid2) == 0) m = 12;
        strcpy(stringala,"THR");   
        if(strcmp(stringala,aminoacid1) == 0) n = 13;
        if(strcmp(stringala,aminoacid2) == 0) m = 13;       
        strcpy(stringala,"SER");   
        if(strcmp(stringala,aminoacid1) == 0) n = 14;
        if(strcmp(stringala,aminoacid2) == 0) m = 14;  
        strcpy(stringala,"PRO");   
        if(strcmp(stringala,aminoacid1) == 0) n = 15;
        if(strcmp(stringala,aminoacid2) == 0) m = 15;         
        strcpy(stringala,"HIS");   
        if(strcmp(stringala,aminoacid1) == 0) n = 16;
        if(strcmp(stringala,aminoacid2) == 0) m = 16;          
        strcpy(stringala,"ARG");   
        if(strcmp(stringala,aminoacid1) == 0) n = 17;
        if(strcmp(stringala,aminoacid2) == 0) m = 17;         
        strcpy(stringala,"LYS");   
        if(strcmp(stringala,aminoacid1) == 0) n = 18;
        if(strcmp(stringala,aminoacid2) == 0) m = 18;          
        strcpy(stringala,"CYS");   
        if(strcmp(stringala,aminoacid1) == 0) n = 19;
        if(strcmp(stringala,aminoacid2) == 0) m = 19;         
        strcpy(stringala,"SEC");   
        if(strcmp(stringala,aminoacid1) == 0) n = 20;
        if(strcmp(stringala,aminoacid2) == 0) m = 20;           
        strcpy(stringala,"GLY");   
        if(strcmp(stringala,aminoacid1) == 0) n = 21;
        if(strcmp(stringala,aminoacid2) == 0) m = 21;
  // ##"####C-N, C-CA, C-C, C-O, N-CA, N-N, N-O, CA-CA, O-O, O-CA;  
        
       // printf("%s\t%s\n",at_char[k],at_char[l]);
       
        r = 0;
        
      if(k != l) { 
          
        if(strcmp(at_char[k],"N")==0 && strcmp(at_char[l],"C")==0) {
        r=1;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");
        
        };
        if(strcmp(at_char[k],"C")==0 && strcmp(at_char[l],"N")==0) {
        r=1;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");
        };
        if(strcmp(at_char[k],"C")==0 && strcmp(at_char[l],"CA")==0) {
            r=2;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"CA")==0 && strcmp(at_char[l],"C")==0) {
            r=2;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"C")==0 && strcmp(at_char[l],"C")==0) {
            r=3;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"C")==0 && strcmp(at_char[l],"O")==0) {
            r=4;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"O")==0 && strcmp(at_char[l],"C")==0) {
            r=4;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };        
        if(strcmp(at_char[k],"N")==0 && strcmp(at_char[l],"CA")==0) {
            r=5;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"CA")==0 && strcmp(at_char[l],"N")==0) {
            r=5;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"N")==0 && strcmp(at_char[l],"N")==0) {
            r=6;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"N")==0 && strcmp(at_char[l],"O")==0) {
            r=7;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"O")==0 && strcmp(at_char[l],"N")==0) {
            r=7;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"CA")==0 && strcmp(at_char[l],"CA")==0) {
            r=8;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"O")==0 && strcmp(at_char[l],"O")==0) {
            r=9;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"O")==0 && strcmp(at_char[l],"CA")==0) {
            r=10;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        if(strcmp(at_char[k],"CA")==0 && strcmp(at_char[l],"O")==0) {
            r=10;
        sprintf(rdf_name,"%s%d%d%s","rdf_",k,l,".dat");            
        };
        
        if(r != 0) {
        
        ranger = 0.0125;
        
        fp = fopen(rdf_name,"w");
    
          for(i=1;i<=200;i++) {
        
             fprintf(fp,"%e\t%e\n",ranger,population[r][n][m][i]);
        
             ranger = ranger + 2.5/200.0;
           
            };
           
          fclose(fp);
        
         };
        
       };
      
     };
     
  };
  
 return 0;   
}

