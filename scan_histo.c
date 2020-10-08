#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>



int main(int argc, char *argv[])
{
    
 char filename[100];
 char string[400];
 int i,n,m,k,l;
 char character_string[50];   
 char *token,*cp;
 char delimiter[] = ",";
 char aminoacid1[3],aminoacid2[3];
 char stringala[3];
 float population[15][22][22][201];
 float population_last[15][22][22];
 float counter3[15][22][22];
 
 float dist[15][22][22],ranger=0.0; 
 int number_of_resid=0,boolean=0;
 char sam[5];
 float phi_dist1[22][22][201];
 float phi_dist2[22][22][201]; 
 float phi_dist3[22][22][201];
 float phi_dist4[22][22][201];
 float counter=0.0,counter2=0.0;
 float ranger_l,ranger_k;
 float d_t[22][22],d_t1;
 float maximum;
 int   maximum_n[422],maximum_m[422];
 char char1[1],char2[10]; 
 int   count;

 FILE *fp2;
 
 FILE *fp = fopen("histogram_matrix","r");
          
 
 n = 1;
 m = 1;
 k = 1;
 
 while(fgets(string,400,fp) != 0) {
 
 
         count = 0;

        for(l=1;l<=70;l++) {
 
           if(isspace(string[l])){

               count++;

           }; 

          };
          
          printf("%d\n",count);
          
            if(count == 4) {
                
     
                    if(strcmp(char1,"#")!=0) 

                     {

                      sscanf(string,"%e%e%e%e%e%e%e%e%e%e%e",
                                                 &ranger,&population[1][n][m][k],&population[2][n][m][k],
                                                 &population[3][n][m][k],&population[4][n][m][k],
                                                 &population[5][n][m][k],&population[6][n][m][k],
                                                 &population[7][n][m][k],&population[8][n][m][k],
                                                 &population[9][n][m][k],&population[10][n][m][k]);
                      
                     printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                                                 ranger,population[1][n][m][k],population[2][n][m][k],
                                                 population[3][n][m][k],population[4][n][m][k],
                                                 population[5][n][m][k],population[6][n][m][k],
                                                 population[7][n][m][k],population[8][n][m][k],
                                                 population[9][n][m][k],population[10][n][m][k]                         
                    ); 
                     
                    printf("%d%s\n",k," k ");
 
                    k ++;
                    
                    if(k>200) {
                        
                        k=1;
                        m++;

                        printf("%d\t%d\n",m,k," m, k ");

                    };
                    
                    if(m>21) {
                        
                        m=1;
                        n++;

                      printf("%d\n",n);
                        
                    };
                    
               };
  
            };      
               
};
 
 fclose(fp);
 
 for(n=1;n<=21;n++){
  
     for(m=1;m<=21;m++){
         
         for(k=150;k<=175;k++){
             
            for(l=1;l<=10;l++){
   
               if(population[l][n][m][k] > 0.0) population_last[l][n][m] = (population[l][n][m][k]);
                
            };
            
         };
     };
 };
 
 for(n=1;n<=21;n++){
  
     for(m=1;m<=21;m++){
         
         for(k=1;k<=200;k++){
          
            for(l=1;l<=10;l++){

                if(population[l][n][m][k] > 0.0) population[l][n][m][k] = -log(population[l][n][m][k]/population_last[l][n][m]);
                
               };

             };   

    };
     
  };
 
 fp = fopen("histogram_matrix_pmf","w");
 fp2 = fopen("histogram_matrix_rdf","w");
 
           for(n=1;n<=21;n++){ 
               
              if(n==1) fprintf(fp,"%s\n","#ALA"); 
              if(n==2) fprintf(fp,"%s\n","#VAL");
              if(n==3) fprintf(fp,"%s\n","#ILE");
              if(n==4) fprintf(fp,"%s\n","#LEU");              
              if(n==5) fprintf(fp,"%s\n","#PHE");
              if(n==6) fprintf(fp,"%s\n","#TYR"); 
              if(n==7) fprintf(fp,"%s\n","#TRP");
              if(n==8) fprintf(fp,"%s\n","#ASN");
              if(n==9) fprintf(fp,"%s\n","#GLN");              
              if(n==10) fprintf(fp,"%s\n","#ASP");              
              if(n==11) fprintf(fp,"%s\n","#GLU"); 
              if(n==12) fprintf(fp,"%s\n","#MET");
              if(n==13) fprintf(fp,"%s\n","#THR");
              if(n==14) fprintf(fp,"%s\n","#SER");              
              if(n==15) fprintf(fp,"%s\n","#PRO");
              if(n==16) fprintf(fp,"%s\n","#HIS"); 
              if(n==17) fprintf(fp,"%s\n","#ARG");
              if(n==18) fprintf(fp,"%s\n","#LYS");
              if(n==19) fprintf(fp,"%s\n","#CYS");              
              if(n==20) fprintf(fp,"%s\n","#SEC");  
              if(n==21) fprintf(fp,"%s\n","#GLY");

            for(m=1;m<=21;m++){
                   
              if(m==1) fprintf(fp,"%s\n","#ALA"); 
              if(m==2) fprintf(fp,"%s\n","#VAL");
              if(m==3) fprintf(fp,"%s\n","#ILE");
              if(m==4) fprintf(fp,"%s\n","#LEU");              
              if(m==5) fprintf(fp,"%s\n","#PHE");
              if(m==6) fprintf(fp,"%s\n","#TYR"); 
              if(m==7) fprintf(fp,"%s\n","#TRP");
              if(m==8) fprintf(fp,"%s\n","#ASN");
              if(m==9) fprintf(fp,"%s\n","#GLN");              
              if(m==10) fprintf(fp,"%s\n","#ASP");              
              if(m==11) fprintf(fp,"%s\n","#GLU");   
              if(m==12) fprintf(fp,"%s\n","#MET");
              if(m==13) fprintf(fp,"%s\n","#THR");
              if(m==14) fprintf(fp,"%s\n","#SER");              
              if(m==15) fprintf(fp,"%s\n","#PRO");
              if(m==16) fprintf(fp,"%s\n","#HIS"); 
              if(m==17) fprintf(fp,"%s\n","#ARG");
              if(m==18) fprintf(fp,"%s\n","#LYS");
              if(m==19) fprintf(fp,"%s\n","#CYS");              
              if(m==20) fprintf(fp,"%s\n","#SEC");  
              if(m==21) fprintf(fp,"%s\n","#GLY");                   
               
                 fprintf(fp,"%s\n","####C-N, C-CA, C-C, C-O, N-CA, N-N, N-O, CA-CA, O-O, O-CA");
                   
                 ranger = 0.0;  
                   
                 for(k=1;k<=200;k++) {  
                   
                   ranger = ranger + 25.0/200.0;
                   
                   if(k < 175.0) {
                 
                   fprintf(fp,"%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n",
                                                 ranger,population[1][n][m][k],population[2][n][m][k],
                                                 population[3][n][m][k],population[4][n][m][k],
                                                 population[5][n][m][k],population[6][n][m][k],
                                                 population[7][n][m][k],population[8][n][m][k],
                                                 population[9][n][m][k],population[10][n][m][k]);
                   
                   } else {
                       
                     
                   fprintf(fp,"%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n",
                                                 ranger,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);                       
                       
                   };
                   
                 };
                 
                 fprintf(fp,"%s\n","       ");
                 
            };
            
           };
 
 fclose(fp);
  
           for(n=1;n<=21;n++){ 
               
              if(n==1) fprintf(fp2,"%s\n","#ALA"); 
              if(n==2) fprintf(fp2,"%s\n","#VAL");
              if(n==3) fprintf(fp2,"%s\n","#ILE");
              if(n==4) fprintf(fp2,"%s\n","#LEU");              
              if(n==5) fprintf(fp2,"%s\n","#PHE");
              if(n==6) fprintf(fp2,"%s\n","#TYR"); 
              if(n==7) fprintf(fp2,"%s\n","#TRP");
              if(n==8) fprintf(fp2,"%s\n","#ASN");
              if(n==9) fprintf(fp2,"%s\n","#GLN");              
              if(n==10) fprintf(fp2,"%s\n","#ASP");              
              if(n==11) fprintf(fp2,"%s\n","#GLU"); 
              if(n==12) fprintf(fp2,"%s\n","#MET");
              if(n==13) fprintf(fp2,"%s\n","#THR");
              if(n==14) fprintf(fp2,"%s\n","#SER");              
              if(n==15) fprintf(fp2,"%s\n","#PRO");
              if(n==16) fprintf(fp2,"%s\n","#HIS"); 
              if(n==17) fprintf(fp2,"%s\n","#ARG");
              if(n==18) fprintf(fp2,"%s\n","#LYS");
              if(n==19) fprintf(fp2,"%s\n","#CYS");              
              if(n==20) fprintf(fp2,"%s\n","#SEC");  
              if(n==21) fprintf(fp2,"%s\n","#GLY");

            for(m=1;m<=21;m++){
                   
              if(m==1) fprintf(fp2,"%s\n","#ALA"); 
              if(m==2) fprintf(fp2,"%s\n","#VAL");
              if(m==3) fprintf(fp2,"%s\n","#ILE");
              if(m==4) fprintf(fp2,"%s\n","#LEU");              
              if(m==5) fprintf(fp2,"%s\n","#PHE");
              if(m==6) fprintf(fp2,"%s\n","#TYR"); 
              if(m==7) fprintf(fp2,"%s\n","#TRP");
              if(m==8) fprintf(fp2,"%s\n","#ASN");
              if(m==9) fprintf(fp2,"%s\n","#GLN");              
              if(m==10) fprintf(fp2,"%s\n","#ASP");              
              if(m==11) fprintf(fp2,"%s\n","#GLU");   
              if(m==12) fprintf(fp2,"%s\n","#MET");
              if(m==13) fprintf(fp2,"%s\n","#THR");
              if(m==14) fprintf(fp2,"%s\n","#SER");              
              if(m==15) fprintf(fp2,"%s\n","#PRO");
              if(m==16) fprintf(fp2,"%s\n","#HIS"); 
              if(m==17) fprintf(fp2,"%s\n","#ARG");
              if(m==18) fprintf(fp2,"%s\n","#LYS");
              if(m==19) fprintf(fp2,"%s\n","#CYS");              
              if(m==20) fprintf(fp2,"%s\n","#SEC");  
              if(m==21) fprintf(fp2,"%s\n","#GLY");                   
               
                 fprintf(fp2,"%s\n","####C-N, C-CA, C-C, C-O, N-CA, N-N, N-O, CA-CA, O-O, O-CA");
                   
                 ranger = 0.0;  
                   
             for(k=1;k<=200;k++){
          
                if(k < 175) { 
                 
                  for(l=1;l<=10;l++){

                    population[l][n][m][k] = exp(-population[l][n][m][k]);
                    
                  };
                  
                } else {
                    
                    population[l][n][m][k] = 1.0;                    
                    
                };
                  
             };
              
                 for(k=1;k<=200;k++) {  
                     
                   ranger = ranger + 25.0/200.0;
                   
                   fprintf(fp2,"%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n",
                                                 ranger,population[1][n][m][k],population[2][n][m][k],
                                                 population[3][n][m][k],population[4][n][m][k],
                                                 population[5][n][m][k],population[6][n][m][k],
                                                 population[7][n][m][k],population[8][n][m][k],
                                                 population[9][n][m][k],population[10][n][m][k]);
                   
                 };
                 
                 fprintf(fp2,"%s\n","       ");
                 
            };
            
           };
 
 fclose(fp2); 
 
 return 0;   
}
