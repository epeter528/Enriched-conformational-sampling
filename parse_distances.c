#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>


int main(int argc, char *argv[])
{
    
 char filename[100];
 char string[1];
 int counter1=0,i,n,m,k;
 int vn,vm;
 char character_string[400];   
 char *token,*cp;
 char delimiter[] = ",";
 char aminoacid1[3],aminoacid2[3];
 char stringala[3];
 float population[15][22][22][201];
 float counter[15][22][22];
 float ranger_1,ranger_2;
 float PI=3.14159265359; 

 float dist[15][22][22],ranger=0.0;
 
// printf("%s\n","file-name with distances .csv");
    
// scanf("%s",&filename);
//   system("rm ./total.csv");
//   system("cat ./*.csv > total.csv");
   
    for(vn=1;vn<=21;vn++){ 
               
            for(vm=1;vm<=21;vm++){
        
                for(i=1;i<=14;i++){      
                   
                 for(k=1;k<=200;k++) {  
                     
                       population[i][vn][vm][k] = 0.0;
                       
                       counter[i][vn][vm] = 0.0;
                       
                     };
       
                };
                   
           };
            
        };  
   
 FILE *fp = fopen("test.csv","r");
 
    while(fgets(string,400,fp)!=0) {
  
    int count = 0;

      for(k=1;k<=18;k++) {
 
         if(isspace(string[k])){

             count++;

         }; 

       };

    printf("%d\n",count);
 
    if(count < 1) {
      
        token = strtok(string,delimiter);
        token = strtok(NULL,delimiter);  
        token = strtok(NULL,delimiter);

        sscanf(token,"%3s",aminoacid1);        
        
        token = strtok(NULL,delimiter);

        sscanf(token,"%3s",aminoacid2);
        
      //  printf("%s\t%s\n",aminoacid1,aminoacid2);
      //    printf("%s\n",stringala);
 
        m = 0;
        n = 0;

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
        strcpy(stringala,"HIS");   
        if(strcmp(stringala,aminoacid1) == 0) n = 20;
        if(strcmp(stringala,aminoacid2) == 0) m = 20;           
        strcpy(stringala,"GLY");   
        if(strcmp(stringala,aminoacid1) == 0) n = 21;
        if(strcmp(stringala,aminoacid2) == 0) m = 21;     
                
     //   printf("%s\t%s\t%d\t%d\n",aminoacid1,aminoacid2,n,m);
        
      if(strcmp("-",aminoacid1)!=0 && strcmp("-",aminoacid1)!=0 && m !=0 && n != 0) { 
        
         for(i=1;i<=14;i++){    
            
             token = strtok(NULL,delimiter);
        
             if(strcmp(token,"-")!=0) sscanf(token,"%f",&dist[i][n][m]);
             if(isnan(dist[i][n][m])==1) dist[i][n][m] = 0.0; 
          }; 
        
        };
          
         for(vn=1;vn<=21;vn++){ 
            
             if(vn == n) {  
               
            for(vm=1;vm<=21;vm++){
               
              if(vm == m) {  
        
                for(i=1;i<=14;i++){      
                     
                 ranger = 0.0;  
                   
                 for(k=1;k<=200;k++) {  
                   
                   ranger = ranger + 10.0/200.0; 
                     
                   if(ranger-0.5*10.0/200.0 <= dist[i][vn][vm] && dist[i][vn][vm] <= ranger+0.5*10.0/300.0) {
                     
                       population[i][vn][vm][k] = population[i][vn][vm][k] + 1.0;
                       
                       counter[i][vn][vm] = counter[i][vn][vm] + 1.0;
                     };
       
                   };
                   
                 };
       
               };
               
             };
       
           };
            
        };         
       
    };

    counter1++;

    if(counter1%10 == 0) printf("%d\n",counter1);

    if(counter1 > 284289570) break;
 
    };
 
 fclose(fp);
 
        for(i=1;i<=14;i++){    
          
           for(n=1;n<=21;n++){ 
            
               for(m=1;m<=21;m++){  
         
                ranger = 0.0;
          
                 for(k=1;k<=200;k++) {  
 
                    ranger = ranger + 10.0/200.0;

                    ranger_2 = ranger;
                    ranger_1 = ranger-10.0/200.0;

                    if(counter[i][n][m] > 0.0) population[i][n][m][k] = population[i][n][m][k]/(4.0/3.0*PI*(pow(ranger_2,3)-pow(ranger_1,3)))/counter[i][n][m]/(4.0/3.0*pow(10.0,3)*PI); 
                    
                 };
               };
           };
           
        };
        
 
 fp = fopen("histogram_matrix","w");
          
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
              if(n==20) fprintf(fp,"%s\n","#HIS");  
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
              if(m==20) fprintf(fp,"%s\n","#HIS");  
              if(m==21) fprintf(fp,"%s\n","#GLY");                   
               
                 fprintf(fp,"%s\n","##C-N C-CA C-CB C-C C-N N-CA N-CB N-C N-N C-CA C-CB C-CG N-CG CA-CG");
                   
                 ranger = 0.0;  
                   
                 for(k=1;k<=200;k++) {  
                   
                   ranger = ranger + 10.0/200.0;
                 
                   fprintf(fp,"%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n",
                                                 ranger,population[1][n][m][k],population[2][n][m][k],
                                                 population[3][n][m][k],population[4][n][m][k],
                                                 population[5][n][m][k],population[6][n][m][k],
                                                 population[7][n][m][k],population[8][n][m][k],
                                                 population[9][n][m][k],population[10][n][m][k],
                                                 population[11][n][m][k],population[12][n][m][k],
                                                 population[13][n][m][k],population[14][n][m][k]);
                   
                 };
                 
            };
            
           };
 
 fclose(fp);
    
 return 0;   
}
