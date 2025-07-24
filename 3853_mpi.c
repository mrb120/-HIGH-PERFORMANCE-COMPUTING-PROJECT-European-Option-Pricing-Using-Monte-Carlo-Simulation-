//EE7218   :High Performance Computing
//Name     :Batuwatta M.R.
//Reg. No. :EG/2020/353
//Project Title :European Option Pricing Using Monte Carlo Simulation

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


#define MAX_DAYS 252
#define NUM_SIMULATIONS 5000000
#define LINE_LENGTH 256
#define NO_OF_COMPANIES 7

double GenerateStandardNormal(double mean,double stnddev,unsigned int*seed) {
    double u1,u2,Z;
    do{
        u1=((double)rand_r(seed))/RAND_MAX;
    } while(u1<=1e-12);
    u2=(double)rand_r(seed)/RAND_MAX;
    // Box_muller method
    Z=sqrt(-2.0*log(u1))*cos(2.0*3.14*u2);
    return mean+Z*stnddev;
}

typedef struct {
    double sigma;
    double S0;
} Result;

Result CalculateSigmaAndS0(char* filename,int company_no){
    Result res;
    FILE *CSV =fopen(filename,"r");
    if (!CSV){
        printf("Error opening CSV file!!");
        exit(1);
    }

    char line[LINE_LENGTH];
    double prices[MAX_DAYS];
    int count=0;
    double S0[NO_OF_COMPANIES];

    fgets(line,sizeof(line),CSV);
    
    while (fgets(line,sizeof(line),CSV) !=NULL)  {
        char temp_line[LINE_LENGTH];
        strcpy(temp_line,line);
        char* rest=temp_line;

        int current_column =0;
        char* token=strtok_r(line,",",&rest);

        while (token !=NULL){
            if(current_column==4 && count ==0){
                S0[company_no]=atof(token);
            }
            if(current_column == 4){
                //printf("Close price: %s\n", token);
                prices[count++]=atof(token);
                break;
            }
        token=strtok_r(NULL,",",&rest);
        current_column++;
        }
     } 
fclose (CSV);

//Daily log returns
double returns[MAX_DAYS];
int retCount=0;

for (int i=1;i<count;i++){
    returns[retCount++]=log(prices[i]/prices[i-1]);
}

// Compute mean of returns
double sum=0;
for (int i=0;i<retCount;i++) {
    sum +=returns[i];
}
double mean=sum/retCount;

    double sq_sum=0;
    for (int i=0;i<retCount;i++) {
        sq_sum +=(returns[i]-mean)*(returns[i]-mean);
    }
    double daily_std=sqrt(sq_sum/(retCount-1));
    res.sigma =daily_std * sqrt(252.0);
    res.S0=S0[company_no];
return res;  
}


int main(int argc,char *argv[]) {
    int numtasks,rank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    
    double mean =0.0;
    double stddev =1.0;
    double T=1.0;
    double r=0.0237;
    char* tickers[NO_OF_COMPANIES]={"TESLA","APPLE","GOOGLE","INTEL","MICROSOFT","NVIDIA","COCACOLA"};

    char *CSV_FILE[NO_OF_COMPANIES]={ "Download Data - STOCK_US_XNAS_TSLA.csv",
                                      "Download Data - STOCK_US_XNAS_AAPL.csv",
                                      "Download Data - STOCK_US_XNAS_GOOG.csv",
                                      "Download Data - STOCK_US_XNAS_INTC.csv",
                                      "Download Data - STOCK_US_XNAS_MSFT.csv",
                                      "Download Data - STOCK_US_XNAS_NVDA.csv",
                                      "Download Data - STOCK_US_XNYS_KO.csv"};
    double S0[NO_OF_COMPANIES];
    double option_prices[NO_OF_COMPANIES];

    double StartTime= MPI_Wtime();;
    

    for (int i=0; i<NO_OF_COMPANIES;i++){
        Result result;
        if (rank==0){
             result=CalculateSigmaAndS0(CSV_FILE[i],i);}

        MPI_Bcast(&result.sigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&result.S0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        int sims_per_process=NUM_SIMULATIONS/numtasks;
        int remaining_sims=NUM_SIMULATIONS % numtasks;
        int my_simulations=sims_per_process;

        if (rank==numtasks-1) {
            my_simulations +=remaining_sims;  
        }
        double TotalPayoff=0.0;
        for (int j=0;j<my_simulations;j++){
            unsigned int seed =(unsigned int)(i^j);
            double K=result.S0 *0.95;
            double z=GenerateStandardNormal(mean,stddev,&seed);
            double S_T=result.S0 *exp((r-0.5*result.sigma*result.sigma)*T+(result.sigma*sqrt(T)*z));
            double Payoff=fmax((S_T-K),0.0);
            TotalPayoff +=Payoff;
        }
        double GlobalTotalPayoff = 0.0;
        MPI_Reduce(&TotalPayoff, &GlobalTotalPayoff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if(rank==0){
        double AveragePayoff = GlobalTotalPayoff /NUM_SIMULATIONS;
        double OptionPrice = AveragePayoff * exp(-r*T);
        option_prices[i] = OptionPrice;
         }
}
    
    double EndTime=MPI_Wtime();
    double exeTime=EndTime-StartTime;
     
     
    if(rank==0){
    for (int i=0;i<NO_OF_COMPANIES;i++){
        printf("%s: Option price: %.2f$\n", tickers[i], option_prices[i]);
    }
    printf("\nTime to execute: %.2fs\n",exeTime);  
    }
    MPI_Finalize();
   return 0;
}