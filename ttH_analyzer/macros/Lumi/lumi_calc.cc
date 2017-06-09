#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>

using namespace std;

int main()
{

    ifstream fin;
    int n_e, n_mu1, n_mu2, n_ee, n_emu1, n_emu2, n_emu3, n_emu4, n_mumu1, n_mumu2, n_mumu3, n_mumu4;
    n_e = n_mu1 = n_mu2 = n_ee = n_emu1 = n_emu2 = n_emu3 = n_emu4 = n_mumu1 = n_mumu2 = n_mumu3 = n_mumu4 = 0;

    double sl_e_run[10000], sl_e_lumi[10000];
    double sl_mu1_run[10000], sl_mu1_lumi[10000];
    double sl_mu2_run[10000], sl_mu2_lumi[10000];
    double dl_ee_run[10000], dl_ee_lumi[10000];
    double dl_emu1_run[10000], dl_emu1_lumi[10000];
    double dl_emu2_run[10000], dl_emu2_lumi[10000];
    double dl_emu3_run[10000], dl_emu3_lumi[10000];
    double dl_emu4_run[10000], dl_emu4_lumi[10000];
    double dl_mumu1_run[10000], dl_mumu1_lumi[10000];
    double dl_mumu2_run[10000], dl_mumu2_lumi[10000];
    double dl_mumu3_run[10000], dl_mumu3_lumi[10000];
    double dl_mumu4_run[10000], dl_mumu4_lumi[10000];


    system("awk '{print $1, $10}' lumi_sl_e.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>sl_e_run[n_e]>>sl_e_lumi[n_e];
        n_e++;
    }
    n_e--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_sl_mu1.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>sl_mu1_run[n_mu1]>>sl_mu1_lumi[n_mu1];
        n_mu1++;
    }
    n_mu1--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_sl_mu2.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>sl_mu2_run[n_mu2]>>sl_mu2_lumi[n_mu2];
        n_mu2++;
    }
    n_mu2--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_ee.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_ee_run[n_ee]>>dl_ee_lumi[n_ee];
        n_ee++;
    }
    n_ee--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_emu1.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_emu1_run[n_emu1]>>dl_emu1_lumi[n_emu1];
        n_emu1++;
    }
    n_emu1--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_emu2.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_emu2_run[n_emu2]>>dl_emu2_lumi[n_emu2];
        n_emu2++;
    }
    n_emu2--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_emu3.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_emu3_run[n_emu3]>>dl_emu3_lumi[n_emu3];
        n_emu3++;
    }
    n_emu3--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_emu4.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_emu4_run[n_emu4]>>dl_emu4_lumi[n_emu4];
        n_emu4++;
    }
    n_emu4--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_mumu1.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_mumu1_run[n_mumu1]>>dl_mumu1_lumi[n_mumu1];
        n_mumu1++;
    }
    n_mumu1--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_mumu2.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_mumu2_run[n_mumu2]>>dl_mumu2_lumi[n_mumu2];
        n_mumu2++;
    }
    n_mumu2--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_mumu3.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_mumu3_run[n_mumu3]>>dl_mumu3_lumi[n_mumu3];
        n_mumu3++;
    }
    n_mumu3--;
    fin.close();
    system("rm -rf out.txt");

    system("awk '{print $1, $10}' lumi_dl_mumu4.txt > out.txt");
    fin.open("out.txt");
    while(!fin.eof()){
        fin>>dl_mumu4_run[n_mumu4]>>dl_mumu4_lumi[n_mumu4];
        n_mumu4++;
    }
    n_mumu4--;
    fin.close();
    system("rm -rf out.txt");

    // Single Electron
    std::cout<<"No. of fills for SL e trigger: "<<n_e<<"\n";
    double sl_e_lumi_sum = 0;
    for(int i=0; i<n_e; i++){
        sl_e_lumi_sum = sl_e_lumi_sum + sl_e_lumi[i];
    }
    std::cout<<"Integrated Luminosity for SL e channel : "<<sl_e_lumi_sum/1000<<"\n\n";

    /*
    std::cout<<"No. of fills for SL mu trigger 1: "<<n_mu1<<"\n";
    double sl_mu1_lumi_sum = 0;
    for(int i=0; i<n_mu1; i++){
        sl_mu1_lumi_sum = sl_mu1_lumi_sum + sl_mu1_lumi[i];
    }
    std::cout<<"Integrated Luminosity for SL mu trigger 1 : "<<sl_mu1_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for SL mu trigger 2: "<<n_mu2<<"\n";
    double sl_mu2_lumi_sum = 0;
    for(int i=0; i<n_mu2; i++){
        sl_mu2_lumi_sum = sl_mu2_lumi_sum + sl_mu2_lumi[i];
    }
    std::cout<<"Integrated Luminosity for SL mu trigger 2 : "<<sl_mu2_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL ee trigger: "<<n_ee<<"\n";
    double dl_ee_lumi_sum = 0;
    for(int i=0; i<n_ee; i++){
        dl_ee_lumi_sum = dl_ee_lumi_sum + dl_ee_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL ee trigger : "<<dl_ee_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL emu trigger 1: "<<n_emu1<<"\n";
    double dl_emu1_lumi_sum = 0;
    for(int i=0; i<n_emu1; i++){
        dl_emu1_lumi_sum = dl_emu1_lumi_sum + dl_emu1_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL emu trigger 1 : "<<dl_emu1_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL emu trigger 2: "<<n_emu2<<"\n";
    double dl_emu2_lumi_sum = 0;
    for(int i=0; i<n_emu2; i++){
        dl_emu2_lumi_sum = dl_emu2_lumi_sum + dl_emu2_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL emu trigger 2 : "<<dl_emu2_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL emu trigger 3: "<<n_emu3<<"\n";
    double dl_emu3_lumi_sum = 0;
    for(int i=0; i<n_emu3; i++){
        dl_emu3_lumi_sum = dl_emu3_lumi_sum + dl_emu3_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL emu trigger 3 : "<<dl_emu3_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL emu trigger 4: "<<n_emu4<<"\n";
    double dl_emu4_lumi_sum = 0;
    for(int i=0; i<n_emu4; i++){
        dl_emu4_lumi_sum = dl_emu4_lumi_sum + dl_emu4_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL emu trigger 4 : "<<dl_emu4_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL mumu trigger 1: "<<n_mumu1<<"\n";
    double dl_mumu1_lumi_sum = 0;
    for(int i=0; i<n_mumu1; i++){
        dl_mumu1_lumi_sum = dl_mumu1_lumi_sum + dl_mumu1_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL mumu trigger 1 : "<<dl_mumu1_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL mumu trigger 2: "<<n_mumu2<<"\n";
    double dl_mumu2_lumi_sum = 0;
    for(int i=0; i<n_mumu2; i++){
        dl_mumu2_lumi_sum = dl_mumu2_lumi_sum + dl_mumu2_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL mumu trigger 2 : "<<dl_mumu2_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL mumu trigger 3: "<<n_mumu3<<"\n";
    double dl_mumu3_lumi_sum = 0;
    for(int i=0; i<n_mumu3; i++){
        dl_mumu3_lumi_sum = dl_mumu3_lumi_sum + dl_mumu3_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL mumu trigger 3 : "<<dl_mumu3_lumi_sum/1000<<"\n\n";

    std::cout<<"No. of fills for DL mumu trigger 4: "<<n_mumu4<<"\n";
    double dl_mumu4_lumi_sum = 0;
    for(int i=0; i<n_mumu4; i++){
        dl_mumu4_lumi_sum = dl_mumu4_lumi_sum + dl_mumu4_lumi[i];
    }
    std::cout<<"Integrated Luminosity for DL mumu trigger 4 : "<<dl_mumu4_lumi_sum/1000<<"\n\n";
    */

    // Single Muon
    double sl_mu_lumi_sum = 0;
    for(int i=0; i<n_mu1; i++){
        sl_mu_lumi_sum = sl_mu_lumi_sum + sl_mu1_lumi[i];
    }
    for(int i=0; i<n_mu2; i++){
        bool match = 0;
        for(int j=0; j<n_mu1; j++){
            if(sl_mu2_run[i] == sl_mu1_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            sl_mu_lumi_sum = sl_mu_lumi_sum + sl_mu2_lumi[i];
            std::cout<<sl_mu2_run[i]<<"\n";
        }
    }
    std::cout<<"Integrated Luminosity for SL mu channel : "<<sl_mu_lumi_sum/1000<<"\n\n";

    // Double Electron
    double dl_ee_lumi_sum = 0;
    for(int i=0; i<n_e; i++){
        dl_ee_lumi_sum = dl_ee_lumi_sum + sl_e_lumi[i];
    }
    for(int i=0; i<n_ee; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(dl_ee_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_ee_lumi_sum = dl_ee_lumi_sum + sl_e_lumi[i];
            std::cout<<dl_ee_run[i]<<"\n";
        }
    }
    std::cout<<"Integrated Luminosity for DL ee channel : "<<dl_ee_lumi_sum/1000<<"\n\n";

    // Electron Muon
    double dl_emu_lumi_sum = 0;
    for(int i=0; i<n_e; i++){
        dl_emu_lumi_sum = dl_emu_lumi_sum + sl_e_lumi[i];
    }
    for(int i=0; i<n_mu1; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(sl_mu1_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_emu_lumi_sum = dl_emu_lumi_sum + sl_mu1_lumi[i];
            std::cout<<sl_mu1_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_mu2; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(sl_mu2_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_emu_lumi_sum = dl_emu_lumi_sum + sl_mu2_lumi[i];
            std::cout<<sl_mu2_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_emu1; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(dl_emu1_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_emu_lumi_sum = dl_emu_lumi_sum + dl_emu1_lumi[i];
            std::cout<<dl_emu1_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_emu2; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(dl_emu2_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_emu_lumi_sum = dl_emu_lumi_sum + dl_emu2_lumi[i];
            std::cout<<dl_emu2_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_emu3; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(dl_emu3_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_emu_lumi_sum = dl_emu_lumi_sum + dl_emu3_lumi[i];
            std::cout<<dl_emu3_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_emu4; i++){
        bool match = 0;
        for(int j=0; j<n_e; j++){
            if(dl_emu4_run[i] == sl_e_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_emu_lumi_sum = dl_emu_lumi_sum + dl_emu4_lumi[i];
            std::cout<<dl_emu4_run[i]<<"\n";
        }
    }
    std::cout<<"Integrated Luminosity for DL emu channel : "<<dl_emu_lumi_sum/1000<<"\n\n";


    // Double Muon
    double dl_mumu_lumi_sum = 0;
    for(int i=0; i<n_mu1; i++){
        dl_mumu_lumi_sum = dl_mumu_lumi_sum + sl_mu1_lumi[i];
    }
    for(int i=0; i<n_mu2; i++){
        bool match = 0;
        for(int j=0; j<n_mu1; j++){
            if(sl_mu2_run[i] == sl_mu1_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_mumu_lumi_sum = dl_mumu_lumi_sum + sl_mu2_lumi[i];
            std::cout<<sl_mu2_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_mumu1; i++){
        bool match = 0;
        for(int j=0; j<n_mu1; j++){
            if(dl_mumu1_run[i] == sl_mu1_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_mumu_lumi_sum = dl_mumu_lumi_sum + dl_mumu1_lumi[i];
            std::cout<<dl_mumu1_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_mumu2; i++){
        bool match = 0;
        for(int j=0; j<n_mu1; j++){
            if(dl_mumu2_run[i] == sl_mu1_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_mumu_lumi_sum = dl_mumu_lumi_sum + dl_mumu2_lumi[i];
            std::cout<<dl_mumu2_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_mumu3; i++){
        bool match = 0;
        for(int j=0; j<n_mu1; j++){
            if(dl_mumu3_run[i] == sl_mu1_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_mumu_lumi_sum = dl_mumu_lumi_sum + dl_mumu3_lumi[i];
            std::cout<<dl_mumu3_run[i]<<"\n";
        }
    }
    for(int i=0; i<n_mumu4; i++){
        bool match = 0;
        for(int j=0; j<n_mu1; j++){
            if(dl_mumu4_run[i] == sl_mu1_run[j]){
                match = 1;
                break;
            }
        }
        if(match==0){
            dl_mumu_lumi_sum = dl_mumu_lumi_sum + dl_mumu4_lumi[i];
            std::cout<<dl_mumu4_run[i]<<"\n";
        }
    }

    std::cout<<"Integrated Luminosity for DL mumu channel : "<<dl_mumu_lumi_sum/1000<<"\n\n";
    









	return 0;
}
