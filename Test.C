const int NFILE=30;
const int NSIMU=5;
const int NHIST=1;
const int NPAD=9;
const int NPMT=120;

const double r_av=1.6; //m
const double h_av=1.2; //m
const double r_bb=1.86; //m
const double h_bb=1.5; //m

const int NPAR=3;
const int MaxMP=40;

TChain *chP = new TChain("DC");
TH1D *hPC[NHIST];
TH1D *hDC[NHIST];
TH1D *hDT[NHIST];
TH1D *hMP[NHIST];
TH1D *hDV[NHIST];
TH1D *hPW[NHIST]; //Prompt dWall
TH1D *hDW[NHIST]; //Delayed dWall
TH1D *hPA[NHIST]; //Prompt Angle
TH1D *hPM[NHIST]; //Prompt MSG

//TH2D *hVTX[NHIST];
//TH2D *hXYP[NHIST];
TH2D *hZRP[NHIST];
//TH2D *hXYD[NHIST];
TH2D *hZRD[NHIST];

double EPMIN=3.5;
double EPMAX=10;
double EDMIN=3.5;
double EDMAX=10;
double DTMIN=0;
double DTMAX=500;
double DVMIN=0; //150; //100; //150; //
double DVMAX=400; //150; //100; //150; //
double MPMIN=0;
double MPMAX=20;
double DWMIN=250;
double DWMAX=2000;

void makehist(int ich, int phase, int ngen){
    int MP;
    double DT[MaxMP], recoX[MaxMP][NPAR], recoR[MaxMP], dVTX[MaxMP];
    double BSEnergy[MaxMP], BSGood[MaxMP], dirKS[MaxMP];
    double dWall[MaxMP], angle[MaxMP], msg[MaxMP];
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    double n50[MaxMP], q50[MaxMP], Ovaq[MaxMP];
    double piLike[MaxMP], effWall[MaxMP];
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TVector3 v_p, v_d;

    chP->Reset();
    for(int isimu=0;isimu<NSIMU;isimu++){
    for(int ifile=0;ifile<NFILE;ifile++){
       // chP->Add(Form("../output/cosmicFastn_skmusic_SK%d_%d_Op01/dc/dc_%05d_0_900.root", phase, ngen, ifile));
        chP->Add(Form("../output/exp000%01d/dc_%05d_0_900.root", isimu, ifile));
    }
    }
    if(chP->GetEntries()<1){
        cout << "error" << endl;
        exit(0);
    }
    chP->SetBranchStatus("*", 0);
    chP->SetBranchStatus("MP", 1); chP->SetBranchAddress("MP", &MP);
    chP->SetBranchStatus("DT", 1); chP->SetBranchAddress("DT", DT);
    chP->SetBranchStatus("recoX", 1); chP->SetBranchAddress("recoX", recoX);
    chP->SetBranchStatus("recoR", 1); chP->SetBranchAddress("recoR", recoR);
    chP->SetBranchStatus("dVTX", 1); chP->SetBranchAddress("dVTX", dVTX);
    chP->SetBranchStatus("BSEnergy", 1); chP->SetBranchAddress("BSEnergy", BSEnergy);
    chP->SetBranchStatus("dWall", 1); chP->SetBranchAddress("dWall", dWall);
    chP->SetBranchStatus("angle", 1); chP->SetBranchAddress("angle", angle);
    chP->SetBranchStatus("msg", 1); chP->SetBranchAddress("msg", msg);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chP->SetBranchStatus("piLike", 1); chP->SetBranchAddress("piLike", piLike);
    chP->SetBranchStatus("n50", 1); chP->SetBranchAddress("n50", n50);
    chP->SetBranchStatus("q50", 1); chP->SetBranchAddress("q50", q50);
    chP->SetBranchStatus("Ovaq", 1); chP->SetBranchAddress("Ovaq", Ovaq);
    chP->SetBranchStatus("effWall", 1); chP->SetBranchAddress("effWall", effWall);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    double eEP=20;
    hPC[ich] = new TH1D(Form("hPC%d",ich), "Prompt ; Energy /MeV; Rate /kHz/0.2MeV", 100, 0, eEP);
    double eED=20;
    hDC[ich] = new TH1D(Form("hDC%d",ich), "Delayed ; Energy /MeV; Rate /kHz/0.2MeV", 100, 0, eED);

    double eDT=1000;
    hDT[ich] = new TH1D(Form("hDT%d",ich), "; #Delta t /#mus; Rate /kHz/10us", 100, 0, eDT);
    hMP[ich] = new TH1D(Form("hMP%d",ich), "; Multiplicity; event/bin", 20, 0, 20);
    hDV[ich] = new TH1D(Form("hDV%d",ich), "; #DeltaVTX /m; Rate /kHz/0.2m", 100, 0, 20);
    hPW[ich] = new TH1D(Form("hPW%d",ich), "; Prompt dWall /m; Rate /kHz", 100, 0, 20);
    hDW[ich] = new TH1D(Form("hDW%d",ich), "; Delayed dWall /m; Rate /kHz", 100, 0, 20);
    hPA[ich] = new TH1D(Form("hPA%d",ich), "; Prompt Angle /deg; Rate /kHz", 30, 0, 90);
    hPM[ich] = new TH1D(Form("hPM%d",ich), "; Prompt MSG; Rate /kHz", 25, 0, 1);
    //hVTX[ich] = new TH2D(Form("hVTX%d",ich), "; R^{2} /m^{2}; Z /m", 100, 0, 4, 100, -2, 2);
    //hXYP[ich] = new TH2D(Form("hXYP%d",ich), "Prompt Candidates; X /m; Y /m", 100, -2, 2, 100, -2, 2);
    hZRP[ich] = new TH2D(Form("hZRP%d",ich), "Prompt Candidates; R^{2} /m^{2}; Z /m", 100, 0, 100, 100, -20, 20);
    hZRD[ich] = new TH2D(Form("hZRD%d",ich), "Delayed Candidates; R^{2} /m^{2}; Z /m", 100, 0, 100, 100, -20, 20);

    //event loop
    int nevent = chP->GetEntries();
    cout << "nevent: " << nevent << endl;
    int countP=0, countD=0;
    vector<int> pFlag;
    for(int i=0;i<nevent;i++){
        chP->GetEntry(i);
        if(i%10000==0) cout << i << endl;
        //int promptFlag=0;
        pFlag.clear();

        //if(DT[0] > 100){ //signal isolation
            //delayed E and DT
            for(int imp=1;imp<MP;imp++){
                v_d.SetXYZ(recoX[imp][0], recoX[imp][1], recoX[imp][2]);

                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(BSEnergy[0]>EPMIN && BSEnergy[0]<EPMAX && dWall[0]>DWMIN

                && ((angle[0] > 38.0 && BSEnergy[0] >= 8.0)
                || (angle[0] > 32.0 && BSEnergy[0] < 8.0)) && angle[0] < 53.0

                && ((BSEnergy[0] < 7 && piLike[0] < 0.47)
              || (BSEnergy[0] >= 7.0 && piLike[0] < 0.37))

                && ((BSEnergy[0] < 4.0 && effWall[0] > 1360.0)
                || (BSEnergy[0] >= 4.0 && BSEnergy[0] < 5.0 && effWall[0] > 950.0)
                || (BSEnergy[0] >= 5.0 && BSEnergy[0] < 6.0 && effWall[0] > 610.0)
                || (BSEnergy[0] >= 6.0 && effWall[0] > 450.0))

                && (q50[0]/n50[0] < 2.0)

                && ((Ovaq[0] > 0.24 && BSEnergy[0] < 5.0)
                || (Ovaq[0] > 0.22 && BSEnergy[0] >= 5.0 && BSEnergy[0] < 6.0)
                || (Ovaq[0] > 0.18 && BSEnergy[0] >= 6.0 && BSEnergy[0] < 8.0)
                || (Ovaq[0] > 0.2 && BSEnergy[0] >= 8.0 && BSEnergy[0] < 10.0)
                || (BSEnergy[0] >= 10.0 && BSEnergy[0] < 12.0 && Ovaq[0] > 0.31)
                || (BSEnergy[0] >= 12.0 && Ovaq[0] > 0.36))

                ){ //prompt E
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

                        if(dVTX[imp]>0 && dVTX[imp]<DVMAX && dWall[imp]>DWMIN){ //dVTX
                        if(BSEnergy[imp]>EDMIN && BSEnergy[imp]<EDMAX){ //delayed E
                            hDT[ich]->Fill(DT[imp]);
                        }
                        if(DT[imp]>DTMIN && DT[imp]<DTMAX){ //dT
                            hPC[ich]->Fill(BSEnergy[imp]);
                           // hPC[ich]->Fill(BSEnergy[0]);
                            hDC[ich]->Fill(BSEnergy[imp]);
                            hDW[ich]->Fill(dWall[imp]*1e-2); //m
                        }
                    }
                }
            }
            //dVTX and promptFlag
            for(int imp=1;imp<MP;imp++){
                v_d.SetXYZ(recoX[imp][0], recoX[imp][1], recoX[imp][2]);
                if(BSEnergy[imp]>EDMIN && BSEnergy[imp]<EDMAX && dWall[imp]>DWMIN){
                    if(DT[imp]>DTMIN && DT[imp]<DTMAX) {

                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(BSEnergy[0]>EPMIN && BSEnergy[0]<EPMAX && dWall[0]>DWMIN

                && ((angle[0] > 38.0 && BSEnergy[0] >= 8.0)
                || (angle[0] > 32.0 && BSEnergy[0] < 8.0)) && angle[0] < 53.0

                && ((BSEnergy[0] < 7 && piLike[0] < 0.47)
                || (BSEnergy[0] >= 7.0 && piLike[0] < 0.37))

                && ((BSEnergy[0] < 4.0 && effWall[0] > 1360.0)
                || (BSEnergy[0] >= 4.0 && BSEnergy[0] < 5.0 && effWall[0] > 950.0)
                || (BSEnergy[0] >= 5.0 && BSEnergy[0] < 6.0 && effWall[0] > 610.0)
                || (BSEnergy[0] >= 6.0 && effWall[0] > 450.0))

                && (q50[0]/n50[0] < 2.0)

                && ((Ovaq[0] > 0.24 && BSEnergy[0] < 5.0)
                || (Ovaq[0] > 0.22 && BSEnergy[0] >= 5.0 && BSEnergy[0] < 6.0)
                || (Ovaq[0] > 0.18 && BSEnergy[0] >= 6.0 && BSEnergy[0] < 8.0)
                || (Ovaq[0] > 0.2 && BSEnergy[0] >= 8.0 && BSEnergy[0] < 10.0)
                || (BSEnergy[0] >= 10.0 && BSEnergy[0] < 12.0 && Ovaq[0] > 0.31)
                || (BSEnergy[0] >= 12.0 && Ovaq[0] > 0.36))

                ){ //prompt E
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                // && TotalPEV[0]<VPMIN  && MP==2
                        //hDV[ich]->Fill(dVTX[imp]); //cm
                        hDV[ich]->Fill(dVTX[imp]*0.01); //m

            if(dVTX[imp]<DVMAX) pFlag.push_back(imp); //dVTX[imp]>0&&
                if(pFlag.size()>0){
                    if(dWall[0]>DWMIN) hPC[ich]->Fill(BSEnergy[0]);
                        hPW[ich]->Fill(dWall[0]*1e-2); //m
                        hPA[ich]->Fill(angle[0]); //m
                        hPM[ich]->Fill(msg[0]); //m
                        hZRP[ich]->Fill(pow(recoR[0]/100.,2), recoX[0][2]/100.); //m                   
                        if(dWall[0]>DWMIN){
                                hMP[ich]->Fill(pFlag.size());
                                countP++;
                        }
                    //hVP[ich]->Fill(TotalPEV[0]);
                    //if(ich==0) cout << "Prompt: " << TrigID[0] << endl;
                }
                //for(const auto &e: pFlag){
                        for(int idx=0;idx<pFlag.size();idx++){
                                int e = pFlag.at(idx);
                                if(BSEnergy[0]>EPMIN && BSEnergy[0]<EPMAX){ // && MP==2
                                        countD++;
                                        //hXYD[ich]->Fill(recoX[e][1]/100., recoX[e][2]/100.); //m
                                        hZRD[ich]->Fill(pow(recoR[e]/100.,2), recoX[e][2]/100.); //cm
                                        //if(ich==0) cout << "Delayed: " << TrigID[e] << endl;          
                                }
                        }
                }
            }
                        }
                }
        }
    //}
    cout << "Prompt: " << endl;
    cout << countP << Form(" in %.1f - %.1f MeV", EPMIN, EPMAX) << endl;
    //cout << countP/(double)entry[ich]*SingRate[ich] << " kHz" << endl;
    //cout << "Delayed: " << endl;
    //cout << countD << Form(" in %.1f - %.1f MeV", EDMIN, EDMAX) << endl;
    //cout << countD/(double)entry[ich]*SingRate[ich] << " kHz" << endl;
}

void drawDC(int ngen = 10000){
    int phase = 7;
    int range[2][2]={{0, 1000}, {2000, 3000}};
    //int range[2][2]={{0, 2000}, {4000, 6000}};
    //int range[2][2]={{0, 1000}, {1000, 2000}};

    for(int ich=0;ich<NHIST;ich++){
        makehist(ich, phase, ngen);
    }
    TH1D *hPCE, *hDCE, *hDTE, *hMPE, *hDVE, *hDWE; //excess
    TH1D *hPWE, *hPAE, *hPME; //excess

    int selecFlag[NPAD]={1, 1, 1, 1, 0, 1, 1, 0, 0};
    double selection[NPAD][2]={
        {EPMIN, EPMAX},
        {EDMIN, EDMAX},
        {DTMIN, DTMAX},
        {DVMIN*0.01, DVMAX*0.01},
        {MPMIN, MPMAX},
        {DWMIN*1e-2, DWMAX*1e-2},
        {DWMIN*1e-2, DWMAX*1e-2},
        {DTMIN, DTMAX},
        {DTMIN, DTMAX}
    };
    double ymin[NPAD], ymax[NPAD];

    int imode = 1;
    double cx[]={900, 1400}, cy[]={1000, 900};
    TCanvas *c1 = new TCanvas("c1", "", cx[imode], cy[imode]);
        c1->Divide(3, 3);
    for(int n=0;n<NPAD;n++){
        c1->cd(n+1);
        gPad->SetLogy();
        gPad->SetGrid();
        if(n==0){
            for(int ich=0;ich<NHIST;ich++){
                hPC[ich]->SetStats(0);
                hPC[ich]->Sumw2();
                //hPC[ich]->Scale(SingRate[id]/entry[id]);

                hPC[ich]->SetLineWidth(2); // Adjust the value to make the line thicker

                hPC[ich]->GetYaxis()->SetTitleFont(132);
                double xmax=hPC[ich]->GetXaxis()->GetXmax();
                double bins=(double)hPC[ich]->GetXaxis()->GetNbins();
                //hPC[ich]->SetTitle(Form("; Prompt Evis /MeV; Rate /kHz /%1.2f MeV", xmax/bins));
                hPC[ich]->GetXaxis()->SetTitleSize(0.06);
                hPC[ich]->GetYaxis()->SetTitleSize(0.06);
                hPC[ich]->GetXaxis()->SetTitleOffset(0.7);
                hPC[ich]->GetYaxis()->SetTitleOffset(0.7);


                hPC[ich]->SetMarkerStyle(8);
                hPC[ich]->SetMarkerColor(ich+3);
                hPC[ich]->SetLineColor(ich+3);
                if(ich==0) hPC[ich]->Draw("pe");
                else hPC[ich]->Draw("pe same");

        // Calculer et afficher l'intégrale de chaque histograme
        double integral = hPC[ich]->Integral();
        std::cout << "L'intégrale de l'histogramme hPC[" << ich << "] est : " <<integral << std::endl;


            }
            if(NHIST>1){
                hPCE = (TH1D*)hPC[0]->Clone("hPCE");
                hPCE->SetStats(0);
                hPCE->Add(hPC[1], -1);
                hPCE->SetMarkerStyle(8);
                hPCE->SetMarkerColor(4);
                hPCE->SetLineColor(4);
                hPCE->Draw("pe same");
            }

            ymin[n] = hPC[0]->GetMinimum();
            ymax[n] = hPC[0]->GetMaximum();
        }
        if(n==1){
            for(int ich=0;ich<NHIST;ich++){
                hDC[ich]->SetStats(0);
                hDC[ich]->Sumw2();
                //hDC[ich]->Scale(SingRate[id]/entry[id]);
                hDC[ich]->SetLineWidth(2); // Adjust the value to make the line thicker
                hDC[ich]->GetYaxis()->SetTitleFont(132);
                double xmax=hDC[ich]->GetXaxis()->GetXmax();
                double bins=(double)hDC[ich]->GetXaxis()->GetNbins();
                //hDC[ich]->SetTitle(Form("; Prompt Evis /MeV; Rate /kHz /%1.2f MeV", xmax/bins));
                //hDC[ich]->SetTitle(Form("; Delayed Evis /MeV; Rate /kHz /%1.2f MeV", xmax/bins));
                hDC[ich]->GetXaxis()->SetTitleSize(0.06);
                hDC[ich]->GetYaxis()->SetTitleSize(0.06);
                hDC[ich]->GetXaxis()->SetTitleOffset(0.7);
                hDC[ich]->GetYaxis()->SetTitleOffset(0.7);

                hDC[ich]->SetMarkerStyle(8);
                hDC[ich]->SetMarkerColor(ich+3);
                hDC[ich]->SetLineColor(ich+3);
                if(ich==0) hDC[ich]->Draw("pe");
                else hDC[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hDCE = (TH1D*)hDC[0]->Clone("hDCE");
                hDCE->SetStats(0);
                hDCE->Add(hDC[1], -1);
                hDCE->SetMarkerStyle(8);
                hDCE->SetMarkerColor(4);
                hDCE->SetLineColor(4);
                hDCE->Draw("pe same");
            }

            ymin[n] = hDC[0]->GetMinimum();
            ymax[n] = hDC[0]->GetMaximum();
        }
        if(n==2){
            for(int ich=0;ich<NHIST;ich++){
                hDT[ich]->SetStats(0);
                hDT[ich]->Sumw2();
                //hDT[ich]->Scale(SingRate[id]/entry[id]);
                hDT[ich]->SetLineWidth(2); // Adjust the value to make the line thicker
                hDT[ich]->GetYaxis()->SetTitleFont(132);
                double xmax=hDT[ich]->GetXaxis()->GetXmax();
                double bins=(double)hDT[ich]->GetXaxis()->GetNbins();
                //hDT[ich]->SetTitle(Form("; Prompt Evis /MeV; Rate /kHz /%1.2f MeV", xmax/bins));
                //hDT[ich]->SetTitle(Form("; #Deltat /#mus; Rate /kHz /%1.2f #mus", xmax/bins));
                hDT[ich]->GetXaxis()->SetTitleSize(0.06);
                hDT[ich]->GetYaxis()->SetTitleSize(0.06);
                hDT[ich]->GetXaxis()->SetTitleOffset(0.7);
                hDT[ich]->GetYaxis()->SetTitleOffset(0.7);

                hDT[ich]->SetMarkerStyle(8);
                hDT[ich]->SetMarkerColor(ich+3);
                hDT[ich]->SetLineColor(ich+3);
                if(ich==0) hDT[ich]->Draw("pe");
                else hDT[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hDTE = (TH1D*)hDT[0]->Clone("hDTE");
                hDTE->SetStats(0);
                hDTE->Add(hDT[1], -1);
                hDTE->SetMarkerStyle(8);
                hDTE->SetMarkerColor(4);
                hDTE->SetLineColor(4);
                hDTE->Draw("pe same");
            }

            ymin[n] = hDT[0]->GetMinimum();
            ymax[n] = hDT[0]->GetMaximum();
        }
        if(n==3){
            for(int ich=0;ich<NHIST;ich++){
                hDV[ich]->SetStats(0);
                hDV[ich]->Sumw2();
                //hDV[ich]->Scale(SingRate[id]/entry[id]);
                hDV[ich]->SetLineWidth(2); // Adjust the value to make the line thicker
                double xmax=hDV[ich]->GetXaxis()->GetXmax();
                double bins=(double)hDV[ich]->GetXaxis()->GetNbins();
                //hDV[ich]->SetTitle(Form("; #DeltaVTX /cm; Rate /kHz /%1.2f cm", xmax/bins));
                hDV[ich]->GetXaxis()->SetTitleSize(0.06);
                hDV[ich]->GetYaxis()->SetTitleSize(0.06);
                hDV[ich]->GetXaxis()->SetTitleOffset(0.7);
                hDV[ich]->GetYaxis()->SetTitleOffset(0.7);

                hDV[ich]->SetMarkerStyle(8);
                hDV[ich]->SetMarkerColor(ich+3);
                hDV[ich]->SetLineColor(ich+3);
                if(ich==0) hDV[ich]->Draw("pe");
                else hDV[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hDVE = (TH1D*)hDV[0]->Clone("hDVE");
                hDVE->SetStats(0);
                hDVE->Add(hDV[1], -1);
                hDVE->SetMarkerStyle(8);
                hDVE->SetMarkerColor(4);
                hDVE->SetLineColor(4);
                hDVE->Draw("pe same");
            }

            ymin[n] = 0; //hDV[0]->GetMinimum();
            ymax[n] = hDV[0]->GetMaximum();
        }
        if(n==4){
            for(int ich=0;ich<NHIST;ich++){
                hMP[ich]->SetStats(0);
                hMP[ich]->GetYaxis()->SetTitleFont(132);
                double xmax=hMP[ich]->GetXaxis()->GetXmax();
                double bins=(double)hMP[ich]->GetXaxis()->GetNbins();
                //hMP[ich]->SetTitle(Form("; Prompt Evis /MeV; Rate /kHz /%1.2f MeV", xmax/bins));
   hMP[ich]->GetXaxis()->SetTitleSize(0.06);
                hMP[ich]->GetYaxis()->SetTitleSize(0.06);
                hMP[ich]->GetXaxis()->SetTitleOffset(0.7);
                hMP[ich]->GetYaxis()->SetTitleOffset(0.7);
                hMP[ich]->SetLineWidth(2); // Adjust the value to make the line thicker

                hMP[ich]->SetMarkerStyle(8);
                hMP[ich]->SetMarkerColor(ich+3);
                hMP[ich]->SetLineColor(ich+3);
                if(ich==0) hMP[ich]->Draw("pe");
                else hMP[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hMPE = (TH1D*)hMP[0]->Clone("hMPE");
                hMPE->SetStats(0);
                hMPE->Add(hMP[1], -1);
                hMPE->SetMarkerStyle(8);
                hMPE->SetMarkerColor(4);
                hMPE->SetLineColor(4);
                hMPE->Draw("pe same");
            }

            ymin[n] = hMP[0]->GetMinimum();
            ymax[n] = hMP[0]->GetMaximum();
        }
        if(n==5){
            for(int ich=0;ich<NHIST;ich++){
                hDW[ich]->SetStats(0);
                hDW[ich]->Sumw2();
                hDW[ich]->SetLineWidth(2); // Adjust the value to make the line thicker
                //hDW[ich]->Scale(SingRate[id]/entry[id]);
                double xmax=hDW[ich]->GetXaxis()->GetXmax();
                double bins=(double)hDW[ich]->GetXaxis()->GetNbins();
                //hDW[ich]->SetTitle(Form("; Veto Charge /p.e.; Rate /kHz /%2.1f p.e.", xmax/bins));
                hDW[ich]->GetXaxis()->SetTitleSize(0.06);
                hDW[ich]->GetYaxis()->SetTitleSize(0.06);
                hDW[ich]->GetXaxis()->SetTitleOffset(0.7);
                hDW[ich]->GetYaxis()->SetTitleOffset(0.7);

                hDW[ich]->SetMarkerStyle(8);
                hDW[ich]->SetMarkerColor(ich+3);
                hDW[ich]->SetLineColor(ich+3);
                if(ich==0) hDW[ich]->Draw("pe");
                else hDW[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hDWE = (TH1D*)hDW[0]->Clone("hDWE");
                hDWE->SetStats(0);
                hDWE->Add(hDW[1], -1);
                hDWE->SetMarkerStyle(8);
                hDWE->SetMarkerColor(4);
                hDWE->SetLineColor(4);
                hDWE->Draw("pe same");
                //hVTX[0]->Draw("colz");
            }

            ymin[n] = hDW[0]->GetMinimum();
            ymax[n] = hDW[0]->GetMaximum();
        }
        if(n==6){
            for(int ich=0;ich<NHIST;ich++){
                hPW[ich]->SetStats(0);
                hPW[ich]->Sumw2();
                //hPW[ich]->Scale(SingRate[id]/entry[id]);
                hPW[ich]->SetLineWidth(2); // Adjust the value to make the line thicker
                double xmax=hPW[ich]->GetXaxis()->GetXmax();
                double bins=(double)hPW[ich]->GetXaxis()->GetNbins();
                //hPW[ich]->SetTitle(Form("; Veto Charge /p.e.; Rate /kHz /%2.1f p.e.", xmax/bins));
                hPW[ich]->GetXaxis()->SetTitleSize(0.06);
               hPW[ich]->GetYaxis()->SetTitleSize(0.06);
                hPW[ich]->GetXaxis()->SetTitleOffset(0.7);
                hPW[ich]->GetYaxis()->SetTitleOffset(0.7);

                hPW[ich]->SetMarkerStyle(8);
                hPW[ich]->SetMarkerColor(ich+3);
                hPW[ich]->SetLineColor(ich+3);
                if(ich==0) hPW[ich]->Draw("pe");
                else hPW[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hPWE = (TH1D*)hPW[0]->Clone("hPWE");
                hPWE->SetStats(0);
                hPWE->Add(hPW[1], -1);
                hPWE->SetMarkerStyle(8);
                hPWE->SetMarkerColor(4);
                hPWE->SetLineColor(4);
                hPWE->Draw("pe same");
                //hVTX[0]->Draw("colz");
            }

            ymin[n] = hPW[0]->GetMinimum();
            ymax[n] = hPW[0]->GetMaximum();
        }
        if(n==7){
            for(int ich=0;ich<NHIST;ich++){
                hPA[ich]->SetStats(0);
                hPA[ich]->Sumw2();
                //hPA[ich]->Scale(SingRate[id]/entry[id]);
                hPA[ich]->SetLineWidth(2); // Adjust the value to make the line thicker
                double xmax=hPA[ich]->GetXaxis()->GetXmax();
                double bins=(double)hPA[ich]->GetXaxis()->GetNbins();
                //hPA[ich]->SetTitle(Form("; Veto Charge /p.e.; Rate /kHz /%2.1f p.e.", xmax/bins));
                hPA[ich]->GetXaxis()->SetTitleSize(0.06);
                hPA[ich]->GetYaxis()->SetTitleSize(0.06);
                hPA[ich]->GetXaxis()->SetTitleOffset(0.7);
                hPA[ich]->GetYaxis()->SetTitleOffset(0.7);

                hPA[ich]->SetMarkerStyle(8);
                hPA[ich]->SetMarkerColor(ich+3);
                hPA[ich]->SetLineColor(ich+3);
                if(ich==0) hPA[ich]->Draw("pe");
                else hPA[ich]->Draw("pe same");
            }
            if(NHIST>1){
                hPAE = (TH1D*)hPA[0]->Clone("hPAE");
                hPAE->SetStats(0);
                                                                                                453,1         79%

                                                                                                                                                                                                 1,1           Top
